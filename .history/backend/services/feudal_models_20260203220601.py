import torch
import torch.nn as nn
import torch.nn.functional as F
import selfies as sf
import numpy as np
import json
import os
import random
import tempfile
import subprocess
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

# Import paths from central config
from config import VINA_EXE, DEFAULT_RECEPTOR_PDBQT, DOCKING_CONFIG

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class Manager(nn.Module):
    def __init__(self, input_dim=5, goal_dim=5):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 128), nn.ReLU(),
            nn.Linear(128, 64), nn.ReLU(),
            nn.Linear(64, goal_dim)
        )
    def forward(self, x): return self.net(x)

class Worker(nn.Module):
    def __init__(self, vocab_size, d=256, cond_dim=6):
        super().__init__()
        self.emb = nn.Embedding(vocab_size, d)
        self.pos = nn.Parameter(torch.randn(1, 512, d))
        self.tf = nn.Transformer(d_model=d, nhead=8, num_encoder_layers=3, 
                                 num_decoder_layers=3, batch_first=True)
        self.film_net = nn.Sequential(
            nn.Linear(cond_dim, 128), nn.SiLU(),
            nn.Linear(128, 256), nn.SiLU(),
            nn.Linear(256, d * 2) 
        )
        self.policy_film = nn.Sequential(
            nn.Linear(cond_dim, 64), nn.SiLU(),
            nn.Linear(64, d * 2)
        )
        self.post_film_ln = nn.LayerNorm(d)
        self.policy = nn.Linear(d, vocab_size)

    def forward(self, src, cond, tgt):
        cond = torch.clamp(cond, -3.0, 3.0) 
        f_params = self.film_net(cond).unsqueeze(1)
        gamma, beta = torch.chunk(f_params, 2, dim=-1)
        s_emb = (self.emb(src) + self.pos[:, :src.size(1)]) * (1 + 0.5 * gamma) + 0.5 * beta
        t_emb = (self.emb(tgt) + self.pos[:, :tgt.size(1)]) * (1 + 0.5 * gamma) + 0.5 * beta
        s_emb = self.post_film_ln(s_emb)
        t_emb = self.post_film_ln(t_emb)
        mask = nn.Transformer.generate_square_subsequent_mask(tgt.size(1)).to(DEVICE)
        h = self.tf(s_emb, t_emb, tgt_mask=mask)
        p_f = self.policy_film(cond).unsqueeze(1)
        pg, pb = torch.chunk(p_f, 2, dim=-1)
        return self.policy(h * (1 + pg) + pb)

class Vocab:
    def __init__(self, path):
        with open(path, 'r') as f: self.all = json.load(f)
        self.c2i = {t: i for i, t in enumerate(self.all)}
        self.i2c = {i: t for i, t in enumerate(self.all)}
        self.n = len(self.all)

class DatasetNormalizer:
    def __init__(self):
        # Default stats from your training
        self.means = np.array([350.0, 3.0, 2.0, 5.0, 80.0])
        self.stds = np.array([100.0, 1.5, 1.5, 3.0, 40.0])

    def normalize(self, raw): return (raw - self.means) / (self.stds + 1e-8)

def get_props_raw(mol):
    if not mol: return np.zeros(5)
    return np.array([Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Descriptors.NumHDonors(mol), 
                     Descriptors.NumHAcceptors(mol), Descriptors.TPSA(mol)], dtype=np.float32)

def generate_candidates(start_smi, manager, worker, vocab, norm, n_samples=16):
    mol = Chem.MolFromSmiles(start_smi)
    if not mol: return [], [0.0]*5
    
    props = get_props_raw(mol)
    z = torch.tensor(norm.normalize(props), dtype=torch.float32, device=DEVICE).unsqueeze(0)
    
    with torch.no_grad():
        raw_goal = manager(z)
        goal = F.normalize(raw_goal, p=2, dim=-1) * 2.0 
    
    goal_numpy = raw_goal.cpu().numpy().tolist()[0]
    cond = torch.cat([goal, torch.tensor([[1.0]], device=DEVICE)], dim=-1).repeat(n_samples, 1)

    try:
        enc = sf.encoder(start_smi)
        tokens = [vocab.c2i.get(t, 0) for t in sf.split_selfies(enc)]
    except: return [], goal_numpy

    src = torch.tensor([tokens], device=DEVICE).repeat(n_samples, 1)
    tgt = torch.tensor([[vocab.c2i.get("<sos>", 0)]], device=DEVICE).repeat(n_samples, 1)

    for _ in range(70):
        with torch.no_grad():
            logits = worker(src, cond, tgt)
        probs = F.softmax(logits[:, -1, :], dim=-1)
        next_tok = torch.multinomial(probs, 1)
        tgt = torch.cat([tgt, next_tok], dim=1)
        if (next_tok == vocab.c2i.get("<eos>", 0)).all(): break

    candidates = []
    tgt_cpu = tgt.cpu().numpy()
    for i in range(n_samples):
        tok_idxs = [t for t in tgt_cpu[i] if t not in [vocab.c2i.get("<sos>",0), vocab.c2i.get("<eos>",0)]]
        tok_strs = [vocab.i2c.get(t, "") for t in tok_idxs]
        try:
            smi = sf.decoder("".join(tok_strs))
            if smi and Chem.MolFromSmiles(smi):
                candidates.append(smi)
        except: pass
        
    return list(set(candidates)), goal_numpy

def run_vina_standalone(smiles):
    print(f"\n[VINA DEBUG] Starting docking for: {smiles}")
    if not smiles: 
        print("[VINA DEBUG] Error: SMILES is empty.")
        return smiles, 0.0, "Empty"
        
    mol = Chem.MolFromSmiles(smiles)
    if not mol: 
        print("[VINA DEBUG] Error: RDKit could not parse SMILES.")
        return smiles, 0.0, "Invalid SMILES"

    try:
        # 1. Check Paths
        if not os.path.exists(VINA_EXE):
            print(f"[VINA DEBUG] CRITICAL: Vina EXE not found at: {VINA_EXE}")
            return smiles, 0.0, "Vina Path Error"
            
        if not os.path.exists(DEFAULT_RECEPTOR_PDBQT):
            print(f"[VINA DEBUG] CRITICAL: Receptor not found at: {DEFAULT_RECEPTOR_PDBQT}")
            return smiles, 0.0, "Receptor Path Error"

        # 2. Prepare Ligand
        print("[VINA DEBUG] Embedding molecule...")
        m = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) == -1:
            print("[VINA DEBUG] Embedding failed.")
            return smiles, 0.0, "Embed Fail"

        prep = MoleculePreparation()
        setup = prep.prepare(m)
        pdbqt_str = PDBQTWriterLegacy().write_string(setup[0])[0]

        # 3. Create Temp Files
        with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False, mode='w') as tmp:
            tmp.write(pdbqt_str)
            in_p = tmp.name
        out_p = in_p.replace(".pdbqt", "_out.pdbqt")
        print(f"[VINA DEBUG] Temp ligand created: {in_p}")

        # 4. Build Command
        center = DOCKING_CONFIG["CENTER"]
        size = DOCKING_CONFIG["SIZE"]

        cmd = [
            VINA_EXE, 
            "--receptor", DEFAULT_RECEPTOR_PDBQT, 
            "--ligand", in_p,
            "--center_x", str(center["x"]), 
            "--center_y", str(center["y"]), 
            "--center_z", str(center["z"]),
            "--size_x", size["x"], 
            "--size_y", size["y"], 
            "--size_z", size["z"],
            "--out", out_p, 
            "--exhaustiveness", "1", # Low for speed in RL
            "--cpu", "1"
        ]
        
        print(f"[VINA DEBUG] Executing: {' '.join(cmd)}")

        # 5. Run Vina
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            print(f"[VINA DEBUG] Vina Error (Code {result.returncode})")
            print(f"[VINA DEBUG] STDERR: {result.stderr}")
            return smiles, 0.0, "Vina Process Error"

        # 6. Parse Output
        score = 0.0
        if os.path.exists(out_p):
            with open(out_p, 'r') as f:
                content = f.read()
                match = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", content)
                if match:
                    score = float(match.group(1))
                    print(f"[VINA DEBUG] Success! Score: {score}")
                else:
                    print("[VINA DEBUG] Could not find score in output file.")
        else:
            print("[VINA DEBUG] Output file was not created by Vina.")
        
        # Cleanup
        for p in [in_p, out_p]:
            if os.path.exists(p): os.remove(p)

        return smiles, score, "Success"

    except subprocess.TimeoutExpired:
        print("[VINA DEBUG] Error: Vina timed out (took longer than 60s)")
        return smiles, 0.0, "Timeout"
    except Exception as e:
        print(f"[VINA DEBUG] Unexpected Exception: {str(e)}")
        return smiles, 0.0, str(e)