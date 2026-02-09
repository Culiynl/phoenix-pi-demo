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
from config import VINA_EXE, DEFAULT_RECEPTOR_PDBQT, DOCKING_CONFIG, DEVICE

# --- AI MODELS ---

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
        self.means = np.array([350.0, 3.0, 2.0, 5.0, 80.0])
        self.stds = np.array([100.0, 1.5, 1.5, 3.0, 40.0])
    def normalize(self, raw): return (raw - self.means) / (self.stds + 1e-8)

# --- UTILS ---

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
            if smi and Chem.MolFromSmiles(smi): candidates.append(smi)
        except: pass
    return list(set(candidates)), goal_numpy

# FIXED: Now accepts 'center' dictionary
def run_vina_standalone(smiles, center):
    if not smiles or not center: return smiles, 0.0, "Input Error"
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return smiles, 0.0, "Invalid SMILES"

    try:
        m = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) == -1: return smiles, 0.0, "Embed Fail"
        prep = MoleculePreparation()
        setup = prep.prepare(m)
        pdbqt_str = PDBQTWriterLegacy().write_string(setup[0])[0]

        with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False, mode='w') as tmp:
            tmp.write(pdbqt_str)
            in_p = tmp.name
        out_p = in_p.replace(".pdbqt", "_out.pdbqt")

        size = DOCKING_CONFIG["SIZE"]
        cmd = [
            VINA_EXE, "--receptor", DEFAULT_RECEPTOR_PDBQT, "--ligand", in_p,
            "--center_x", str(center["x"]), "--center_y", str(center["y"]), "--center_z", str(center["z"]),
            "--size_x", size["x"], "--size_y", size["y"], "--size_z", size["z"],
            "--out", out_p, "--exhaustiveness", "1", "--cpu", "1"
        ]
        
        subprocess.run(cmd, capture_output=True, timeout=30)
        score = 0.0
        if os.path.exists(out_p):
            with open(out_p, 'r') as f:
                content = f.read()
                match = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", content)
                if match: score = float(match.group(1))
        
        for p in [in_p, out_p]:
            if os.path.exists(p): os.remove(p)
        return smiles, score, "Success"
    except Exception as e:
        return smiles, 0.0, str(e)