import os
import subprocess
import json
import pandas as pd
import numpy as np
import re
from Bio.PDB import PDBList, PDBParser, Select, PDBIO

# Safe Import for PyMol
try:
    import pymol
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False

# Safe Import for ChEMBL
try:
    from chembl_webresource_client.new_client import new_client
    CHEMBL_AVAILABLE = True
except Exception:
    CHEMBL_AVAILABLE = False

from config import PDB_DIR, P2RANK_EXEC, IMG_DIR, AIZYNTH_PYTHON_EXE, AIZYNTH_RUNNER_SCRIPT, AIZYNTH_CONFIG, BASE_DIR
from state import log_msg

class ProteinManager:
    @staticmethod
    def download_and_clean_pdb(pdb_id: str):
        log_msg(f"   [PDB] Downloading {pdb_id}...")
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(pdb_id, pdir=PDB_DIR, file_format='pdb')
        
        expected_path = os.path.join(PDB_DIR, f"pdb{pdb_id.lower()}.ent")
        if not os.path.exists(expected_path):
             expected_path = os.path.join(PDB_DIR, f"{pdb_id.lower()}.pdb")
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, expected_path)
        io = PDBIO()
        io.set_structure(structure)
        
        class NotWater(Select):
            def accept_residue(self, r): return r.get_resname() != "HOH"
            
        clean_path = os.path.join(PDB_DIR, f"{pdb_id}_clean.pdb")
        io.save(clean_path, NotWater())
        return clean_path

    @staticmethod
    def find_pocket_p2rank(pdb_path: str, pdb_id: str):
        if not os.path.exists(P2RANK_EXEC):
            log_msg("[P2Rank] Executable not found. Using Fallback.")
            return ProteinManager.find_pocket_center_native(pdb_path)

        out_dir = os.path.join(PDB_DIR, "p2rank_out")
        os.makedirs(out_dir, exist_ok=True)
        
        cmd = [P2RANK_EXEC, "predict", "-f", pdb_path, "-o", out_dir, "-threads", "1", "-visualizations", "0"]
        try:
            use_shell = (os.name == 'nt')
            subprocess.run(cmd, check=True, capture_output=True, text=True, shell=use_shell)
            
            base_name = os.path.basename(pdb_path)
            csv_name = f"{base_name}_predictions.csv"
            csv_path = os.path.join(out_dir, csv_name)
            
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                df.columns = df.columns.str.strip()
                if not df.empty:
                    row = df.iloc[0]
                    log_msg(f"   [P2Rank] Pocket found: Score {row['score']}")
                    return (float(row['center_x']), float(row['center_y']), float(row['center_z']))
        except Exception as e:
            log_msg(f"[P2Rank Error] {e}")

        return ProteinManager.find_pocket_center_native(pdb_path)

    @staticmethod
    def find_pocket_center_native(pdb_path: str):
        parser = PDBParser(QUIET=True)
        s = parser.get_structure("x", pdb_path)
        atoms = [a.get_vector() for a in s.get_atoms()]
        c = sum(atoms, start=np.array([0,0,0])) / len(atoms)
        return (c[0], c[1], c[2])

    @staticmethod
    def generate_pocket_viz(pdb_path: str, center: tuple, pdb_id: str):
        if not PYMOL_AVAILABLE: return None
        try:
            pymol.cmd.reinitialize()
            pymol.cmd.load(pdb_path, "prot")
            pymol.cmd.hide("all")
            pymol.cmd.show("cartoon", "prot")
            pymol.cmd.color("slate", "prot")
            pymol.cmd.pseudoatom("p", pos=[center[0], center[1], center[2]])
            pymol.cmd.show("spheres", "p")
            pymol.cmd.set("sphere_scale", 3.0, "p")
            pymol.cmd.color("red", "p")
            pymol.cmd.zoom("p", buffer=15)
            fname = f"{pdb_id}_viz.png"
            pymol.cmd.png(os.path.join(IMG_DIR, fname), width=600, height=400, ray=1)
            return f"http://localhost:8000/static/{fname}"
        except: return None

class RetroManager:
    @staticmethod
    def run_retrosynthesis(smiles: str):
        print(f"\n🚀 [BACKEND] STARTING RETROSYNTHESIS FOR: {smiles}", flush=True)
        
        if not os.path.exists(AIZYNTH_PYTHON_EXE): return {"error": "AiZynth Python Env Missing"}
        if not os.path.exists(AIZYNTH_RUNNER_SCRIPT): return {"error": "Runner Script Missing"}

        cmd = [AIZYNTH_PYTHON_EXE, AIZYNTH_RUNNER_SCRIPT, "--smiles", smiles, "--config", AIZYNTH_CONFIG]
        
        try:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=BASE_DIR)
            while True:
                line = process.stderr.readline()
                if not line and process.poll() is not None: break
                if line: print(line.strip(), flush=True)

            stdout, _ = process.communicate()
            if process.returncode != 0: return {"error": "Subprocess failed"}
            return json.loads(stdout)
        except Exception as e:
            return {"error": str(e)}

def fetch_chembl_data(target_name: str):
    if not CHEMBL_AVAILABLE: return None
    search_term = target_name
    t = new_client.target.filter(target_synonym__icontains=search_term).filter(target_type='SINGLE PROTEIN')
    
    if not t and "(" in search_term:
        match = re.search(r'\((.*?)\)', search_term)
        if match:
            t = new_client.target.filter(target_synonym__icontains=match.group(1)).filter(target_type='SINGLE PROTEIN')
            
    if not t: return None
    
    acts = new_client.activity.filter(target_chembl_id=t[0]['target_chembl_id'], standard_type="IC50", standard_value__isnull=False)
    data = []
    for a in acts:
        if len(data) > 200: break
        if a.get('canonical_smiles'):
            data.append({"smiles": a['canonical_smiles'], "standard_value": float(a['standard_value'])})
    return {"compounds": data}