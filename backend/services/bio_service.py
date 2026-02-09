import os
import subprocess
import json
import tempfile
import pandas as pd
import numpy as np
import re
import sys
from Bio.PDB import PDBList, PDBParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import Descriptors, QED
from .sascorer import calculateScore

# --- SAFE IMPORTS ---
try:
    import pymol
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False

try:
    from chembl_webresource_client.new_client import new_client
    CHEMBL_AVAILABLE = True
except Exception:
    CHEMBL_AVAILABLE = False

# --- CONFIG & STATE ---
from config import (
    BASE_DIR, PDB_DIR, P2RANK_EXEC, IMG_DIR, 
    AIZYNTH_ENV_PYTHON, AIZYNTH_CLI, AIZYNTH_CONFIG
)
from state import log_msg


class ProteinManager:
    @staticmethod
    def download_and_clean_pdb(pdb_id: str, mission_id: str = None):
        """Downloads a PDB file, removes water, and saves it. If mission_id provided, saves to project dir."""
        from config import PROJECTS_DIR, PDB_DIR
        log_msg(f"   [PDB] Downloading {pdb_id}...")
        pdbl = PDBList()
        
        # Download to a temporary location first if we want to move it
        temp_dir = tempfile.gettempdir()
        pdbl.retrieve_pdb_file(pdb_id, pdir=temp_dir, file_format='pdb')
        
        expected_raw = os.path.join(temp_dir, f"pdb{pdb_id.lower()}.ent")
        if not os.path.exists(expected_raw):
             # Try other formats if .ent fails
             expected_raw = os.path.join(temp_dir, f"{pdb_id.lower()}.pdb")
        
        if not os.path.exists(expected_raw):
            raise Exception(f"Failed to find downloaded PDB for {pdb_id}")

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, expected_raw)
        io = PDBIO()
        io.set_structure(structure)
        
        class NotWater(Select):
            def accept_residue(self, r): return r.get_resname() not in ["HOH", "WAT"]
            
        # Determine output path
        if mission_id:
            out_dir = os.path.join(PROJECTS_DIR, mission_id)
            os.makedirs(out_dir, exist_ok=True)
            clean_path = os.path.join(out_dir, f"{pdb_id}_clean.pdb")
        else:
            clean_path = os.path.join(PDB_DIR, f"{pdb_id}_clean.pdb")
            
        io.save(clean_path, NotWater())
        return clean_path

    @staticmethod
    def find_pocket_p2rank(pdb_path: str, pdb_id: str):
        # Fallback if P2Rank isn't configured or found
        if not os.path.exists(P2RANK_EXEC):
            log_msg("[P2Rank] Executable not found. Using Center of Mass Fallback.")
            return ProteinManager.find_pocket_center_native(pdb_path)

        out_dir = os.path.join(PDB_DIR, "p2rank_out")
        os.makedirs(out_dir, exist_ok=True)
        
        # Determine shell usage based on OS
        use_shell = (os.name == 'nt')
        
        cmd = [P2RANK_EXEC, "predict", "-f", pdb_path, "-o", out_dir, "-threads", "1", "-visualizations", "0"]
        try:
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
            # Ray trace slightly for better look
            pymol.cmd.png(os.path.join(IMG_DIR, fname), width=600, height=400, ray=1)
            return f"http://localhost:8000/static/{fname}"
        except: return None


import uuid
class RetroManager:
    @staticmethod
    def run_retrosynthesis(smiles: str):
        print(f"\n[BACKEND] STARTING RETROSYNTHESIS FOR: {smiles}", flush=True)

        # 1. Validation
        if not os.path.exists(AIZYNTH_ENV_PYTHON):
            return {"error": f"AiZynth Python not found at {AIZYNTH_ENV_PYTHON}"}
        
        abs_config_path = os.path.abspath(AIZYNTH_CONFIG)
        if not os.path.exists(abs_config_path):
            return {"error": f"AiZynth Config not found at {abs_config_path}"}

        # 2. Define Paths for this specific run
        run_id = str(uuid.uuid4())[:8]
        json_path = os.path.join(tempfile.gettempdir(), f"retro_{run_id}.json")
        from config import DATA_DIR 
        
        # 3. CONSTRUCT THE WORKER SCRIPT
        worker_code = f"""
import subprocess
import json
import os
import sys

try:
    from aizynthfinder.reactiontree import ReactionTree
except ImportError:
    print("CRITICAL: aizynthfinder not installed in this environment.")
    sys.exit(1)

CONFIG_FILE = r"{abs_config_path}"
TARGET_SMILES = "{smiles}"
OUTPUT_JSON = r"{json_path}"
IMAGE_FOLDER = r"{IMG_DIR}"
AIZYNTH_CLI = r"{AIZYNTH_CLI}"

command = [
    AIZYNTH_CLI,
    "--config", CONFIG_FILE,
    "--smiles", TARGET_SMILES,
    "--output", OUTPUT_JSON
]

try:
    subprocess.run(command, check=True)
except subprocess.CalledProcessError as e:
    sys.exit(1)

if not os.path.exists(OUTPUT_JSON):
    sys.exit(1)

with open(OUTPUT_JSON, 'r') as f:
    data = json.load(f)

routes = []
if isinstance(data, dict) and 'routes' in data:
    routes = data['routes']
elif isinstance(data, list):
    routes = data

processed_routes = []
for i, route_dict in enumerate(routes[:5]):
    try:
        tree = ReactionTree.from_dict(route_dict)
        img_filename = f"retro_{run_id}_{{i}}.png"
        img_path = os.path.join(r"{DATA_DIR}", img_filename) 
        
        try:
            img = tree.to_image(in_stock_colors={{True: "#99ff99", False: "#ffcc99"}})
            img.save(img_path)
            fname = img_filename
        except:
            fname = ""

        steps = []
        for reaction in tree.reactions():
            rxn_reactants = []
            for mol in reaction.reactants:
                if hasattr(mol, "smiles"): rxn_reactants.append(mol.smiles)
                elif isinstance(mol, tuple) and len(mol) > 0 and hasattr(mol[0], "smiles"):
                    rxn_reactants.append(mol[0].smiles)
                else: rxn_reactants.append(str(mol))
            
            steps.append({{
                "reaction_smarts": getattr(reaction, "smarts", "N/A"),
                "reactants": " + ".join(rxn_reactants),
                # Add individual SMILES for rendering in UI
                "reactant_smiles": [m.smiles if hasattr(m, "smiles") else str(m) for m in reaction.reactants]
            }})

        processed_routes.append({{
            "image": fname,
            "steps": steps,
            "score": round(float(route_dict.get("scores", {{}}).get("state_score", 0)), 3)
        }})
    except Exception as e:
        print(f"Error: {{e}}")

with open(OUTPUT_JSON, "w") as f:
    json.dump({{"routes": processed_routes}}, f)
"""

        script_path = os.path.join(tempfile.gettempdir(), f"worker_{run_id}.py")
        try:
            with open(script_path, 'w') as f:
                f.write(worker_code)

            process = subprocess.Popen(
                [AIZYNTH_ENV_PYTHON, script_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=BASE_DIR
            )
            
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                return {"error": "AiZynth process failed."}

            if os.path.exists(json_path):
                with open(json_path, 'r') as f:
                    results = json.load(f)
                
                for route in results.get("routes", []):
                    if route.get("image"):
                        route["image_url"] = f"http://localhost:8000/static/{route['image']}"
                
                return results
            return {"error": "No results generated"}

        finally:
            if os.path.exists(script_path): os.remove(script_path)
            if os.path.exists(json_path): os.remove(json_path)



def fetch_chembl_csv(target_name: str, mission_id: str):
    """Searches ChEMBL, extracts IC50 data, and saves as CSV in project directory."""
    if not CHEMBL_AVAILABLE: 
        return {"error": "ChEMBL client not available."}
        
    from config import PROJECTS_DIR
    out_dir = os.path.join(PROJECTS_DIR, mission_id)
    os.makedirs(out_dir, exist_ok=True)
    csv_filename = f"{mission_id}_chembl_data.csv"
    csv_path = os.path.join(out_dir, csv_filename)

    log_msg(f"   [ChEMBL] Searching for {target_name}...")
    
    # 1. Find Target
    targets = new_client.target.filter(target_synonym__icontains=target_name).filter(target_type='SINGLE PROTEIN')
    if not targets:
        targets = new_client.target.filter(pref_name__icontains=target_name)
    
    if not targets:
        return {"error": f"No targets found in ChEMBL for '{target_name}'"}
    
    target_id = targets[0]['target_chembl_id']
    target_name_found = targets[0].get('pref_name', target_name)
    log_msg(f"   [ChEMBL] Found target: {target_name_found} ({target_id})")

    # 2. Get Bioactivities
    acts = new_client.activity.filter(
        target_chembl_id=target_id, 
        standard_type="IC50", 
        standard_units="nM", 
        standard_value__isnull=False
    ).only(['canonical_smiles', 'standard_value', 'standard_type', 'standard_units', 'molecule_chembl_id'])
    
    data = []
    for a in acts:
        if len(data) >= 1000: break
        smiles = a.get('canonical_smiles')
        if not smiles: continue
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol: continue
            
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            qed_val = QED.qed(mol)
            sa_score = calculateScore(mol)
            
            data.append({
                "molecule_id": a.get('molecule_chembl_id'),
                "smiles": smiles,
                "value": float(a['standard_value']),
                "type": a['standard_type'],
                "units": a['standard_units'],
                "mw": round(mw, 2),
                "logp": round(logp, 2),
                "qed": round(qed_val, 3),
                "sascore": round(sa_score, 3)
            })
        except Exception:
            continue
    
    if not data:
        return {"error": "No IC50 data found for this target."}

    # 3. Save to CSV
    df = pd.DataFrame(data)
    df['pIC50'] = -np.log10(df['value'] * 1e-9)
    df.to_csv(csv_path, index=False)
    
    def get_hist(vals, n_bins=15):
        h, e = np.histogram(vals.dropna(), bins=n_bins)
        return [{"bin": round(float(e[i]), 2), "count": int(count)} for i, count in enumerate(h)]

    # Calculate distributions for all major properties
    pic50_bins = get_hist(df['pIC50'])
    mw_bins = get_hist(df['mw'])
    logp_bins = get_hist(df['logp'])
    qed_bins = get_hist(df['qed'])
    sascore_bins = get_hist(df['sascore'])

    # Top 50 candidates for dynamic UI filtering
    top_candidates = df.sort_values('pIC50', ascending=False).head(50).to_dict('records')

    log_msg(f"   [ChEMBL] Saved {len(df)} records with full property distributions.")
    
    return {
        "csv_file": csv_filename,
        "count": len(df),
        "target_id": target_id,
        "target_name": target_name_found,
        "pIC50_mean": float(df['pIC50'].mean()),
        "pIC50_max": float(df['pIC50'].max()),
        "mw_mean": float(df['mw'].mean()),
        "logp_mean": float(df['logp'].mean()),
        "pic50_bins": pic50_bins,
        "mw_bins": mw_bins,
        "logp_bins": logp_bins,
        "qed_bins": qed_bins,
        "sascore_bins": sascore_bins,
        "top_candidates": top_candidates
    }