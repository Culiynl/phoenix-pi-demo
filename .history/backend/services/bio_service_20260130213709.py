import os
import subprocess
import json
import tempfile
import textwrap
import pandas as pd
import numpy as np
import re
import sys
from Bio.PDB import PDBList, PDBParser, PDBIO, Select

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
from aizynthfinder.env_call_test import OUTPUT_JSON
from config import (
    BASE_DIR, PDB_DIR, P2RANK_EXEC, IMG_DIR, 
    AIZYNTH_ENV_PYTHON, AIZYNTH_CLI, AIZYNTH_CONFIG
)
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


class RetroManager:
    @staticmethod
    def run_retrosynthesis(smiles: str):
        print(f"\n[BACKEND] STARTING RETROSYNTHESIS FOR: {smiles}", flush=True)

        # 1. Validation
        if not os.path.exists(AIZYNTH_ENV_PYTHON):
            return {"error": f"AiZynth Python not found at {AIZYNTH_ENV_PYTHON}"}
        if not os.path.exists(AIZYNTH_CONFIG):
            return {"error": f"AiZynth Config not found at {AIZYNTH_CONFIG}"}

        # 2. Define Paths for this specific run
        json_fd, json_path = tempfile.mkstemp(suffix=".json")
        os.close(json_fd) 
        
        # 3. CONSTRUCT THE WORKER SCRIPT
        worker_code = f"""
import subprocess
import json
import os
import sys

# Try importing, catch error if environment is wrong
try:
    from aizynthfinder.reactiontree import ReactionTree
except ImportError:
    print("CRITICAL: aizynthfinder not installed in this environment.")
    sys.exit(1)

CONFIG_FILE = r"{AIZYNTH_CONFIG}"
TARGET_SMILES = "{smiles}"
OUTPUT_JSON = r"{json_path}"
IMAGE_FOLDER = r"{IMG_DIR}"
AIZYNTH_CLI = r"{AIZYNTH_CLI}"

print(f"--- 1. Running AiZynthFinder Search ---")
print(f"Target: {{TARGET_SMILES}}")

command = [
    AIZYNTH_CLI,
    "--config", CONFIG_FILE,
    "--smiles", TARGET_SMILES,
    "--output", OUTPUT_JSON
]

try:
    # Run the CLI tool
    subprocess.run(command, check=True)
    print(">>> Search complete.")
except subprocess.CalledProcessError as e:
    print(f"Error running CLI: {{e}}")
    sys.exit(1)

print(f"--- 2. Generating Images ---")

if not os.path.exists(OUTPUT_JSON):
    print(f"Error: {{OUTPUT_JSON}} was not found.")
    sys.exit(1)

# Ensure image folder exists
os.makedirs(IMAGE_FOLDER, exist_ok=True)

with open(OUTPUT_JSON, 'r') as f:
    data = json.load(f)

# Handle list vs dict output format
routes = data['routes'] if 'routes' in data else data
routes = data['routes'] if isinstance(data, dict) and 'routes' in data else data

if not isinstance(routes, list):
    print("No valid routes found.")
    routes = []

from aizynthfinder.reactiontree import ReactionTree

processed_routes = []
for i, route_dict in enumerate(routes[:5]):
    try:
        print(f"  Processing Route {i+1}...")
        tree = ReactionTree.from_dict(route_dict)

        # --- SAVE IMAGE ---
        img = tree.to_image(in_stock_colors={True: "#99ff99", False: "#ffcc99"})
        fname = f'retro_{i}_{os.path.basename(OUTPUT_JSON)}.png'
        full_path = os.path.join(IMAGE_FOLDER, fname)
        img.save(full_path)

        # --- EXTRACT STEPS ---
        steps = []
        for reaction in tree.reactions():
            reactants = [mol.smiles for mol in reaction.reactants]
            steps.append({
                "reaction_smarts": reaction.smarts,
                "reactants": " + ".join(reactants)
            })

        processed_routes.append({
            "image": fname,
            "steps": steps
        })

        print(f"    -> Saved route {i+1}")

    except Exception as e:
        print(f"    -> Error on Route {i+1 if 'i' in locals() else '?'}: {e}")


# Overwrite output JSON with cleaned structure
with open(OUTPUT_JSON, "w") as f:
    json.dump({"routes": processed_routes}, f)

print("Done.")
"""

        # 4. Write Temporary Script
        script_fd, script_path = tempfile.mkstemp(suffix=".py")
        try:
            with os.fdopen(script_fd, 'w') as f:   # ← set mode HERE instead
                f.write(worker_code)

            print(f"Spawning AiZynth Process...", flush=True)
            
            # 5. Run Subprocess
            process = subprocess.Popen(
                [AIZYNTH_ENV_PYTHON, script_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=BASE_DIR
            )
            
            # Stream logs
            while True:
                line = process.stdout.readline()
                if not line and process.poll() is not None:
                    break
                if line:
                    print(f"[AiZynth] {line.strip()}", flush=True)
            
            stderr = process.stderr.read()
            if stderr:
                print(f"[AiZynth Error] {stderr}")

            if process.returncode != 0:
                return {"error": "AiZynth process failed. Check console logs."}

            # 6. Parse Results
            if os.path.exists(json_path):
                with open(json_path, 'r') as f:
                    results = json.load(f)

                base_name = os.path.basename(json_path)

                # --- FIX IS HERE ---
                if isinstance(results, dict) and "routes" in results:
                    routes = results["routes"]
                elif isinstance(results, list):
                    routes = results
                else:
                    return {"error": "Unexpected AiZynth output format"}

                clean_routes = []
                route_index = -1 

                for route_index, route in enumerate(routes[:5]):
                    img_name = f"retro_{route_index}_{base_name}.png"
                    route['image_url'] = f"http://localhost:8000/static/{img_name}"
                    clean_routes.append(route)

                return {
                    "is_solved": len(clean_routes) > 0,
                    "routes": clean_routes
                }


                return {"routes": routes}
            else:
                return {"error": "No output JSON generated."}


        except Exception as e:
            print(f"❌ EXCEPTION: {e}", flush=True)
            return {"error": str(e)}
        finally:
            # Cleanup
            if os.path.exists(script_path):
                os.remove(script_path)
            if os.path.exists(json_path):
                os.remove(json_path)


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