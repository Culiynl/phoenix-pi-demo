import os
import subprocess
import json
import tempfile
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

from config import (
    BASE_DIR, PDB_DIR, P2RANK_EXEC, IMG_DIR, 
    AIZYNTH_ENV_PYTHON, AIZYNTH_CLI, AIZYNTH_CONFIG
)
from state import log_msg

class RetroManager:
    @staticmethod
    def run_retrosynthesis(smiles: str):
        print(f"\n🚀 [BACKEND] STARTING RETROSYNTHESIS FOR: {smiles}", flush=True)

        # 1. Validation
        if not os.path.exists(AIZYNTH_ENV_PYTHON):
            return {"error": f"AiZynth Python not found at {AIZYNTH_ENV_PYTHON}"}
        if not os.path.exists(AIZYNTH_CONFIG):
            return {"error": f"AiZynth Config not found at {AIZYNTH_CONFIG}"}

        # 2. Define Paths for this specific run
        # We use a temp JSON file to capture the output data
        json_fd, json_path = tempfile.mkstemp(suffix=".json")
        os.close(json_fd) # Close file descriptor so subprocess can write to it
        
        # We point images to the public static folder
        # We assume standard naming: route_1.png, etc.
        # Clean up old images for this specific run if needed, or handle naming unique
        # For simplicity, we stick to your logic
        
        # 3. CONSTRUCT THE WORKER SCRIPT
        # We inject the variables from our API request into the script string
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
print(f"Found {{len(routes)}} routes.")

generated_images = []

for i, route_dict in enumerate(routes[:5]):
    try:
        print(f"  Visualizing Route {{i+1}}...")
        tree = ReactionTree.from_dict(route_dict)
        img = tree.to_image(in_stock_colors={{True: "#99ff99", False: "#ffcc99"}})
        
        # Create a unique filename based on smiles hash or timestamp ideally, 
        # but here we simply name them specifically
        fname = f"retro_{{i}}_{{os.path.basename(OUTPUT_JSON)}}.png"
        full_path = os.path.join(IMAGE_FOLDER, fname)
        
        img.save(full_path)
        generated_images.append(fname)
        print(f"    -> Saved: {{full_path}}")
    except Exception as e:
        print(f"    -> Error on Route {{i+1}}: {{e}}")

# Write a summary back to a new JSON or just print
print("Done.")
"""

        # 4. Write Temporary Script
        script_fd, script_path = tempfile.mkstemp(suffix=".py", mode="w+t")
        try:
            with os.fdopen(script_fd, 'w') as f:
                f.write(worker_code)

            print(f"⏳ Spawning AiZynth Process...", flush=True)
            
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
                
                # Add image links to the result so frontend can display them
                # The worker script saved them to IMG_DIR
                # We need to find which images correspond to this run
                # (Simple heuristic: matching the temp filename)
                base_name = os.path.basename(json_path)
                
                routes = results.get('routes', results)
                for i, route in enumerate(routes[:5]):
                    img_name = f"retro_{i}_{base_name}.png"
                    # Add a field for the API to serve
                    route['image_url'] = f"http://localhost:8000/static/{img_name}"

                return results
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