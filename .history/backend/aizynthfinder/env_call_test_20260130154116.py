import subprocess
import tempfile
import os
import textwrap
import sys

# ------------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------------
AIZYNTH_PYTHON = r"F:\PyMol\envs\aizynth-env\python.exe"
CONFIG_FILE = "config.yml"
TARGET_SMILES = "O=C(O)C=C1C(CCC2=CC=CC=C2)=CC3=C1OCC3C"
OUTPUT_JSON = "results.json"
IMAGE_FOLDER = "visualizations"

# ------------------------------------------------------------------
# WORKER SCRIPT (RUNS INSIDE AIZYNTH ENV)
# ------------------------------------------------------------------
worker_code = f"""
import subprocess
import json
import os
from aizynthfinder.reactiontree import ReactionTree

CONFIG_FILE = r"{CONFIG_FILE}"
TARGET_SMILES = "{TARGET_SMILES}"
OUTPUT_JSON = r"{OUTPUT_JSON}"
IMAGE_FOLDER = r"{IMAGE_FOLDER}"

print(f"--- 1. Running AiZynthFinder Search ---")
print(f"Target: {{TARGET_SMILES}}")

# CALL AS PYTHON MODULE INSTEAD OF CLI EXECUTABLE
command = [
    "python", "-m", "aizynthfinder",
    "--config", CONFIG_FILE,
    "--smiles", TARGET_SMILES,
    "--output", OUTPUT_JSON
]

try:
    subprocess.run(command, check=True)
    print(">>> Search complete.")
except subprocess.CalledProcessError as e:
    print(f"Error running AiZynthFinder: {{e}}")
    exit()

print(f"--- 2. Generating Images ---")

if not os.path.exists(OUTPUT_JSON):
    print(f"Error: {{OUTPUT_JSON}} was not found.")
    exit()

os.makedirs(IMAGE_FOLDER, exist_ok=True)

with open(OUTPUT_JSON, 'r') as f:
    data = json.load(f)

routes = data['routes'] if 'routes' in data else data
print(f"Found {{len(routes)}} routes.")

for i, route_dict in enumerate(routes[:5]):
    try:
        print(f"  Visualizing Route {{i+1}}...")
        tree = ReactionTree.from_dict(route_dict)
        img = tree.to_image(in_stock_colors={{True: "#99ff99", False: "#ffcc99"}})
        filename = os.path.join(IMAGE_FOLDER, f"route_{{i+1}}.png")
        img.save(filename)
        print(f"    -> Saved: {{filename}}")
    except Exception as e:
        print(f"    -> Error on Route {{i+1}}: {{e}}")

print(f"Done! Images saved to '{{IMAGE_FOLDER}}'")
"""

# ------------------------------------------------------------------
# WRITE TEMP WORKER FILE
# ------------------------------------------------------------------
with tempfile.NamedTemporaryFile(delete=False, suffix=".py", mode="w") as temp_script:
    temp_script.write(textwrap.dedent(worker_code))
    temp_script_path = temp_script.name

print(f"Created temporary worker script at: {temp_script_path}")

# ------------------------------------------------------------------
# RUN USING THE TARGET ENVIRONMENT PYTHON
# ------------------------------------------------------------------
try:
    print("Launching AiZynthFinder environment Python...\n")
    subprocess.run([AIZYNTH_PYTHON, temp_script_path], check=True)
except subprocess.CalledProcessError as e:
    print(f"Pipeline failed: {e}")
    sys.exit(1)
finally:
    os.remove(temp_script_path)
    print("Temporary worker script removed.")
