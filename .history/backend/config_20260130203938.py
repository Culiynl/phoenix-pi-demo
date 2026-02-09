import os
import sys

# --- BASE PATHS ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

# --- DIRECTORIES ---
# Images will be saved here so the API can serve them
IMG_DIR = os.path.join(BASE_DIR, "data", "images") 
CACHE_DIR = os.path.join(BASE_DIR, "data", "gnn_training_cache")
PDB_DIR = os.path.join(BASE_DIR, "data", "pdb_files")
MODEL_DIR = os.path.join(BASE_DIR, "data", "saved_models")

# --- AIZYNTHFINDER CONFIGURATION ---
# The folder containing your config.yml
AIZYNTH_FOLDER = os.path.join(BASE_DIR, "aizynthfinder") 
AIZYNTH_CONFIG = os.path.join(AIZYNTH_FOLDER, "config.yml")

# Executables (Use raw strings 'r' for Windows paths)
AIZYNTH_ENV_PYTHON = r"F:\PyMol\envs\aizynth-env\python.exe"
AIZYNTH_CLI = r"F:\PyMol\envs\aizynth-env\Scripts\aizynthcli.exe"

# --- OTHER CONFIGS ---
GOOGLE_API_KEY = "YOUR_KEY_HERE"
ENTREZ_EMAIL = "your.email@example.com"
P2RANK_DIR = os.path.join(BASE_DIR, "p2rank_2.5.1")
P2RANK_EXEC = os.path.join(P2RANK_DIR, "prank.bat") if os.name == 'nt' else os.path.join(P2RANK_DIR, "prank")

# Create directories
for d in [CACHE_DIR, PDB_DIR, IMG_DIR, MODEL_DIR, AIZYNTH_FOLDER]:
    os.makedirs(d, exist_ok=True)

# Patching
try:
    import transformers.utils.import_utils
    transformers.utils.import_utils.check_torch_load_is_safe = lambda: None
except: pass