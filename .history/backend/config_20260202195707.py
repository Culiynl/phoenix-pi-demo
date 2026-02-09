import os
import sys
import logging

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

# --- DIRECTORIES ---
DATA_DIR = os.path.join(BASE_DIR, "data")
PROJECTS_DIR = os.path.join(DATA_DIR, "projects")
IMG_DIR = os.path.join(DATA_DIR, "images")

RECEPTOR_DIR = os.path.join(DATA_DIR, "receptor")
VINA_EXE = os.path.join(BASE_DIR, "vina", "vina_1.2.7_win.exe") # Adjust name based on OS

DEFAULT_RECEPTOR_PDBQT = os.path.join(RECEPTOR_DIR, "6od6_protein_dry.pdbqt")
DEFAULT_RECEPTOR_PDB = os.path.join(RECEPTOR_DIR, "6od6_clean.pdb")
DEFAULT_POCKET_PDB = os.path.join(RECEPTOR_DIR, "6od6_pocket_clean.pdb")

DOCKING_CONFIG = {
    "SIZE": {"x": "22", "y": "22", "z": "22"},
    "CENTER": {"x": 12.34, "y": 45.67, "z": 89.10} # Will be auto-calculated if None
}

for d in [DATA_DIR, PROJECTS_DIR, IMG_DIR, RECEPTOR_DIR]:
    os.makedirs(d, exist_ok=True)

# --- CONFIGURATION ---
GOOGLE_API_KEY = "AIzaSyBIAIbvC5EKpQVRsGnuPGAjVQD_wkg1yEc" 
ENTREZ_EMAIL = "your.email@example.com"

# --- DIRECTORIES ---
CACHE_DIR = os.path.join(BASE_DIR, "data", "gnn_training_cache")
PDB_DIR = os.path.join(BASE_DIR, "data", "pdb_files")
IMG_DIR = os.path.join(BASE_DIR, "data", "images")
MODEL_DIR = os.path.join(BASE_DIR, "data", "saved_models")

# --- EXTERNAL TOOLS ---
P2RANK_DIR = os.path.join(BASE_DIR, "p2rank_2.5.1")
P2RANK_EXEC = os.path.join(P2RANK_DIR, "prank.bat") if os.name == 'nt' else os.path.join(P2RANK_DIR, "prank")

# --- AIZYNTHFINDER CONFIGURATION ---
# The folder containing your config.yml
AIZYNTH_FOLDER = os.path.join(BASE_DIR, "aizynthfinder") 
AIZYNTH_CONFIG = os.path.join(AIZYNTH_FOLDER, "config.yml")

# Executables (Use raw strings 'r' for Windows paths)
AIZYNTH_ENV_PYTHON = r"F:\PyMol\envs\aizynth-env\python.exe"
AIZYNTH_CLI = r"F:\PyMol\envs\aizynth-env\Scripts\aizynthcli.exe"

# Create directories
for d in [CACHE_DIR, PDB_DIR, IMG_DIR, MODEL_DIR, AIZYNTH_FOLDER]:
    os.makedirs(d, exist_ok=True)

# --- PATCH TRANSFORMERS ---
try:
    import transformers.utils.import_utils
    import transformers.modeling_utils
    transformers.utils.import_utils.check_torch_load_is_safe = lambda: None
    transformers.modeling_utils.check_torch_load_is_safe = lambda: None
except ImportError:
    pass