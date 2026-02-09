import os
import sys
import logging
import torch

# --- 1. BASE PATH SETUP ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# --- 2. MAIN DIRECTORIES ---
DATA_DIR = os.path.join(BASE_DIR, "data")
PROJECTS_DIR = os.path.join(DATA_DIR, "projects")
IMG_DIR = os.path.join(DATA_DIR, "images")
CACHE_DIR = os.path.join(DATA_DIR, "gnn_training_cache")
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")
MODEL_DIR = os.path.join(DATA_DIR, "saved_models")
MODEL_CHECKPOINTS = os.path.join(DATA_DIR, "checkpoints")

# --- 3. DOCKING & VINA CONFIG ---
RECEPTOR_DIR = os.path.join(DATA_DIR, "receptor")
VINA_EXE = os.path.join(BASE_DIR, "vina", "vina_1.2.7_win.exe")  # Adjust if needed

DEFAULT_RECEPTOR_PDBQT = os.path.join(RECEPTOR_DIR, "6od6_protein_dry.pdbqt")
DEFAULT_RECEPTOR_PDB = os.path.join(RECEPTOR_DIR, "6od6_clean.pdb")
DEFAULT_POCKET_PDB = os.path.join(RECEPTOR_DIR, "6od6_pocket_clean.pdb")

DOCKING_CONFIG = {
    "SIZE": {"x": "22", "y": "22", "z": "22"},
    "CENTER": {"x": None, "y": None, "z": None} 
}

P2RANK_DIR = os.path.join(BASE_DIR, "p2rank_2.5.1")
P2RANK_EXEC = os.path.join(P2RANK_DIR, "prank.bat") if os.name == 'nt' else os.path.join(P2RANK_DIR, "prank")
# Directory to store PDBs and P2Rank results
PROTEIN_STORAGE = os.path.join(DATA_DIR, "proteins")

# AIZynthFinder
AIZYNTH_FOLDER = os.path.join(BASE_DIR, "aizynthfinder") 
AIZYNTH_CONFIG = os.path.join(AIZYNTH_FOLDER, "config.yml")
AIZYNTH_ENV_PYTHON = r"F:\PyMol\envs\aizynth-env\python.exe"
AIZYNTH_CLI = r"F:\PyMol\envs\aizynth-env\Scripts\aizynthcli.exe"

# --- 5. API KEYS & SECRETS ---
GOOGLE_API_KEY = "AIzaSyBIAIbvC5EKpQVRsGnuPGAjVQD_wkg1yEc" 
ENTREZ_EMAIL = "your.email@example.com"

CHECKPOINT_PATH = os.path.join(MODEL_CHECKPOINTS, "feudal_model.pt")
VOCAB_PATH = os.path.join(MODEL_CHECKPOINTS, "master_vocab.json")

# --- 6. INITIALIZATION ---
# Create all necessary directories
for d in [DATA_DIR, PROJECTS_DIR, IMG_DIR, CACHE_DIR, PDB_DIR, MODEL_DIR, RECEPTOR_DIR, AIZYNTH_FOLDER]:
    os.makedirs(d, exist_ok=True)

# --- 7. PATCHES (Transformers/PyTorch) ---
try:
    import transformers.utils.import_utils
    import transformers.modeling_utils
    transformers.utils.import_utils.check_torch_load_is_safe = lambda: None
    transformers.modeling_utils.check_torch_load_is_safe = lambda: None
except ImportError:
    pass