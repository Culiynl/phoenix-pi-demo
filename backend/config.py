import os
import sys
import logging

# Define Base Directory
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
import torch
# Determine device
if torch.cuda.is_available() and torch.cuda.device_count() > 0:
    DEVICE = torch.device("cuda")
else:
    DEVICE = torch.device("cpu")

# --- 2. MAIN DIRECTORIES ---
# Cloud Run filesystem is read-only except for /tmp
if os.environ.get("K_SERVICE"):
    DATA_DIR = "/tmp/data"
else:
    DATA_DIR = os.path.join(BASE_DIR, "data")

PROJECTS_DIR = os.path.join(DATA_DIR, "projects")
IMG_DIR = os.path.join(DATA_DIR, "images")
CACHE_DIR = os.path.join(DATA_DIR, "gnn_training_cache")
PDB_DIR = os.path.join(DATA_DIR, "pdb_files")
MODEL_DIR = os.path.join(DATA_DIR, "saved_models")
MODEL_CHECKPOINTS = os.path.join(DATA_DIR, "checkpoints")

# --- 3. DOCKING & VINA CONFIG ---
# --- 3. DOCKING & VINA CONFIG ---
RECEPTOR_DIR = os.path.join(DATA_DIR, "receptor")

# Detect Vina path based on OS or Environment
if os.environ.get("K_SERVICE"):
    # Cloud Run (Linux) - assumes apt-get install autodock-vina
    VINA_EXE = "/usr/bin/vina"
elif os.name == 'nt':
    # Local Windows
    VINA_EXE = os.path.join(BASE_DIR, "vina", "vina_1.2.7_win.exe")
else:
    # Local Linux/Mac
    VINA_EXE = "vina"

DEFAULT_RECEPTOR_PDBQT = os.path.join(RECEPTOR_DIR, "6od6_protein_dry.pdbqt")
DEFAULT_RECEPTOR_PDB = os.path.join(RECEPTOR_DIR, "6od6_clean.pdb")
DEFAULT_POCKET_PDB = os.path.join(RECEPTOR_DIR, "6od6_pocket_clean.pdb")

DOCKING_CONFIG = {
    "SIZE": {"x": "22", "y": "22", "z": "22"},
    "CENTER": {"x": None, "y": None, "z": None} 
}

P2RANK_DIR = os.path.join(BASE_DIR, "p2rank_2.5.1")

# Default to None if not found, preventing crash
P2RANK_EXEC = None

if os.name == 'nt':
    # Local Windows
    potential_exec = os.path.join(P2RANK_DIR, "prank.bat")
    if os.path.exists(potential_exec):
        P2RANK_EXEC = potential_exec
else:
    # Linux (Cloud Run)
    # Only set if actually present (unlikely unless we COPY it in Dockerfile)
    potential_exec = os.path.join(P2RANK_DIR, "prank")
    if os.path.exists(potential_exec):
        P2RANK_EXEC = potential_exec
        try:
            os.chmod(P2RANK_EXEC, 0o755)
        except OSError:
            pass

PROTEIN_STORAGE = os.path.join(DATA_DIR, "proteins")
for d in [PROTEIN_STORAGE]:
    os.makedirs(d, exist_ok=True)

# AIZynthFinder
AIZYNTH_FOLDER = os.path.join(BASE_DIR, "aizynthfinder") 
AIZYNTH_CONFIG = os.path.join(AIZYNTH_FOLDER, "config.yml")
AIZYNTH_ENV_PYTHON = r"F:\PyMol\envs\aizynth-env\python.exe"
AIZYNTH_CLI = r"F:\PyMol\envs\aizynth-env\Scripts\aizynthcli.exe"

from dotenv import load_dotenv
import os
load_dotenv() # This loads the .env file locally
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
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