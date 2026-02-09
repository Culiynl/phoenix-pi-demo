import os
import sys
import logging

# --- 1. CRITICAL PATH FIX ---
# Ensures 'train_llm_with_pic50.py' is importable
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

# --- CONFIGURATION ---
GOOGLE_API_KEY = "AIzaSyBIAIbvC5EKpQVRsGnuPGAjVQD_wkg1yEc"  # <--- REPLACE
ENTREZ_EMAIL = "your.email@example.com"

# --- DIRECTORIES ---
CACHE_DIR = os.path.join(BASE_DIR, "data", "gnn_training_cache")
PDB_DIR = os.path.join(BASE_DIR, "data", "pdb_files")
IMG_DIR = os.path.join(BASE_DIR, "data", "images")
MODEL_DIR = os.path.join(BASE_DIR, "data", "saved_models")

# --- EXTERNAL TOOLS ---
P2RANK_DIR = os.path.join(BASE_DIR, "p2rank_2.5.1")
P2RANK_EXEC = os.path.join(P2RANK_DIR, "prank.bat") if os.name == 'nt' else os.path.join(P2RANK_DIR, "prank")

AIZYNTH_PYTHON_EXE = r"F:\PyMol\envs\aizynth-env\python.exe"
AIZYNTH_CONFIG = os.path.join(BASE_DIR, "aizynthfinder", "config.yml")
AIZYNTH_RUNNER_SCRIPT = os.path.join(BASE_DIR, "aizynthfinder", "retro_runner.py")

# Create directories
for d in [CACHE_DIR, PDB_DIR, IMG_DIR, MODEL_DIR]:
    os.makedirs(d, exist_ok=True)

# --- PATCH TRANSFORMERS ---
try:
    import transformers.utils.import_utils
    import transformers.modeling_utils
    transformers.utils.import_utils.check_torch_load_is_safe = lambda: None
    transformers.modeling_utils.check_torch_load_is_safe = lambda: None
except ImportError:
    pass