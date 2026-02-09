import torch
import os
import json
from config import (
    CHECKPOINT_PATH, VOCAB_PATH, DEVICE, PROJECTS_DIR, 
    DEFAULT_RECEPTOR_PDBQT, DEFAULT_POCKET_PDB, DOCKING_CONFIG
)
from .feudal_models import (
    Manager, Worker, Vocab, DatasetNormalizer, 
    generate_candidates, run_vina_standalone
)

class FeudalOptimizationManager:
    _instance = None

    def __init__(self):
        self.vocab = Vocab(VOCAB_PATH)
        self.norm = DatasetNormalizer()
        self.manager = Manager().to(DEVICE)
        self.worker = Worker(self.vocab.n).to(DEVICE)
        
        # Determine vocab size mismatch protection
        ckpt = torch.load(CHECKPOINT_PATH, map_location=DEVICE)
        required_vocab_size = ckpt['worker_state']['emb.weight'].shape[0]
        if self.vocab.n != required_vocab_size:
             # Re-init worker if vocab mismatch detected
             self.worker = Worker(required_vocab_size).to(DEVICE)

        self.manager.load_state_dict(ckpt['manager_state'])
        self.worker.load_state_dict(ckpt['worker_state'])
        self.manager.eval()
        self.worker.eval()

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def calculate_pocket_center(self, pocket_pdb_path):
        """Helper to find the dynamic center of the pocket."""
        coords = []
        if not os.path.exists(pocket_pdb_path):
            print(f"[Feudal] Warning: Pocket file not found at {pocket_pdb_path}. Using config default.")
            return DOCKING_CONFIG["CENTER"]
        
        with open(pocket_pdb_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except ValueError: continue
        
        if not coords: return DOCKING_CONFIG["CENTER"]
        n = len(coords)
        return {
            "x": sum(c[0] for c in coords) / n,
            "y": sum(c[1] for c in coords) / n,
            "z": sum(c[2] for c in coords) / n
        }

    def run_optimization(self, project_id: str, start_smi: str, steps=3, batch_size=8):
        log_dir = os.path.join(PROJECTS_DIR, project_id, "feudal_logs")
        os.makedirs(log_dir, exist_ok=True)
        
        # 1. CALCULATE DYNAMIC CENTER ONCE
        center = self.calculate_pocket_center(DEFAULT_POCKET_PDB)
        print(f"[Feudal] Optimization starting. Dynamic Center: {center}")
        
        history = []
        curr_smi = start_smi
        
        # 2. INITIAL DOCKING
        _, init_score, _ = run_vina_standalone(start_smi, center)
        curr_score = init_score
        
        history.append({
            "step": 0, "smiles": curr_smi, "score": curr_score, "goal": None
        })

        for s in range(steps):
            # 3. GENERATE CANDIDATES
            candidates, goal = generate_candidates(
                curr_smi, self.manager, self.worker, self.vocab, self.norm, n_samples=batch_size
            )
            
            best_step_smi, best_step_score = curr_smi, curr_score
            
            # 4. DOCK ALL CANDIDATES AT THE CALCULATED CENTER
            for cand in candidates:
                _, score, status = run_vina_standalone(cand, center)
                if score < best_step_score:
                    best_step_score = score
                    best_step_smi = cand
            
            curr_smi, curr_score = best_step_smi, best_step_score
            history.append({
                "step": s + 1, "smiles": curr_smi, "score": curr_score, "goal": goal
            })

        # Save results to project folder
        log_path = os.path.join(log_dir, f"run_{project_id}.json")
        with open(log_path, "w") as f:
            json.dump(history, f)
            
        return history, log_path

    def get_history(self, project_id: str):
        """Helper for Pareto ranking to load the saved JSON."""
        log_path = os.path.join(PROJECTS_DIR, project_id, "feudal_logs", f"run_{project_id}.json")
        if os.path.exists(log_path):
            with open(log_path, "r") as f:
                return json.load(f)
        return None