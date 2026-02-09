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

    def log_event(self, project_id, message, type="info"):
        """Writes a log message to a file that the frontend polls."""
        log_path = os.path.join(PROJECTS_DIR, project_id, "feudal_logs", "progress.json")
        entry = {"time": os.popen('time /t').read().strip(), "msg": message, "type": type}
        
        try:
            logs = []
            if os.path.exists(log_path):
                with open(log_path, "r") as f:
                    logs = json.load(f)
            logs.append(entry)
            with open(log_path, "w") as f:
                json.dump(logs, f)
        except: pass

    def run_optimization(self, project_id, start_smi, steps=3, batch_size=8):
        log_dir = os.path.join(PROJECTS_DIR, project_id, "feudal_logs")
        os.makedirs(log_dir, exist_ok=True)
        
        # Clear old progress logs
        with open(os.path.join(log_dir, "progress.json"), "w") as f: json.dump([], f)

        self.log_event(project_id, "Calculating dynamic pocket center...", "warn")
        center = self.calculate_pocket_center(DEFAULT_POCKET_PDB)
        self.log_event(project_id, f"Center identified: {center['x']}, {center['y']}", "success")

        history = []
        curr_smi = start_smi
        
        self.log_event(project_id, "Running initial docking...", "info")
        _, init_score, _ = run_vina_standalone(start_smi, center)
        curr_score = init_score
        
        history.append({"step": 0, "smiles": curr_smi, "score": curr_score, "goal": None})

        for s in range(steps):
            self.log_event(project_id, f"Step {s+1}: Generating candidates...", "info")
            candidates, goal = generate_candidates(
                curr_smi, self.manager, self.worker, self.vocab, self.norm, n_samples=batch_size
            )
            
            self.log_event(project_id, f"Step {s+1}: Docking {len(candidates)} molecules...", "info")
            best_step_smi, best_step_score = curr_smi, curr_score
            
            for i, cand in enumerate(candidates):
                _, score, status = run_vina_standalone(cand, center)
                if score < best_step_score:
                    best_step_score = score
                    best_step_smi = cand
            
            if best_step_score < curr_score:
                self.log_event(project_id, f"Improvement found: {best_step_score} kcal/mol", "success")
            else:
                self.log_event(project_id, f"No improvement in Step {s+1}", "warn")

            curr_smi, curr_score = best_step_smi, best_step_score
            history.append({"step": s+1, "smiles": curr_smi, "score": curr_score, "goal": goal})

            # Save intermediate history so graph updates
            with open(os.path.join(log_dir, f"run_{project_id}.json"), "w") as f:
                json.dump(history, f)

        self.log_event(project_id, "Optimization complete.", "success")

    def get_history(self, project_id: str):
        """Helper for Pareto ranking to load the saved JSON."""
        log_path = os.path.join(PROJECTS_DIR, project_id, "feudal_logs", f"run_{project_id}.json")
        if os.path.exists(log_path):
            with open(log_path, "r") as f:
                return json.load(f)
        return None