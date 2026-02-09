import torch
import os
import json
import pandas as pd
from config import (
    CHECKPOINT_PATH, VOCAB_PATH, DEVICE, PROJECTS_DIR, VINA_EXE, 
    DEFAULT_RECEPTOR_PDBQT, DEFAULT_POCKET_PDB
)
# Import the classes from your provided script
from .feudal_models import Manager, Worker, Vocab, DatasetNormalizer, generate_candidates, run_vina_standalone

class FeudalOptimizationManager:
    _instance = None

    def __init__(self):
        self.vocab = Vocab(VOCAB_PATH)
        self.norm = DatasetNormalizer()
        self.manager = Manager().to(DEVICE)
        self.worker = Worker(self.vocab.n).to(DEVICE)
        
        ckpt = torch.load(CHECKPOINT_PATH, map_location=DEVICE)
        self.manager.load_state_dict(ckpt['manager_state'])
        self.worker.load_state_dict(ckpt['worker_state'])
        self.manager.eval()
        self.worker.eval()

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def run_optimization(self, project_id: str, start_smi: str, steps=3, batch_size=8):
        # Setup specific log directory
        log_dir = os.path.join(PROJECTS_DIR, project_id, "feudal_logs")
        os.makedirs(log_dir, exist_ok=True)
        
        history = []
        curr_smi = start_smi
        
        # 1. Initial Score
        _, init_score, _ = run_vina_standalone(start_smi)
        curr_score = init_score
        
        history.append({
            "step": 0, "smiles": curr_smi, "score": curr_score, "goal": None
        })

        for s in range(steps):
            # Manager/Worker Inference
            candidates, goal = generate_candidates(
                curr_smi, self.manager, self.worker, self.vocab, self.norm, n_samples=batch_size
            )
            
            # Docking candidates
            best_step_smi, best_step_score = curr_smi, curr_score
            for cand in candidates:
                _, score, status = run_vina_standalone(cand)
                if score < best_step_score:
                    best_step_score = score
                    best_step_smi = cand
            
            curr_smi, curr_score = best_step_smi, best_step_score
            history.append({
                "step": s + 1, "smiles": curr_smi, "score": curr_score, "goal": goal
            })

        # Save results
        log_path = os.path.join(log_dir, f"run_{project_id}.json")
        with open(log_path, "w") as f:
            json.dump(history, f)
            
        return history, log_path