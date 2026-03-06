import torch
import os
import json
from config import (
    DEFAULT_RECEPTOR_PDBQT, DEFAULT_POCKET_PDB, DOCKING_CONFIG,
    VOCAB_PATH, DEVICE, CHECKPOINT_PATH, PROJECTS_DIR
)
from services.local_storage_service import _clean_serializable
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
                json.dump(_clean_serializable(logs), f)
        except: pass

    def run_optimization(self, project_id, start_smi, steps=5, batch_size=10, beam_width=2, min_hac=15):
        # Defensive casting for float inputs
        steps, batch_size, beam_width, min_hac = int(steps), int(batch_size), int(beam_width), int(min_hac)
        
        from concurrent.futures import ProcessPoolExecutor
        from rdkit import Chem
        import multiprocessing
        
        log_dir = os.path.join(PROJECTS_DIR, project_id, "feudal_logs")
        os.makedirs(log_dir, exist_ok=True)
        
        # Clear old progress logs
        with open(os.path.join(log_dir, "progress.json"), "w") as f: json.dump([], f)

        self.log_event(project_id, "Calculating dynamic pocket center...", "warn")
        center = self.calculate_pocket_center(DEFAULT_POCKET_PDB)
        self.log_event(project_id, f"Center identified: {center['x']}, {center['y']}", "success")

        # ⚡ DOCKING CACHE: Saves ~30% of docking calls by remembering SMILES
        docking_cache = {}
        
        self.log_event(project_id, f"Running initial docking for {start_smi[:30]}...", "info")
        _, init_score, status = run_vina_standalone(start_smi, center)
        
        if status != "Success":
            self.log_event(project_id, f"Initial docking failed: {status}", "error")
            return

        docking_cache[start_smi] = (init_score, status)
        
        # --- Initialize History & Beam ---
        initial_history = [{
            "step": 0,
            "smiles": start_smi,
            "score": init_score,
            "goal": None, 
            "status": "Initial"
        }]
        
        # Beam entry: (smiles, score, path_history)
        beam = [(start_smi, init_score, initial_history)]

        # Determine optimal workers
        max_workers = max(1, (os.cpu_count() or 4) - 1)
        
        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            for s in range(steps):
                self.log_event(project_id, f"Step {s+1}/{steps}: Generating candidates (Beam Width: {len(beam)})...", "info")
                
                new_beam_candidates = []
                all_futures = {}
                
                # Automatically carry over current beam (Elitism)
                for smi, score, history in beam:
                    new_beam_candidates.append((smi, score, history))
                
                # Branch out from each member of the beam
                for smi, score, history in beam:
                    candidates, step_goal = generate_candidates(
                        smi, self.manager, self.worker, self.vocab, self.norm, n_samples=batch_size
                    )
                    
                    if not candidates: continue
                    
                    for cand in candidates:
                        if cand in docking_cache:
                            # ⚡ FAST PATH: Instantly use cached score!
                            c_score, c_st = docking_cache[cand]
                            if c_st == "Success":
                                new_history = history + [{
                                    "step": s + 1, "smiles": cand, "score": c_score,
                                    "goal": step_goal, "status": "Beam (Cached)"
                                }]
                                new_beam_candidates.append((cand, c_score, new_history))
                        else:
                            # ⏳ SLOW PATH: Needs Vina evaluation
                            f = pool.submit(run_vina_standalone, cand, center)
                            all_futures[f] = (history, step_goal)
                
                if all_futures:
                    self.log_event(project_id, f"Step {s+1}: Docking {len(all_futures)} new candidates...", "info")
                    for f in all_futures:
                        parent_history, step_goal = all_futures[f]
                        try:
                            # Unpack correctly: smiles, score, status
                            c_smi, c_score, c_st = f.result(timeout=60)
                            docking_cache[c_smi] = (c_score, c_st)
                            
                            if c_st == "Success":
                                # 📏 HAC FILTER: Ensure molecule isn't too tiny
                                mol = Chem.MolFromSmiles(c_smi)
                                hac = mol.GetNumHeavyAtoms() if mol else 0
                                if hac < min_hac:
                                    self.log_event(project_id, f"Candidate {c_smi[:15]} rejected: HAC {hac} < min {min_hac}", "warn")
                                    continue

                                new_history = parent_history + [{
                                    "step": s + 1, "smiles": c_smi, "score": c_score,
                                    "goal": step_goal, "status": "Beam"
                                }]
                                new_beam_candidates.append((c_smi, c_score, new_history))
                        except Exception as e:
                            self.log_event(project_id, f"Docking failure: {str(e)}", "warn")

                # Sort by Best Score (lowest energy)
                new_beam_candidates.sort(key=lambda x: x[1])
                
                # Deduplicate SMILES and prune to BEAM_WIDTH
                unique_beam = []
                seen_smiles = set()
                for item in new_beam_candidates:
                    if item[0] not in seen_smiles:
                        seen_smiles.add(item[0])
                        unique_beam.append(item)
                        if len(unique_beam) == beam_width:
                            break
                            
                beam = unique_beam
                
                # Update progress with the best of this step
                best_this_step = beam[0]
                self.log_event(project_id, f"Step {s+1} Best: {best_this_step[1]:.2f} kcal/mol", "success")

                # Save intermediate history for UI updates
                with open(os.path.join(log_dir, f"run_{project_id}.json"), "w") as f:
                    json.dump(_clean_serializable(best_this_step[2]), f)

        self.log_event(project_id, "Beam Search Optimization complete.", "success")

    def get_history(self, project_id: str):
        """Helper for Pareto ranking to load the saved JSON."""
        log_path = os.path.join(PROJECTS_DIR, project_id, "feudal_logs", f"run_{project_id}.json")
        if os.path.exists(log_path):
            with open(log_path, "r") as f:
                return json.load(f)
        return None