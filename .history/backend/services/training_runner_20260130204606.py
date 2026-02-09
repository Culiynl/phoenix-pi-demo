# import os
# import torch
# import pandas as pd
# import numpy as np
# import optuna
# import torch.nn.functional as F
# from torch_geometric.loader import DataLoader
# from torch_geometric.data import Batch
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import r2_score
# import tempfile

# from config import CACHE_DIR, MODEL_DIR, BASE_DIR
# from state import training_state, log_msg, reset_state
# from services.bio_service import ProteinManager

# # --- IMPORT TRAINING MODULE ---
# # (Relies on config.py having set sys.path)
# try:
#     from train_llm_with_pic50 import (
#         PIC50Oracle, PocketLigandBindingModel, run_vina,
#         get_plip_interactions, create_interaction_graph, 
#         get_decoration_pocket_graph, get_decoration_features,
#         DEVICE, PROTEIN_FEATURE_DIM, INTERACTION_EMBED_DIM
#     )
#     import train_llm_with_pic50 # Import module to set globals
# except ImportError as e:
#     print(f"❌ CRITICAL ERROR: Could not import 'train_llm_with_pic50.py': {e}")
#     # Define fallbacks if needed or raise
#     raise e

# class GNNCacheDataset(torch.utils.data.Dataset):
#     def __init__(self, f): self.f = f
#     def __len__(self): return len(self.f)
#     def __getitem__(self, i): return torch.load(self.f[i], map_location='cpu')

# def collate_fn(batch):
#     out = {}
#     for k in ['l_embeds', 's_mask', 'a_tags', 'd_ids', 'l_mask', 'y']:
#         out[k] = torch.stack([b[k] for b in batch])
#     out['p_data'] = Batch.from_data_list([b['p_data'] for b in batch])
#     out['int_graph'] = Batch.from_data_list([b['int_graph'] for b in batch])
#     out['deco_graph'] = Batch.from_data_list([b['deco_graph'] for b in batch])
#     return out

# def run_training_task(csv_path: str, pdb_id: str):
#     reset_state()
#     model_save_path = os.path.join(MODEL_DIR, f"{pdb_id}.pt")
    
#     try:
#         log_msg(f">>> STARTING FOR PDB: {pdb_id}")
        
#         # 1. PREP & POCKET
#         training_state["status_text"] = "Prep & Pocket Detection..."
#         clean_pdb = ProteinManager.download_and_clean_pdb(pdb_id)
#         center = ProteinManager.find_pocket_p2rank(clean_pdb, pdb_id)
#         training_state["current_pocket_center"] = center
        
#         img = ProteinManager.generate_pocket_viz(clean_pdb, center, pdb_id)
#         if img: 
#             training_state["viz_image"] = img
#             log_msg("Visualization generated.")

#         # Set Global Variables in Train Script for Vina
#         train_llm_with_pic50.RECEPTOR_PDB_DRY = clean_pdb
#         train_llm_with_pic50.POCKET_CENTER = center

#         # 2. LOAD DATA
#         training_state["status_text"] = "Processing CSV..."
#         df = pd.read_csv(csv_path)
#         df.columns = [c.lower().strip() for c in df.columns]
        
#         s_col = next((c for c in df.columns if 'smile' in c), None)
#         v_col = next((c for c in df.columns if any(x in c for x in ['standard_value','pic50','ic50','value'])), None)
        
#         if not s_col or not v_col: raise Exception(f"CSV missing columns. Found: {df.columns}")
#         df = df.dropna(subset=[s_col, v_col])
#         training_state["total_ligands"] = len(df)
#         log_msg(f"Loaded {len(df)} ligands.")

#         # 3. DOCK & CACHE
#         oracle = PIC50Oracle(None)
#         files = []
#         MAX_ITEMS = 100 
        
#         for i, row in df.iterrows():
#             if len(files) >= MAX_ITEMS: break
#             cf = os.path.join(CACHE_DIR, f"{pdb_id}_{i}.pt")
            
#             try:
#                 raw_val = float(row[v_col])
#                 val = 9.0 - np.log10(raw_val * 1e-9) if raw_val > 50 else raw_val
#             except: continue

#             if not os.path.exists(cf):
#                 training_state["status_text"] = f"Docking {i+1}/{min(len(df), MAX_ITEMS)}"
#                 training_state["progress"] = int((len(files)/MAX_ITEMS)*30)
                
#                 try:
#                     smi = row[s_col]
#                     score, pose = run_vina(smi, center=center)
#                     if not pose: continue
                    
#                     inputs = oracle.lig_tok(smi, return_tensors="pt", truncation=True, padding="max_length", max_length=150)
#                     with torch.no_grad():
#                         l_emb = oracle.lig_mod(**inputs.to(DEVICE)).last_hidden_state[0]
                    
#                     with open(clean_pdb,'r') as f: prot = f.readlines()
#                     with open(pose,'r') as f: lig = [l for l in f if "ATOM" in l or "HETATM" in l]
                    
#                     tmp = os.path.join(tempfile.gettempdir(), f"tmp_{pdb_id}_{i}.pdb")
#                     with open(tmp,'w') as f: f.writelines(prot + lig + ["END\n"])
                    
#                     plip = get_plip_interactions(tmp)
#                     s_mask, a_tags, d_ids = get_decoration_features(smi, 150)
                    
#                     l_emb = l_emb[:150]
#                     l_mask = inputs['attention_mask'][0][:150]
                    
#                     data = {
#                         'p_data': oracle.pocket_graph_template,
#                         'l_embeds': l_emb.cpu(),
#                         's_mask': s_mask, 'a_tags': a_tags, 'd_ids': d_ids,
#                         'l_mask': l_mask.cpu(),
#                         'int_graph': create_interaction_graph(plip, smi),
#                         'deco_graph': get_decoration_pocket_graph(tmp, smi, d_ids, oracle.full_protein_info),
#                         'y': torch.tensor([val], dtype=torch.float)
#                     }
#                     torch.save(data, cf)
#                     if os.path.exists(pose): os.remove(pose)
#                     if os.path.exists(tmp): os.remove(tmp)
#                 except Exception as e:
#                     continue
            
#             files.append(cf)
            
#         if len(files) < 5: 
#             log_msg("Not enough valid data."); training_state["is_running"] = False; return

#         # 4. OPTUNA
#         training_state["status_text"] = "Optimizing Hyperparameters..."
#         train_f, val_f = train_test_split(files, test_size=0.2)
        
#         def objective(trial):
#             lr = trial.suggest_float("lr", 1e-4, 1e-3, log=True)
#             h_dim = trial.suggest_categorical("h_dim", [128, 256])
            
#             model = PocketLigandBindingModel(
#                 p_in_dim=PROTEIN_FEATURE_DIM, l_in_dim=768, 
#                 p_gcn_dims=[h_dim]*2, h_dim=h_dim, n_heads=4, drop=0.2, 
#                 num_attn_blocks=1, interaction_embed_dim=INTERACTION_EMBED_DIM, deco_pocket_gnn_dim=64
#             ).to(DEVICE)
            
#             opt = torch.optim.Adam(model.parameters(), lr=lr)
#             loader = DataLoader(GNNCacheDataset(train_f), batch_size=16, shuffle=True, collate_fn=collate_fn)
            
#             for _ in range(2): 
#                 model.train()
#                 for b in loader:
#                     opt.zero_grad()
#                     p = model(b['p_data'].to(DEVICE), b['l_embeds'].to(DEVICE), b['s_mask'].to(DEVICE), b['a_tags'].to(DEVICE), b['d_ids'].to(DEVICE), b['l_mask'].to(DEVICE), b['int_graph'].to(DEVICE), b['deco_graph'].to(DEVICE))
#                     loss = F.mse_loss(p.squeeze(), b['y'].to(DEVICE).squeeze())
#                     loss.backward(); opt.step()
#             return 0.5 

#         study = optuna.create_study(direction="maximize")
#         study.optimize(objective, n_trials=3)
#         bp = study.best_params
        
#         # 5. FINAL TRAIN
#         training_state["status_text"] = "Final Training..."
#         final_model = PocketLigandBindingModel(
#             p_in_dim=PROTEIN_FEATURE_DIM, l_in_dim=768, 
#             p_gcn_dims=[bp['h_dim']]*2, h_dim=bp['h_dim'], n_heads=4, drop=0.2, 
#             num_attn_blocks=1, interaction_embed_dim=INTERACTION_EMBED_DIM, deco_pocket_gnn_dim=64
#         ).to(DEVICE)
        
#         opt = torch.optim.Adam(final_model.parameters(), lr=bp['lr'])
#         loader = DataLoader(GNNCacheDataset(train_f), batch_size=16, shuffle=True, collate_fn=collate_fn)
#         v_loader = DataLoader(GNNCacheDataset(val_f), batch_size=16, collate_fn=collate_fn)
        
#         for ep in range(1, 11):
#             training_state["progress"] = 50 + int(ep*5)
#             final_model.train()
#             for b in loader:
#                 opt.zero_grad()
#                 p = final_model(b['p_data'].to(DEVICE), b['l_embeds'].to(DEVICE), b['s_mask'].to(DEVICE), b['a_tags'].to(DEVICE), b['d_ids'].to(DEVICE), b['l_mask'].to(DEVICE), b['int_graph'].to(DEVICE), b['deco_graph'].to(DEVICE))
#                 loss = F.mse_loss(p.squeeze(), b['y'].to(DEVICE).squeeze())
#                 loss.backward(); opt.step()
            
#             final_model.eval()
#             yt, yp = [], []
#             with torch.no_grad():
#                 for b in v_loader:
#                     p = final_model(b['p_data'].to(DEVICE), b['l_embeds'].to(DEVICE), b['s_mask'].to(DEVICE), b['a_tags'].to(DEVICE), b['d_ids'].to(DEVICE), b['l_mask'].to(DEVICE), b['int_graph'].to(DEVICE), b['deco_graph'].to(DEVICE))
#                     yt.extend(b['y'].tolist()); yp.extend(p.squeeze().tolist())
            
#             r2 = r2_score(yt, yp)
#             training_state["best_r2"] = max(training_state["best_r2"], r2)
#             log_msg(f"Epoch {ep} | Val R2: {r2:.3f}")
            
#             if r2 >= training_state["best_r2"]:
#                 torch.save(final_model.state_dict(), model_save_path)

#         log_msg("Done.")
#         training_state["progress"] = 100
#         training_state["status_text"] = "Complete"
        
#     except Exception as e:
#         log_msg(f"Error: {e}")
#         import traceback; traceback.print_exc()
#     finally:
#         training_state["is_running"] = False