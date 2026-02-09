from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, QED, Descriptors
from .sascorer import calculateScore
import pandas as pd
import numpy as np
import os

def calculate_pareto_winner(history, weights):
    if not history: return None
    
    # Starting molecule for similarity reference
    start_mol = Chem.MolFromSmiles(history[0]['smiles'])
    start_fp = AllChem.GetMorganFingerprintAsBitVect(start_mol, 2) if start_mol else None

    processed = []
    for item in history:
        mol = Chem.MolFromSmiles(item['smiles'])
        if not mol: continue
        
        # Calculate Metrics
        qed_val = QED.qed(mol)
        sa_val = calculateScore(mol)
        
        # Calculate Tanimoto Similarity
        curr_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        similarity = DataStructs.TanimotoSimilarity(start_fp, curr_fp) if start_fp else 0.0

        # Penalize failed docking (0.0)
        affinity_rank = -item['score'] if item['score'] < 0 else -15.0

        processed.append({
            **item,
            "qed": round(qed_val, 3),
            "sa": round(sa_val, 3),
            "similarity": round(similarity, 3),
            "affinity_rank": affinity_rank
        })

    # Weighted Scoring
    for p in processed:
        p['final_score'] = (
            (p['affinity_rank'] * weights.get('affinity', 1.0)) +
            (p['qed'] * weights.get('qed', 1.0) * 5.0) +
            (p['similarity'] * weights.get('similarity', 1.0) * 5.0) +
            ((10 - p['sa']) * weights.get('sa', 1.0) * 0.5)
        )
    
    winner = sorted(processed, key=lambda x: x['final_score'], reverse=True)[0]
    # Return full history so the graph can update
    winner['full_history'] = history 
    return winner

def find_best_molecule_from_csv(csv_path, weights):
    """
    Ranks molecules in a CSV by weighted properties and returns the best SMILES.
    Weights can include: mw, logp, qed, sascore, pIC50
    """
    if not os.path.exists(csv_path):
        return None
        
    df = pd.read_csv(csv_path)
    if df.empty:
        return None
        
    # Ensure all required columns exist or fill with defaults
    for col in ['mw', 'logp', 'qed', 'sascore', 'pIC50']:
        if col not in df.columns:
            df[col] = 0.0
            
    # Calculate weighted score (higher is better)
    # Note: for mw and sascore, lower is generally better, so we use negative weights or invert
    df['weighted_score'] = (
        (df['qed'] * weights.get('qed', 1.0) * 10) +
        (df['pIC50'] * weights.get('pIC50', 1.0)) +
        (df['logp'] * weights.get('logp', 0.0)) + # Often neutral unless specific target
        ((500 - df['mw']) * weights.get('mw', 0.01)) +
        ((10 - df['sascore']) * weights.get('sascore', 0.5))
    )
    
    best_row = df.sort_values(by='weighted_score', ascending=False).iloc[0]
    return {
        "smiles": best_row['smiles'],
        "score": float(best_row['weighted_score']),
        "molecule_id": best_row.get('molecule_id', 'Unknown'),
        "properties": {
            "mw": float(best_row['mw']),
            "logp": float(best_row['logp']),
            "qed": float(best_row['qed']),
            "sascore": float(best_row['sascore']),
            "pIC50": float(best_row['pIC50'])
        }
    }