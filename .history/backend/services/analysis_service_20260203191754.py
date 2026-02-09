from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, QED
from .sascorer import calculateScore

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