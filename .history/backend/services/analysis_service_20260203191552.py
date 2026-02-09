from rdkit import Chem, DataStructs
from rdkit.Chem import QED, Descriptors, AllChem
from .sascorer import calculateScore

def calculate_pareto_winner(history, weights):
    """
    weights = {"affinity": float, "qed": float, "sa": float, "similarity": float}
    """
    if not history: return None
    
    # 1. Identify the starting molecule for similarity comparison
    start_smi = history[0]['smiles']
    start_mol = Chem.MolFromSmiles(start_smi)
    start_fp = AllChem.GetMorganFingerprintAsBitVect(start_mol, 2) if start_mol else None

    processed = []
    for item in history:
        mol = Chem.MolFromSmiles(item['smiles'])
        if not mol: continue
        
        # Calculate scientific metrics
        qed_val = QED.qed(mol)
        sa_val = calculateScore(mol)
        
        # Calculate Tanimoto Similarity (0.0 to 1.0)
        curr_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        similarity = DataStructs.TanimotoSimilarity(start_fp, curr_fp) if start_fp else 0.0

        # Normalize score (lower kcal/mol is better, so we negate for the "higher is better" Pareto)
        # If score is 0 (failed docking), we penalize it heavily (-1.0)
        affinity_score = -item['score'] if item['score'] < 0 else -10.0 

        processed.append({
            **item,
            "qed": round(qed_val, 3),
            "sa": round(sa_val, 3),
            "similarity": round(similarity, 3),
            "affinity_rank": affinity_score
        })

    # 2. Ranking Logic
    for p in processed:
        # We want: High affinity (negative score), High QED, Low SA, and User-defined Similarity
        # SA is inverted (10 - SA) because lower SA is better
        p['final_score'] = (
            (p['affinity_rank'] * weights.get('affinity', 1.0)) +
            (p['qed'] * weights.get('qed', 1.0) * 5) + # Boost QED to match affinity scale
            (p['similarity'] * weights.get('similarity', 1.0) * 5) +
            ((10 - p['sa']) * weights.get('sa', 1.0) * 0.5)
        )
    
    # Return the one with the highest weighted score
    return sorted(processed, key=lambda x: x['final_score'], reverse=True)[0]