from rdkit.Chem import QED, Descriptors, Chem
from sascorer import calculateScore

def calculate_pareto_winner(results_list, weights):
    """
    weights = {"affinity": 1.0, "qed": 1.0, "sa": 1.0}
    """
    processed = []
    for item in results_list:
        mol = Chem.MolFromSmiles(item['smiles'])
        if not mol: continue
        
        qed_val = QED.qed(mol)
        sa_val = calculateScore(mol)
        # Affinity improvement is negative (lower is better), so we negate for ranking
        affinity = -item['score'] 

        processed.append({
            **item,
            "qed": qed_val,
            "sa": sa_val,
            "norm_affinity": affinity # simplified normalization logic here
        })

    # Apply weights and return best
    for p in processed:
        p['final_score'] = (p['norm_affinity'] * weights['affinity']) + \
                           (p['qed'] * weights['qed']) + \
                           ((1 - (p['sa']/10)) * weights['sa'])
    
    return sorted(processed, key=lambda x: x['final_score'], reverse=True)[0]