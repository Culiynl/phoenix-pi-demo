import os
import sys
import math
import subprocess
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from config import VINA_EXE, DEFAULT_RECEPTOR_PDBQT, DEFAULT_RECEPTOR_PDB, DEFAULT_POCKET_PDB

# Check PLIP installation
try:
    from plip.structure.preparation import PDBComplex
except ImportError:
    print("[WARNING] PLIP is not installed. Structural analysis will be limited.")
    PDBComplex = None

# Analysis Constraints
HEAVY_ATOM_CUTOFF = 10
KEY_RESIDUES = ["ASP32", "ASP228"]

# PLIP-based Thresholds (Angstroms)
THRESHOLDS = {
    "HYDROPH_DIST_MAX": 4.0,
    "HBOND_DIST_MAX": 4.1,
    "HALOGEN_DIST_MAX": 4.0,
    "SALTBRIDGE_DIST_MAX": 5.5
}

def get_heavy_atom_count(smiles):
    if not isinstance(smiles, str): return 0
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol.GetNumHeavyAtoms() if mol else 0
    except:
        return 0

def pdbqt_to_pdb(pdbqt_path, pdb_path, chain="L", resnum=999):
    element_counts = defaultdict(int)
    if not os.path.exists(pdbqt_path): return
    with open(pdbqt_path, 'r') as f_in, open(pdb_path, 'w') as f_out:
        atom_index = 1
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                try: x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                except ValueError: continue

                raw_atom = line[12:16].strip()
                element = line[77:79].strip()
                if not element: element = raw_atom[0]
                
                element_counts[element] += 1
                unique_name = f"{element}{element_counts[element]}"
                
                f_out.write(f"HETATM{atom_index:5d} {unique_name:<4} LIG {chain}{resnum:>4}    "
                            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2}\n")
                atom_index += 1
        f_out.write("END\n")

def merge_complex(protein_pdb, ligand_pdb, output_pdb):
    with open(output_pdb, 'w') as f_out:
        if os.path.exists(protein_pdb):
            with open(protein_pdb, 'r') as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")) and "LIG" not in line:
                        f_out.write(line)
            f_out.write("TER\n")
        
        if os.path.exists(ligand_pdb):
            with open(ligand_pdb, 'r') as f:
                for line in f:
                    if line.startswith("HETATM"):
                        f_out.write(line)
        f_out.write("END\n")

def run_plip(pdb_file):
    if PDBComplex is None:
        return {"error": "PLIP not installed"}
        
    try:
        mol = PDBComplex()
        mol.load_pdb(pdb_file, as_string=False)
        mol.analyze()
        
        results = []
        if not mol.interaction_sets:
            return results

        for bs_id, site in mol.interaction_sets.items():
            interaction_types = [
                ("Hydrophobic", ["hydrophobic_contacts"]), ("H-Bond", ["hbonds", "hbonds_pdon", "hbonds_ldon"]),
                ("Salt Bridge", ["saltbridges", "saltbridges_pneg", "saltbridges_lneg"]), ("Pi-Stacking", ["pistacking"]),
                ("Halogen", ["halogen_bonds"]), ("Water Bridge", ["water_bridges"]), ("Metal", ["metal_complexes"])
            ]
            
            for label, attrs in interaction_types:
                interactions = []
                for attr in attrs:
                    val = getattr(site, attr, [])
                    if val: interactions.extend(val)

                interactions = list(set(interactions))
                for i in interactions:
                    l_name = "UNK"
                    if hasattr(i, 'ligand_atom') and i.ligand_atom: l_name = i.ligand_atom.name
                    elif hasattr(i, 'ligand_atoms') and i.ligand_atoms: l_name = i.ligand_atoms[0].name
                    rt, rn = getattr(i, 'restype', 'UNK'), getattr(i, 'resnr', '?')
                    dist = getattr(i, 'distance', getattr(i, 'distance_ad', 0.0))
                    results.append({
                        "type": label,
                        "ligand_atom": l_name,
                        "residue": f"{rt}{rn}",
                        "distance": round(float(dist), 2)
                    })
        return results
    except Exception as e:
        return {"error": str(e)}

def parse_pdb_coordinates(pdb_path, target_residues=None, is_ligand=False):
    atoms = []
    if not os.path.exists(pdb_path): return atoms
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                try: res_seq = int(line[22:26].strip())
                except ValueError: continue
                atom_name, element = line[12:16].strip(), line[76:78].strip()
                if not element: element = atom_name[0]
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                res_id = f"{res_name}{res_seq}"

                if is_ligand and "LIG" in line:
                    atoms.append({'atom': atom_name, 'element': element, 'x': x, 'y': y, 'z': z})
                elif not is_ligand and "LIG" not in line:
                    if target_residues and res_id in target_residues:
                        atoms.append({'atom': atom_name, 'element': element, 'res_id': res_id, 'res_name': res_name, 'x': x, 'y': y, 'z': z})
    return atoms

def classify_interaction(l_atom, l_elem, p_atom, p_elem, p_resname, distance):
    interactions = []
    l, p = l_elem.upper(), p_elem.upper()
    
    # Matching User-Provided Logic Exactly
    if p_resname in ['ASP', 'GLU'] and p_atom in ['OD1', 'OD2', 'OE1', 'OE2'] and l == 'N' and distance <= THRESHOLDS["SALTBRIDGE_DIST_MAX"]:
        interactions.append("Salt Bridge")
    if l in ['O', 'N', 'F'] and p in ['O', 'N'] and distance <= THRESHOLDS["HBOND_DIST_MAX"]:
        interactions.append("H-Bond")
    if l == 'C' and p == 'C' and distance <= THRESHOLDS["HYDROPH_DIST_MAX"]:
        interactions.append("Hydrophobic")
        
    if distance <= 4.0 and not interactions:
        interactions.append("van der Waals")
        
    return interactions

def analyze_structural_interactions(complex_pdb, ligand_pdb):
    """Combines PLIP and manual residue distance analysis."""
    plip_results = run_plip(complex_pdb)
    
    prot_atoms = parse_pdb_coordinates(complex_pdb, target_residues=KEY_RESIDUES, is_ligand=False)
    lig_atoms = parse_pdb_coordinates(ligand_pdb, is_ligand=True)
    
    detailed_interactions = []
    for res_id in KEY_RESIDUES:
        res_atoms = [pa for pa in prot_atoms if pa['res_id'] == res_id]
        if not res_atoms: continue

        for l_atom in lig_atoms:
            for p_atom in res_atoms:
                dist = math.dist((l_atom['x'], l_atom['y'], l_atom['z']), (p_atom['x'], p_atom['y'], p_atom['z']))
                if dist <= THRESHOLDS["SALTBRIDGE_DIST_MAX"]:
                    itypes = classify_interaction(l_atom['atom'], l_atom['element'], p_atom['atom'], p_atom['element'], p_atom['res_name'], dist)
                    for itype in itypes:
                        detailed_interactions.append({
                            'lig_atom': l_atom['atom'], 
                            'residue': p_atom['res_id'],
                            'prot_atom': p_atom['atom'], 
                            'distance': round(dist, 3), 
                            'interaction': itype
                        })
    return {
        "plip": plip_results,
        "key_residues": detailed_interactions
    }
