import os
import subprocess
import traceback
import uuid
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

# Import config from the parent directory
# This works because api.py (the entry point) is in the parent directory
from config import (
    PROJECTS_DIR, VINA_EXE, DEFAULT_RECEPTOR_PDBQT, 
    DEFAULT_RECEPTOR_PDB, DEFAULT_POCKET_PDB, DOCKING_CONFIG
)

def calculate_pocket_center(pocket_pdb_path):
    coords = []
    # Return default if file doesn't exist
    if not os.path.exists(pocket_pdb_path): 
        return DOCKING_CONFIG["CENTER"]
    
    with open(pocket_pdb_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    # PDB Coords: X(30-38), Y(38-46), Z(46-54)
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError: continue
    
    if not coords: return DOCKING_CONFIG["CENTER"]
    n = len(coords)
    return {
        "x": sum(c[0] for c in coords) / n,
        "y": sum(c[1] for c in coords) / n,
        "z": sum(c[2] for c in coords) / n
    }

def pdbqt_to_pdb(pdbqt_path, pdb_path):
    """Converts Vina PDBQT to PDB while PRESERVING MODEL TAGS for multi-pose support."""
    element_counts = defaultdict(int)
    with open(pdbqt_path, 'r') as f_in, open(pdb_path, 'w') as f_out:
        atom_index = 1
        for line in f_in:
            # Preserve MODEL tags so Mol* knows there are multiple poses
            if line.startswith("MODEL"):
                f_out.write(line)
                element_counts = defaultdict(int) # Reset counts per model
                atom_index = 1
                continue
            
            if line.startswith("ENDMDL"):
                f_out.write(line)
                continue

            if line.startswith(("ATOM", "HETATM")):
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    raw_atom = line[12:16].strip()
                    element = line[77:79].strip() or raw_atom[0]
                    element_counts[element] += 1
                    unique_name = f"{element}{element_counts[element]}"
                    
                    f_out.write(
                        f"HETATM{atom_index:5d} {unique_name:<4} LIG L   1    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2}\n"
                    )
                    atom_index += 1
                except: continue
        f_out.write("END\n")

def merge_complex(protein_pdb, ligand_pdb, output_pdb):
    """Merges receptor and docked ligand into one PDB."""
    with open(output_pdb, 'w') as f_out:
        # Write Protein
        if os.path.exists(protein_pdb):
            with open(protein_pdb, 'r') as f:
                for line in f:
                    # Exclude existing ligands if needed, or just keep protein
                    if line.startswith(("ATOM", "HETATM")) and "LIG" not in line:
                        f_out.write(line)
            f_out.write("TER\n")
        
        # Write Ligand
        if os.path.exists(ligand_pdb):
            with open(ligand_pdb, 'r') as f:
                for line in f:
                    if line.startswith("HETATM"):
                        f_out.write(line)
        f_out.write("END\n")

class DockingManager:
    @staticmethod
    def run_docking(project_id: str, smiles: str):
        run_id = str(uuid.uuid4())[:8]
        output_dir = os.path.join(PROJECTS_DIR, project_id, "docking", run_id)
        os.makedirs(output_dir, exist_ok=True)

        ligand_pdbqt = os.path.join(output_dir, "ligand_vina.pdbqt")
        ligand_pdb = os.path.join(output_dir, "ligand_multi.pdb")

        try:
            # 1. Prepare Ligand
            mol = Chem.MolFromSmiles(smiles)
            m = Chem.AddHs(mol)
            AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
            prep = MoleculePreparation()
            mol_setups = prep.prepare(m)
            with open(os.path.join(output_dir, "ligand_prep.pdbqt"), 'w') as f:
                f.write(PDBQTWriterLegacy().write_string(mol_setups[0])[0])

            # 2. Run Vina
            if center_override:
                center = center_override
                print(f"[Docking] Using P2Rank Center Override: {center}")
            else:
                center = calculate_pocket_center(DEFAULT_POCKET_PDB)
            cmd = [
                VINA_EXE, "--receptor", DEFAULT_RECEPTOR_PDBQT,
                "--ligand", os.path.join(output_dir, "ligand_prep.pdbqt"),
                "--center_x", str(center["x"]), "--center_y", str(center["y"]), "--center_z", str(center["z"]),
                "--size_x", DOCKING_CONFIG["SIZE"]["x"], "--size_y", DOCKING_CONFIG["SIZE"]["y"], "--size_z", DOCKING_CONFIG["SIZE"]["z"],
                "--out", ligand_pdbqt
            ]
            subprocess.run(cmd, capture_output=True, check=True)

            scores = []
            if os.path.exists(ligand_pdbqt):
                with open(ligand_pdbqt, 'r') as f:
                    for line in f:
                        if "REMARK VINA RESULT:" in line:
                            # Line looks like: REMARK VINA RESULT:    -7.2      0.000      0.000
                            parts = line.split()
                            if len(parts) > 3:
                                scores.append(float(parts[3]))
                
                pdbqt_to_pdb(ligand_pdbqt, ligand_pdb)

            return {
                "success": True,
                "scores": scores,  # Return the REAL list of scores
                "ligand_url": f"/static/projects/{project_id}/docking/{run_id}/ligand_multi.pdb",
                "receptor_url": f"/static/receptor/6od6_clean.pdb",
                "run_id": run_id
            }
        except Exception as e:
            return {"success": False, "error": str(e)}