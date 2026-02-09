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
    """Converts Vina PDBQT to PDB with unique atom names for PLIP."""
    element_counts = defaultdict(int)
    with open(pdbqt_path, 'r') as f_in, open(pdb_path, 'w') as f_out:
        atom_index = 1
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except: continue

                raw_atom = line[12:16].strip()
                # Vina puts element in 77-79, if empty use atom name
                element = line[77:79].strip() 
                if not element: element = raw_atom[0]
                
                element_counts[element] += 1
                unique_name = f"{element}{element_counts[element]}"
                
                # Write standard PDB line
                f_out.write(
                    f"HETATM{atom_index:5d} {unique_name:<4} LIG L   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element:>2}\n"
                )
                atom_index += 1
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
        # 1. Setup Project Paths
        run_id = str(uuid.uuid4())[:8]
        # output: data/projects/{project_id}/docking/{run_id}
        output_dir = os.path.join(PROJECTS_DIR, project_id, "docking", run_id)
        os.makedirs(output_dir, exist_ok=True)

        ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
        vina_out_pdbqt = os.path.join(output_dir, "vina_out.pdbqt")
        ligand_pdb = os.path.join(output_dir, "ligand.pdb")
        complex_pdb = os.path.join(output_dir, "complex.pdb")

        try:
            # 2. Prepare Ligand (Meeko)
            mol = Chem.MolFromSmiles(smiles)
            if not mol: raise ValueError("Invalid SMILES string")
            
            m = Chem.AddHs(mol)
            AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
            
            prep = MoleculePreparation()
            mol_setups = prep.prepare(m)
            with open(ligand_pdbqt, 'w') as f:
                f.write(PDBQTWriterLegacy().write_string(mol_setups[0])[0])

            # 3. Calculate Center 
            center = calculate_pocket_center(DEFAULT_POCKET_PDB)

            # 4. Run Vina
            # Check if Vina exists
            if not os.path.exists(VINA_EXE):
                raise FileNotFoundError(f"Vina executable not found at: {VINA_EXE}")

            cmd = [
                VINA_EXE,
                "--receptor", DEFAULT_RECEPTOR_PDBQT,
                "--ligand", ligand_pdbqt,
                "--center_x", str(center["x"]), 
                "--center_y", str(center["y"]), 
                "--center_z", str(center["z"]),
                "--size_x", DOCKING_CONFIG["SIZE"]["x"], 
                "--size_y", DOCKING_CONFIG["SIZE"]["y"], 
                "--size_z", DOCKING_CONFIG["SIZE"]["z"],
                "--out", vina_out_pdbqt,
                "--cpu", "4", "--exhaustiveness", "8"
            ]
            
            subprocess.run(cmd, capture_output=True, check=True)

            # 5. Extract Score & Convert
            score = 0.0
            if os.path.exists(vina_out_pdbqt):
                with open(vina_out_pdbqt, 'r') as f:
                    for line in f:
                        if "REMARK VINA RESULT:" in line:
                            parts = line.split()
                            if len(parts) > 3:
                                score = float(parts[3])
                            break
                
                pdbqt_to_pdb(vina_out_pdbqt, ligand_pdb)
                merge_complex(DEFAULT_RECEPTOR_PDB, ligand_pdb, complex_pdb)

            # 6. Return Data 
            # Returns a URL relative to the FastAPI static mount
            relative_url = f"/static/projects/{project_id}/docking/{run_id}/complex.pdb"
            
            return {
                "success": True,
                "score": score,
                "file_url": relative_url,
                "run_id": run_id
            }

        except subprocess.CalledProcessError as e:
            print(f"Vina STDERR: {e.stderr.decode('utf-8') if e.stderr else 'No stderr'}")
            return {"success": False, "error": "Docking executable failed"}
        except Exception as e:
            print(f"Docking Error: {e}")
            traceback.print_exc()
            return {"success": False, "error": str(e)}