import os
import subprocess
import requests
import pandas as pd
import json
from config import PROTEIN_STORAGE, P2RANK_EXEC

class PDBManager:
    @staticmethod
    def fetch_from_rcsb(pdb_id: str):
        pdb_id = pdb_id.lower()
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        target_path = os.path.join(PROTEIN_STORAGE, f"{pdb_id}.pdb")
        
        response = requests.get(url)
        if response.status_code == 200:
            with open(target_path, "wb") as f:
                f.write(response.content)
            return target_path
        return None

    @staticmethod
    def run_p2rank(pdb_filename: str, project_id: str, log_callback):
        pdb_path = os.path.join(PROTEIN_STORAGE, pdb_filename)
        output_dir = os.path.join(PROTEIN_STORAGE, "p2rank_results", project_id)
        os.makedirs(output_dir, exist_ok=True)

        # Command: prank.bat predict -f <file> -o <out_dir>
        cmd = [P2RANK_EXEC, "predict", "-f", pdb_path, "-o", output_dir, "-v"]
        
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
            text=True, bufsize=1, universal_newlines=True
        )

        for line in process.stdout:
            log_callback(line.strip())
        
        process.wait()
        return output_dir

    @staticmethod
    def get_p2rank_results(pdb_filename: str, project_id: str):
        # P2Rank names output as: {pdb_filename}_predictions.csv
        csv_path = os.path.join(PROTEIN_STORAGE, "p2rank_results", project_id, f"{pdb_filename}_predictions.csv")
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            # Clean column names (strip spaces)
            df.columns = df.columns.str.strip()
            return df.to_dict(orient="records")
        return None