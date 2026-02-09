import os
import subprocess
import requests
import pandas as pd
import json
from config import PROTEIN_STORAGE, P2RANK_EXEC

class PDBManager:
    # services/pdb_service.py update
    @staticmethod
    def fetch_from_rcsb(pdb_id: str):
        pdb_id = pdb_id.upper() # RCSB likes Uppercase
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        # Ensure this directory exists in config.py: PROTEIN_STORAGE = data/proteins
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
        # output_dir from run_p2rank was: PROTEIN_STORAGE / "p2rank_results" / project_id
        base_dir = os.path.join(PROTEIN_STORAGE, "p2rank_results", project_id)
        
        # P2Rank naming possibilities:
        # 1. 6OD6.pdb_predictions.csv
        # 2. 6OD6_predictions.csv
        possible_names = [
            f"{pdb_filename}_predictions.csv",
            f"{pdb_filename.replace('.pdb', '')}_predictions.csv"
        ]

        csv_path = None
        for name in possible_names:
            target = os.path.join(base_dir, name)
            if os.path.exists(target):
                csv_path = target
                break

        if csv_path:
            try:
                print(f"[P2Rank] Loading results from: {csv_path}")
                df = pd.read_csv(csv_path)
                # P2Rank CSVs have spaces in headers like "  rank"
                df.columns = df.columns.str.strip()
                
                # Convert numeric columns safely
                return df.to_dict(orient="records")
            except Exception as e:
                print(f"[P2Rank] Error reading CSV: {e}")
                return None
        
        return None
    
    @staticmethod
    def search_rcsb(query_text: str):
        """Searches RCSB PDB for keyword and returns metadata."""
        # RCSB Search API v2 URL
        search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        
        # Search Query Structure
        query = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {"value": query_text}
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 25},
                "results_content_type": ["experimental"]
            }
        }
        
        response = requests.post(search_url, json=query)
        if response.status_code != 200:
            return []
            
        results = response.json().get("result_set", [])
        pdb_ids = [res["identifier"] for res in results]
        
        if not pdb_ids:
            return []

        # Fetch Detailed Metadata for these IDs
        data_url = f"https://data.rcsb.org/rest/v1/core/entry/"
        metadata_list = []
        
        for p_id in pdb_ids:
            meta_res = requests.get(data_url + p_id)
            if meta_res.status_code == 200:
                d = meta_res.json()
                metadata_list.append({
                    "id": p_id,
                    "title": d.get("struct", {}).get("title", "No title"),
                    "resolution": d.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                    # Fixed organism path
                    "organism": d.get("rcsb_entity_source_organism", [{}])[0].get("rcsb_taxonomy_lineage", [{}])[-1].get("name", "Unknown"),
                    "date": d.get("rcsb_accession_info", {}).get("deposit_date", "").split("T")[0],
                    "link": f"https://www.rcsb.org/structure/{p_id}"
                })
        
        return metadata_list