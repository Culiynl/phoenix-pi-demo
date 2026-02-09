import glob
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
        pdb_id = pdb_id.upper()
        formats = ['pdb', 'cif']
        
        for fmt in formats:
            url = f"https://files.rcsb.org/download/{pdb_id}.{fmt}"
            target_path = os.path.join(PROTEIN_STORAGE, f"{pdb_id}.{fmt}")
            
            response = requests.get(url)
            if response.status_code == 200:
                with open(target_path, "wb") as f:
                    f.write(response.content)
                print(f"[PDB] Successfully downloaded {pdb_id} as {fmt}")
                return {"filename": f"{pdb_id}.{fmt}", "format": fmt}
                
        print(f"[PDB] Error: {pdb_id} not found in PDB or CIF format.")
        return None

    # services/pdb_service.py

    @staticmethod
    def run_p2rank(pdb_filename: str, project_id: str):
        # Safety Check
        if not project_id or project_id == "undefined":
            print("[P2Rank] CRITICAL ERROR: project_id is undefined")
            return

        pdb_path = os.path.normpath(os.path.join(PROTEIN_STORAGE, pdb_filename))
        output_dir = os.path.normpath(os.path.join(PROTEIN_STORAGE, "p2rank_results", project_id))
        log_file = os.path.join(output_dir, "p2_progress.json")
        os.makedirs(output_dir, exist_ok=True)

        # Initialize log file
        with open(log_file, "w") as f:
            json.dump([{"time": "", "msg": f"P2Rank started for {pdb_filename}...", "type": "info"}], f)

        def write_ui_log(message):
            try:
                with open(log_file, "r") as f: logs = json.load(f)
                logs.append({"time": "", "msg": message, "type": "info"})
                with open(log_file, "w") as f: json.dump(logs, f)
            except: pass

        # Use 'call' for batch files or ensure CMD environment
        cmd = [P2RANK_EXEC, "predict", "-f", pdb_path, "-o", output_dir, "-v"]
        P2RANK_DIR = os.path.dirname(P2RANK_EXEC)
        # We use subprocess.PIPE and bufsize=1 (line buffered)
        process = subprocess.Popen(
            cmd, 
            cwd=P2RANK_DIR,
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True, 
            bufsize=1, 
            universal_newlines=True
        )

        # Use iter to read lines as they appear (prevents buffering freeze)
        for line in iter(process.stdout.readline, ''):
            clean_line = line.strip()
            if clean_line:
                print(f"[P2Rank CLI]: {clean_line}")
                write_ui_log(clean_line)
        
        process.stdout.close()
        process.wait()

        if process.returncode == 0:
            write_ui_log("P2Rank analysis complete!")
        else:
            write_ui_log(f"P2Rank process exited with error code {process.returncode}")

    import glob

    @staticmethod
    def get_p2rank_results(pdb_filename: str, project_id: str):
        # 1. Define the base project results folder
        base_dir = os.path.join(PROTEIN_STORAGE, "p2rank_results", project_id)
        
        if not os.path.exists(base_dir):
            return None

        # 2. Get the stem (e.g., '6OD6' from '6OD6.pdb')
        protein_stem = pdb_filename.split('.')[0]
        
        # 3. ROOT PROBLEM FIX: Recursive search
        # We look inside project_id/**/ (any subfolder) for a CSV containing the protein stem
        search_pattern = os.path.join(base_dir, "**", f"*{protein_stem}*predictions.csv")
        print(f"[P2Rank API] Looking for nested results: {search_pattern}")
        
        # recursive=True allows glob to look inside those 'predict_6OD6' subfolders
        found_files = glob.glob(search_pattern, recursive=True)

        if not found_files:
            print(f"[P2Rank API] CSV not found yet in subfolders of {base_dir}")
            return None

        # 4. Pick the newest file if multiple exist
        csv_path = max(found_files, key=os.path.getmtime)
        
        try:
            print(f"[P2Rank API] Found file at: {csv_path}")
            df = pd.read_csv(csv_path)
            # Standardize column names (remove leading/trailing spaces)
            df.columns = df.columns.str.strip()
            
            # 5. Extract and format for UI
            pockets = []
            for _, row in df.iterrows():
                pockets.append({
                    "rank": int(row['rank']),
                    "score": float(row['score']),
                    "probability": float(row['probability']),
                    "center_x": float(row['center_x']),
                    "center_y": float(row['center_y']),
                    "center_z": float(row['center_z']),
                    "residue_ids": str(row['residue_ids'])
                })
            return pockets
        except Exception as e:
            print(f"[P2Rank API] Critical Error parsing CSV: {e}")
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