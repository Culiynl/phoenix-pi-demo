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

    import glob

@staticmethod
def get_p2rank_results(pdb_filename: str, project_id: str):
    # The base directory we told P2Rank to save to
    base_dir = os.path.join(PROTEIN_STORAGE, "p2rank_results", project_id)
    
    if not os.path.exists(base_dir):
        print(f"[P2Rank] Directory does not exist yet: {base_dir}")
        return None

    # Get the ID (e.g., 6OD6)
    base_stem = pdb_filename.split('.')[0].upper()
    
    # RECURSIVE SEARCH: Look in base_dir and ALL subdirectories (like /visualizations/)
    # The '**' means search all subfolders
    search_pattern = os.path.join(base_dir, "**", f"*{base_stem}*predictions.csv")
    found_files = glob.glob(search_pattern, recursive=True)

    if not found_files:
        # Log exactly where we looked to help debugging
        print(f"[P2Rank] Still waiting for CSV. Searched: {search_pattern}")
        return None

    # Use the first valid predictions file found
    csv_path = found_files[0]
    
    try:
        print(f"[P2Rank] File found! Loading: {csv_path}")
        df = pd.read_csv(csv_path)
        df.columns = df.columns.str.strip()
        
        # Clean data for JSON
        results = []
        for _, row in df.iterrows():
            results.append({
                "rank": int(row['rank']),
                "score": float(row['score']),
                "probability": float(row['probability']),
                "center_x": float(row['center_x']),
                "center_y": float(row['center_y']),
                "center_z": float(row['center_z']),
                "residue_ids": str(row['residue_ids'])
            })
        return results
    except Exception as e:
        print(f"[P2Rank] Error parsing CSV: {e}")
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