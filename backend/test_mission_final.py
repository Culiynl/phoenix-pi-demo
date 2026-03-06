import requests
import json
import os

BASE_URL = "http://localhost:8000"
MISSION_ID = "mission_dynamicalzheimers_1cf87696"

def test_endpoints():
    print(f"Testing Mission Data for {MISSION_ID}...")
    resp = requests.get(f"{BASE_URL}/api/mission/{MISSION_ID}")
    if resp.status_code == 200:
        data = resp.json()
        print(f"  - final_summary present: {'final_summary' in data}")
        print(f"  - final_report present: {'final_report' in data}")
        print(f"  - retro_data present: {'retro_data' in data}")
        print(f"  - best_candidate_profile.interactions present: {('interactions' in data.get('best_candidate_profile', {}))}")
    else:
        print(f"  - Failed to fetch mission data: {resp.status_code}")

    print("\nTesting PDF Export...")
    resp = requests.get(f"{BASE_URL}/api/mission/{MISSION_ID}/export/pdf")
    if resp.status_code == 200:
        print(f"  - PDF Export Success (Content Length: {len(resp.content)})")
    else:
        print(f"  - PDF Export Failed: {resp.status_code}")

    print("\nTesting PDB Download...")
    resp = requests.get(f"{BASE_URL}/api/mission/{MISSION_ID}/download/pdb")
    if resp.status_code == 200:
        print(f"  - PDB Download Success (Content Length: {len(resp.content)})")
    elif resp.status_code == 404:
        print("  - PDB Download Failed (404 as expected if complex.pdb missing)")
    else:
        print(f"  - PDB Download Error: {resp.status_code}")

if __name__ == "__main__":
    test_endpoints()
