import sys
import os
import asyncio

# Add backend to path
sys.path.append(os.path.abspath("."))

from services.bio_service import RetroManager
from services.local_storage_service import update_mission_field

async def main():
    mission_id = "mission_dynamicalzheimers_1cf87696"
    smiles = "O=C(O)c1coc2cc(OCc3ccc(F)cc3)ccc2c1=O"
    
    print(f"Updating retrosynthesis for {mission_id}...")
    # Since run_retrosynthesis is a static method and not actually async in my previous view (wait, let me check)
    # In orchestrator_service.py: results = await run_in_threadpool(RetroManager.run_retrosynthesis, smiles)
    # So it's a sync method.
    
    results = RetroManager.run_retrosynthesis(smiles)
    if "error" in results:
        print(f"Error: {results['error']}")
        return

    update_mission_field(mission_id, "retro_data", results)
    print("Successfully updated mission data.")

if __name__ == "__main__":
    asyncio.run(main())
