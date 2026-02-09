"""
Heavy Compute Service
Runs on Compute Engine VM to handle large ML models (AiZynthFinder, P2Rank, etc.)
"""

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import os

app = FastAPI(title="PHOENIX-PI Heavy Compute Service")


class RetroRequest(BaseModel):
    smiles: str


class DockingRequest(BaseModel):
    pdb_id: str
    ligand_smiles: str


@app.post("/retro")
async def retrosynthesis(request: RetroRequest):
    """Run AiZynthFinder retrosynthesis analysis."""
    try:
        from services.bio_service import RetroManager
        result = RetroManager.run_retrosynthesis(request.smiles)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/docking")
async def molecular_docking(request: DockingRequest):
    """Run P2Rank molecular docking."""
    try:
        from services.bio_service import DockingManager
        result = DockingManager.run_docking(request.pdb_id, request.ligand_smiles)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "service": "heavy-compute"}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)
