import os
import shutil
import threading
from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from config import IMG_DIR, BASE_DIR
from schemas import DiseaseQuery, TargetRequest, RetroRequest
from state import training_state

# Import Services
from services.llm_service import review_literature
from services.bio_service import RetroManager, fetch_chembl_data
from services.training_runner import run_training_task

app = FastAPI()
app.mount("/static", StaticFiles(directory=IMG_DIR), name="static")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/api/review")
async def api_review_literature(query: DiseaseQuery):
    try:
        return review_literature(query)
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/chembl/download")
async def api_get_chembl_data(req: TargetRequest):
    data = fetch_chembl_data(req.target_name)
    if not data:
        raise HTTPException(404, "Target not found or ChEMBL unavailable")
    return data

@app.post("/api/retrosynthesis")
async def run_retro(req: RetroRequest):
    try:
        routes = RetroManager.run_retrosynthesis(req.smiles)

        # Ensure routes is always a list
        if routes is None:
            routes = []
        if not isinstance(routes, list):
            routes = [routes]

        return {
            "is_solved": len(routes) > 0,
            "routes": routes
        }

    except Exception as e:
        print(f"❌ Retrosynthesis failed: {e}")
        return {
            "is_solved": False,
            "routes": [],
            "detail": str(e)
        }


@app.post("/api/train")
async def start_training(pdb_id: str = Form(...), file: UploadFile = File(...)):
    if training_state["is_running"]: 
        raise HTTPException(400, "Training is already in progress.")
    
    upload_path = os.path.join(BASE_DIR, "data", f"upload_{file.filename}")
    with open(upload_path, "wb+") as f:
        shutil.copyfileobj(file.file, f)
    
    # Run training in background thread
    threading.Thread(target=run_training_task, args=(upload_path, pdb_id)).start()
    return {"status": "started"}

@app.get("/api/training/status")
async def get_status():
    return training_state

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)