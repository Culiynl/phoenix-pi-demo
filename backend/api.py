import os
import sys
# Force flush stdout/stderr for Cloud Run logging
sys.stdout.reconfigure(line_buffering=True)
os.environ["CUDA_VISIBLE_DEVICES"] = "" # Disable CUDA for Cloud Run
print("PHOENIX-PI Backend Starting...", file=sys.stdout)
import shutil
import json
import io
import threading
from pydantic import BaseModel
from typing import Optional, List
from fastapi import FastAPI, HTTPException, UploadFile, File, Form, BackgroundTasks
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from config import IMG_DIR, BASE_DIR, DATA_DIR, PROJECTS_DIR 
from schemas import (
    DiseaseQuery, TargetRequest, RetroRequest, 
    DockingRequest, FeudalRequest, ParetoRequest
)
class MissionStartRequest(BaseModel):
    target_name: Optional[str] = None
    description: Optional[str] = None
    research_mode: Optional[str] = "hybrid" # arxiv, internal, or hybrid

from state import training_state
# Import Services
# Services will be imported lazily inside endpoints to improve Cloud Run startup time
from services.llm_service import review_literature

app = FastAPI()

# Add CORS middleware to allow frontend requests
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Background Orchestrator Tasks
orchestrator_instance = None
active_missions = {}

app.mount("/static", StaticFiles(directory=DATA_DIR), name="static")



@app.get("/api/research/arxiv")
async def api_search_arxiv(q: str):
    from services.llm_service import search_arxiv_papers, summarize_literature_gemini
    
    try:
        papers = search_arxiv_papers(q)
        # Optional: Auto-summarize if requested, or let frontend trigger it separately
        # summary = summarize_literature_gemini(papers)
        return {"papers": papers}
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/research/synthesize")
async def api_synthesize_research(payload: dict):
    from services.llm_service import summarize_literature_gemini
    papers = payload.get("papers", [])
    summary = summarize_literature_gemini(papers)
    return {"summary": summary}

@app.post("/api/chembl/download")
async def api_get_chembl_data(req: TargetRequest):
    from services.bio_service import fetch_chembl_data
    data = fetch_chembl_data(req.target_name)
    if not data:
        raise HTTPException(404, "Target not found or ChEMBL unavailable")
    return data

@app.post("/api/retrosynthesis")
async def run_retro(req: RetroRequest):
    # run_in_threadpool prevents the AI task from freezing the whole FastAPI server
    from starlette.concurrency import run_in_threadpool
    from services.bio_service import RetroManager
    try:
        results = await run_in_threadpool(RetroManager.run_retrosynthesis, req.smiles)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/train")
async def start_training(pdb_id: str = Form(...), file: UploadFile = File(...)):
    from services.training_runner import run_training_task
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

@app.post("/api/docking/run") 
async def api_run_docking(req: DockingRequest):
    from services.docking_service import DockingManager
    # Pass req.center_override into the manager
    result = DockingManager.run_docking(
        req.project_id, 
        req.smiles, 
        center_override=req.center_override
    )
    
    if not result["success"]:
        raise HTTPException(status_code=500, detail=result.get("error"))
    return result

@app.get("/api/feudal/status/{project_id}")
async def get_feudal_status(project_id: str):
    log_path = os.path.join(PROJECTS_DIR, project_id, "feudal_logs", "progress.json")
    hist_path = os.path.join(PROJECTS_DIR, project_id, "feudal_logs", f"run_{project_id}.json")
    
    data = {"logs": [], "history": []}
    
    if os.path.exists(log_path):
        with open(log_path, "r") as f: data["logs"] = json.load(f)
    
    if os.path.exists(hist_path):
        with open(hist_path, "r") as f: data["history"] = json.load(f)
            
    return data

@app.post("/api/feudal/optimize")
async def start_feudal_task(req: FeudalRequest, background_tasks: BackgroundTasks):
    from services.feudal_service import FeudalOptimizationManager
    feudal = FeudalOptimizationManager.get_instance()
    # We run this as a background task because it takes 1-2 minutes
    background_tasks.add_task(feudal.run_optimization, req.project_id, req.smiles, req.steps, req.batch_size)
    return {"status": "started", "message": "Optimization loop running in background"}

@app.post("/api/feudal/pareto-rank")
async def get_pareto(req: ParetoRequest):
    # Construct path to the optimization log generated by the feudal task
    log_path = os.path.join(PROJECTS_DIR, req.project_id, "feudal_logs", f"run_{req.project_id}.json")
    
    if not os.path.exists(log_path):
        raise HTTPException(status_code=404, detail="Optimization logs not found. Run optimization first.")
    
    try:
        with open(log_path, "r") as f:
            history = json.load(f)
        
        from services.analysis_service import calculate_pareto_winner
        # Calculate the winner using the weights provided by the UI
        winner = calculate_pareto_winner(history, req.weights)
        return winner
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Ranking failed: {str(e)}")
from services.pdb_service import PDBManager

@app.get("/api/pdb/fetch/{pdb_id}")
async def fetch_pdb(pdb_id: str):
    # result is now a dict: {"filename": "9BA8.cif", "format": "cif"}
    result = PDBManager.fetch_from_rcsb(pdb_id)
    
    if not result:
        raise HTTPException(404, "PDB or CIF not found in RCSB")
    
    # Return the data directly from the dictionary
    return {
        "status": "success", 
        "filename": result["filename"], 
        "format": result["format"]
    }

@app.post("/api/p2rank/run")
async def run_p2rank_api(pdb_filename: str, project_id: str, background_tasks: BackgroundTasks):
    # Pass the task to the background
    background_tasks.add_task(PDBManager.run_p2rank, pdb_filename, project_id)
    return {"status": "started"}
@app.get("/api/p2rank/status/{project_id}")
async def get_p2rank_status(project_id: str):
    log_path = os.path.join(DATA_DIR, "proteins", "p2rank_results", project_id, "p2_progress.json")
    if os.path.exists(log_path):
        with open(log_path, "r") as f:
            return json.load(f)
    return []
@app.get("/api/pdb/search")
async def search_pdb(q: str):
    results = PDBManager.search_rcsb(q)
    return results
@app.get("/api/p2rank/results")
async def get_p2rank_results(pdb_filename: str, project_id: str):
    results = PDBManager.get_p2rank_results(pdb_filename, project_id)
    if not results: raise HTTPException(404, "Results not found")
    return results
from rdkit import Chem
from rdkit.Chem import Draw
from fastapi.responses import Response
import io
@app.post("/api/mission/settings")
async def update_mission_settings(mission_id: str = Form(...), auto_correct: bool = Form(...)):
    """Toggles mission settings like auto-confirm."""
    from services.firebase_service import update_mission_field
    update_mission_field(mission_id, "auto_correct", auto_correct)
    return {"status": "updated", "auto_correct": auto_correct}

@app.get("/api/render")
async def render_molecule(smiles: str):
    if not smiles:
        raise HTTPException(400, "SMILES string required")
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        # Return a fallback "not found" image or empty
        mol = Chem.MolFromSmiles("C") # Just a methane placeholder or return error
        if not mol: raise HTTPException(400, "Invalid SMILES")
        
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return Response(content=buf.getvalue(), media_type="image/png")

@app.post("/api/mission/resume")
async def resume_mission(background_tasks: BackgroundTasks, mission_id: str = Form(...)):
    """Resumes an autonomous research mission."""
    from services.orchestrator_service import ResearchOrchestrator
    # Note: In a production app, we'd restore chat history from a database.
    # For this local demo, we'll try to find the active instance or create a new one.
    # We'll use a simple global registry for active missions.
    orchestrator = active_missions.get(mission_id)
    if not orchestrator:
        orchestrator = ResearchOrchestrator(mission_id)
        active_missions[mission_id] = orchestrator
    
    background_tasks.add_task(orchestrator.resume_mission)
    return {"status": "resuming", "mission_id": mission_id}

@app.post("/api/mission/start")
async def start_mission(background_tasks: BackgroundTasks, req: MissionStartRequest):
    """
    Starts an autonomous research mission for the given target.
    """
    target_name = req.target_name
    description = req.description
    from services.firebase_service import create_mission
    
    # Generate a clean mission ID
    from services.llm_service import generate_mission_title
    if (not target_name or target_name.strip() == "") and description:
        target_name = generate_mission_title(description)
    
    safe_target = "".join(c for c in target_name if c.isalnum()).lower()
    if not safe_target: safe_target = "mission"
    mission_id = f"mission_{safe_target}_{os.urandom(4).hex()}"
    
    # Initialize mission in Firestore
    try:
        create_mission(mission_id, target_name, research_mode=req.research_mode)
    except Exception as e:
        print(f"Firestore init failed: {e}")
    
    from services.orchestrator_service import ResearchOrchestrator
    orchestrator = ResearchOrchestrator(mission_id)
    active_missions[mission_id] = orchestrator
    background_tasks.add_task(orchestrator.start_mission, target_name, description)
    
    return {"status": "started", "mission_id": mission_id}

@app.post("/api/mission/{mission_id}/chat")
async def chat_mission(mission_id: str, payload: dict, background_tasks: BackgroundTasks):
    from services.firebase_service import update_mission_log
    from datetime import datetime
    
    message = payload.get("message", "").strip()
    if not message:
        return {"status": "ignored"}
        
    # Log user message to Firestore
    update_mission_log(mission_id, {
        "type": "user",
        "content": message,
        "timestamp": datetime.now().isoformat()
    })
    
    # Logic: If user says "confirm", "proceed", "yes", "go", etc., trigger resume
    trigger_words = ["confirm", "proceed", "yes", "go", "continue", "next"]
    if any(word in message.lower() for word in trigger_words):
        from services.orchestrator_service import ResearchOrchestrator
        orchestrator = active_missions.get(mission_id)
        if not orchestrator:
            orchestrator = ResearchOrchestrator(mission_id)
            active_missions[mission_id] = orchestrator
        
        background_tasks.add_task(orchestrator.resume_mission)
        return {"status": "resuming", "message": "Triggered resumption via chat."}

    return {"status": "received"}

@app.post("/api/mission/{mission_id}/upload")
async def upload_mission_file(mission_id: str, file: UploadFile = File(...), file_type: str = Form("image")):
    """Upload images or PDFs for multimodal analysis."""
    from config import PROJECTS_DIR
    
    upload_dir = os.path.join(PROJECTS_DIR, mission_id, "uploads")
    os.makedirs(upload_dir, exist_ok=True)
    
    file_path = os.path.join(upload_dir, file.filename)
    with open(file_path, "wb") as f:
        f.write(await file.read())
    
    return {
        "status": "uploaded",
        "filename": file.filename,
        "file_type": file_type,
        "path": file_path
    }

@app.post("/api/upload")
async def upload_file(file: UploadFile = File(...)):
    file_location = f"data/uploads/{file.filename}"
    os.makedirs("data/uploads", exist_ok=True)
    with open(file_location, "wb+") as file_object:
        file_object.write(file.file.read())
    return {"info": f"file '{file.filename}' saved at '{file_location}'"}

@app.get("/api/mission/{mission_id}/log")
async def get_mission_log(mission_id: str):
    # In a real implementation this would fetch from Firestore
    # For now we'll just return a mock or empty list if the instance isn't available
    return {"log": [], "status": "unknown"} # Placeholder until Firestore is connected


@app.get("/api/mission/download_csv/{mission_id}/{filename}")
async def download_mission_csv(mission_id: str, filename: str):
    csv_path = os.path.join(PROJECTS_DIR, mission_id, filename)
    if not os.path.exists(csv_path):
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(path=csv_path, filename=filename, media_type='text/csv')

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)