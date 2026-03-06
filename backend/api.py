import os
import sys
import asyncio
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
from fastapi import FastAPI, HTTPException, UploadFile, File, Form, BackgroundTasks, Request, Response
from fastapi.responses import FileResponse, JSONResponse
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
print("PHOENIX-PI Backend API v2.2 (Resumption Fixed) Initialized")

# Add CORS middleware to allow frontend requests
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:8000"],
    allow_methods=["*"],
    allow_headers=["*"],
    allow_credentials=True, 
)

# Deep Diagnostics for ImportError Troubleshooting
import importlib
import services.local_storage_service as ls
importlib.reload(ls)
print(f"--- DIAGNOSTICS ---")
print(f"LocalStorage File: {ls.__file__}")
import inspect
try:
    from services.local_storage_service import link_pdb_to_mission
    print("Function link_pdb_to_mission found in memory.")
except ImportError:
    print("Function link_pdb_to_mission NOT found in memory.")
print(f"-------------------")

# Background Orchestrator Tasks
# Format: { mission_id: {"orchestrator": obj, "task": asyncio.Task} }
active_missions = {}

from fastapi.responses import JSONResponse, Response

@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    import traceback
    traceback.print_exc()
    return JSONResponse(
        status_code=500,
        content={"message": f"Internal Server Error: {str(exc)}"},
        headers={
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "*",
            "Access-Control-Allow-Headers": "*",
        }
    )

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
async def fetch_pdb(pdb_id: str, mission_id: Optional[str] = None):
    # Pass mission_id if we want to save it directly to the project folder
    result = PDBManager.fetch_from_rcsb(pdb_id, mission_id=mission_id)
    
    if not result:
        raise HTTPException(404, "PDB or CIF not found in RCSB")
    
    return {
        "status": "success", 
        "filename": result["filename"],
        "path": result.get("path")
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
import io
import asyncio
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdChemReactions
from concurrent.futures import ThreadPoolExecutor
@app.post("/api/mission/settings")
async def update_mission_settings(mission_id: str = Form(...), auto_correct: bool = Form(...)):
    """Toggles mission settings like auto-confirm."""
    from services.local_storage_service import update_mission_field
    from state import log_msg
    log_msg(f"   [API] Updating mission settings for {mission_id}: auto_correct={auto_correct} (type: {type(auto_correct)})")
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

@app.get("/api/render/reaction")
async def render_reaction(reactants: str, product: str):
    """Renders a chemical reaction scheme: Reactant1 + Reactant2 -> Product."""
    try:
        from rdkit.Chem import rdChemReactions
        reactant_smis = reactants.split(",")
        
        rxn = rdChemReactions.ChemicalReaction()
        for s in reactant_smis:
            m = Chem.MolFromSmiles(s.strip())
            if m: rxn.AddReactant(m)
        
        p = Chem.MolFromSmiles(product.strip())
        if p: rxn.AddProduct(p)
        
        # Draw the reaction
        d2d = Draw.MolDraw2DCairo(800, 300)
        d2d.DrawReaction(rxn)
        d2d.FinishDrawing()
        
        buf = io.BytesIO()
        buf.write(d2d.GetDrawingText())
        return Response(content=buf.getvalue(), media_type="image/png")
    except Exception as e:
        raise HTTPException(500, str(e))

@app.post("/api/mission/resume")
async def resume_mission(background_tasks: BackgroundTasks, mission_id: str = Form(...)):
    """Resumes an autonomous research mission."""
    from services.orchestrator_service import ResearchOrchestrator
    # Note: In a production app, we'd restore chat history from a database.
    # For this local demo, we'll try to find the active instance or create a new one.
    # We'll use a simple global registry for active missions.
    orchestrator_data = active_missions.get(mission_id)
    if not orchestrator_data:
        orchestrator = ResearchOrchestrator(mission_id)
        active_missions[mission_id] = {"orchestrator": orchestrator, "task": None}
    else:
        orchestrator = orchestrator_data["orchestrator"]
    
    # Cancel existing task if any
    if active_missions[mission_id]["task"] and not active_missions[mission_id]["task"].done():
        active_missions[mission_id]["task"].cancel()
        
    task = asyncio.create_task(orchestrator.resume_mission())
    active_missions[mission_id]["task"] = task
    return {"status": "resuming", "mission_id": mission_id}

@app.post("/api/mission/start")
async def start_mission(background_tasks: BackgroundTasks, req: MissionStartRequest):
    """
    Starts an autonomous research mission for the given target.
    """
    target_name = req.target_name
    description = req.description
    from services.local_storage_service import create_mission
    
    # Generate a clean mission ID
    from services.llm_service import generate_mission_title
    if (not target_name or target_name.strip() == "") and description:
        target_name = generate_mission_title(description)
    
    safe_target = "".join(c for c in target_name if c.isalnum()).lower()
    if not safe_target: safe_target = "mission"
    mission_id = f"mission_{safe_target}_{os.urandom(4).hex()}"
    
    # Initialize mission in Local Storage
    try:
        create_mission(mission_id, target_name, research_mode=req.research_mode, description=description)
    except Exception as e:
        print(f"Local storage init failed: {e}")
    
    from services.orchestrator_service import ResearchOrchestrator
    orchestrator = ResearchOrchestrator(mission_id)
    
    # Stop all other active missions if starting a new one (as requested)
    for mid, data in active_missions.items():
        if data["task"] and not data["task"].done():
            print(f"Cancelling Mission: {mid}")
            data["task"].cancel()

    task = asyncio.create_task(orchestrator.start_mission(target_name, description))
    active_missions[mission_id] = {"orchestrator": orchestrator, "task": task}
    
    return {"status": "started", "mission_id": mission_id}

@app.post("/api/mission/{mission_id}/chat")
async def chat_mission(mission_id: str, payload: dict, background_tasks: BackgroundTasks):
    from services.local_storage_service import update_mission_log
    from datetime import datetime
    
    try:
        message = payload.get("message", "").strip()
        if not message:
            return JSONResponse(content={"status": "ignored"}, headers={"Access-Control-Allow-Origin": "*"})
            
        print(f"[API] Chat received for {mission_id}: {message}")
        
        # Log user message to Local Storage
        update_mission_log(mission_id, {
            "type": "user",
            "content": message,
            "timestamp": datetime.now().isoformat()
        })
        
        # Logic: If user says "confirm", "proceed", "yes", "go", etc., trigger resume
        trigger_words = ["confirm", "proceed", "yes", "go", "continue", "next"]
        if any(word in message.lower() for word in trigger_words):
            print(f"[API] Triggering resumption for {mission_id}")
            from services.orchestrator_service import ResearchOrchestrator
            
            orchestrator_data = active_missions.get(mission_id)
            if not orchestrator_data:
                print(f"[API] Initializing new orchestrator for {mission_id}")
                orchestrator = ResearchOrchestrator(mission_id)
                active_missions[mission_id] = {"orchestrator": orchestrator, "task": None}
            else:
                orchestrator = orchestrator_data["orchestrator"]
            
            if active_missions[mission_id]["task"] and not active_missions[mission_id]["task"].done():
                print(f"[API] Cancelling old task for {mission_id}")
                active_missions[mission_id]["task"].cancel()
                
            # Use asyncio.create_task to ensure it runs in the same event loop
            task = asyncio.create_task(orchestrator.resume_mission())
            active_missions[mission_id]["task"] = task
            
            return JSONResponse(
                content={"status": "resuming", "message": "Triggered resumption via chat."},
                headers={"Access-Control-Allow-Origin": "*"}
            )

        return JSONResponse(content={"status": "received"}, headers={"Access-Control-Allow-Origin": "*"})
    except Exception as e:
        import traceback
        traceback.print_exc()
        return JSONResponse(
            status_code=500,
            content={"status": "error", "message": str(e)},
            headers={"Access-Control-Allow-Origin": "*"}
        )

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
    from services.local_storage_service import get_mission
    mission = get_mission(mission_id)
    if mission:
        return {"log": mission.get("thought_log", []), "status": mission.get("status", "unknown")}
    return {"log": [], "status": "unknown"}


@app.get("/api/mission/download_csv/{mission_id}/{filename}")
async def download_mission_csv(mission_id: str, filename: str):
    csv_path = os.path.join(PROJECTS_DIR, mission_id, filename)
    if not os.path.exists(csv_path):
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(path=csv_path, filename=filename, media_type='text/csv')

@app.get("/api/mission/{mission_id}")
async def get_mission_data(mission_id: str):
    from services.local_storage_service import get_mission
    mission = get_mission(mission_id)
    if not mission:
        raise HTTPException(404, "Mission not found")
    return mission

@app.get("/api/mission/export/{mission_id}")
async def export_mission(mission_id: str):
    from services.local_storage_service import get_mission, _clean_serializable
    mission = get_mission(mission_id)
    if not mission:
        raise HTTPException(404, "Mission not found")
        
    filename = f"{mission_id}.json"
    # Create a temporary file or stream
    import io
    import json
    clean_mission = _clean_serializable(mission)
    stream = io.StringIO(json.dumps(clean_mission, indent=2))
    return Response(content=stream.getvalue(), media_type="application/json", headers={"Content-Disposition": f"attachment; filename={filename}"})

@app.post("/api/mission/import")
async def import_mission(file: UploadFile = File(...)):
    from services.local_storage_service import _save_mission
    try:
        content = await file.read()
        data = json.loads(content)
        # Validate minimal fields
        if "id" not in data:
            raise HTTPException(400, "Invalid mission file (missing id)")
            
        _save_mission(data["id"], data)
        return {"status": "imported", "id": data["id"]}
    except Exception as e:
        raise HTTPException(500, f"Import failed: {str(e)}")

@app.get("/api/projects")
async def list_projects():
    from services.local_storage_service import get_all_projects
    return get_all_projects()

@app.delete("/api/projects/{project_id}")
async def delete_project_api(project_id: str):
    from services.local_storage_service import delete_project
    success = delete_project(project_id)
    if not success:
        raise HTTPException(status_code=404, detail="Project not found")
    return {"status": "deleted", "id": project_id}

@app.get("/api/mission/{mission_id}/download/pdb")
async def download_mission_pdb(mission_id: str):
    """Serves the docked complex PDB file if it exists."""
    pdb_path = os.path.join(PROJECTS_DIR, mission_id, "analysis", "complex.pdb")
    if not os.path.exists(pdb_path):
        # Fallback to receptor if complex isn't generated yet
        pdb_path = os.path.join(PROJECTS_DIR, mission_id, "docking", "best_pose.pdbqt") # or similar
        if not os.path.exists(pdb_path):
            raise HTTPException(status_code=404, detail="Docked complex PDB not found. Run interaction analysis first.")
    
    return FileResponse(path=pdb_path, filename=f"{mission_id}_docked_complex.pdb", media_type='application/octet-stream')

@app.get("/api/mission/{mission_id}/export/pdf")
async def export_mission_pdf(mission_id: str):
    """Generates and serves a comprehensive mission report PDF."""
    from services.local_storage_service import get_mission
    from reportlab.lib.pagesizes import LETTER
    from reportlab.pdfgen import canvas
    from reportlab.lib import colors
    from reportlab.lib.units import inch
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, HRFlowable
    from rdkit import Chem
    from rdkit.Chem import Draw
    import os
    
    mission = get_mission(mission_id)
    if not mission:
        raise HTTPException(404, "Mission not found")
        
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=LETTER, topMargin=0.5*inch)
    styles = getSampleStyleSheet()
    
    # Custom Styles
    accent_color = colors.HexColor("#7f1d1d")  # Dark Red
    subtitle_style = ParagraphStyle(
        'PhoenixSubtitle',
        parent=styles['Normal'],
        fontSize=9,
        textColor=colors.HexColor("#000000"),
        spaceAfter=4,
        letterSpacing=2
    )
    label_style = ParagraphStyle(
        'PhoenixLabel',
        parent=styles['Normal'],
        fontSize=10,
        textColor=colors.HexColor("#bbbbbb"),
        spaceAfter=3
    )
    title_style = ParagraphStyle(
        'PhoenixTitle',
        parent=styles['Title'],
        fontSize=30,
        leading=34,
        textColor=accent_color,
        spaceAfter=4,
        fontName='Helvetica-Bold'
    )
    heading_style = ParagraphStyle(
        'PhoenixHeading',
        parent=styles['Heading2'],
        fontSize=13,
        fontName='Helvetica-Bold',
        textColor=accent_color,
        spaceBefore=20,
        spaceAfter=8,
        borderPadding=(0, 0, 4, 0)
    )
    
    story = []

    logo_path = os.path.join(os.path.dirname(BASE_DIR), "PHOENIX_LOGO.png")

    # ── HEADER BLOCK ─────────────────────────────────────────────────────────
    # Logo on the left, stacked title on the right
    if os.path.exists(logo_path):
        logo_img = Image(logo_path, width=1.0*inch, height=1.0*inch)
    else:
        logo_img = Paragraph("☽", title_style)

    title_block = [
        # Paragraph("PHOENIX", subtitle_style),
        Paragraph("FINAL REPORT", title_style),
        # Paragraph(f"Mission: {mission.get('target_name', 'Unknown Target')}", subtitle_style),
    ]

    header_table = Table(
        [[logo_img, title_block]],
        colWidths=[1.2*inch, 5.3*inch]
    )
    header_table.setStyle(TableStyle([
        ('VALIGN',  (0, 0), (-1, -1), 'MIDDLE'),
        ('ALIGN',   (0, 0), (0, 0),  'CENTER'),
        ('ALIGN',   (1, 0), (1, 0),  'LEFT'),
        ('LEFTPADDING',  (1, 0), (1, 0), 10),
        ('RIGHTPADDING', (0, 0), (-1,-1), 0),
    ]))
    story.append(header_table)
    story.append(Spacer(1, 6))
    story.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor("#7f1d1d"), spaceBefore=2, spaceAfter=4))

    # Meta line
    from datetime import datetime
    generated_at = datetime.now().strftime("%B %d, %Y  %H:%M UTC")
    meta_data = [
        [Paragraph(f"<b>Target:</b>  {mission.get('target_name', 'N/A')}", label_style),
         Paragraph(f"<b>Mission ID:</b>  {mission_id}", label_style),
         Paragraph(f"<b>Generated:</b>  {generated_at}", label_style)]
    ]
    meta_table = Table(meta_data, colWidths=[2.2*inch, 2.5*inch, 2*inch])
    meta_table.setStyle(TableStyle([
        ('ALIGN', (0,0), (-1,-1), 'LEFT'),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('FONTSIZE', (0,0), (-1,-1), 9),
    ]))
    story.append(meta_table)
    story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor("#333333"), spaceBefore=6, spaceAfter=16))

    
    # Summary (AI Generated)
    summary = mission.get("final_summary") or mission.get("final_report") or "Mission summary pending completion."
    story.append(Paragraph("Intelligence Executive Summary", heading_style))
    # 1. ADD THESE IMPORTS at the top of your function (or file)
    import re
    import html

    # 2. DEFINE A BULLET STYLE (Add this right after your 'heading_style' definition)
    bullet_style = ParagraphStyle(
        'PhoenixBullet',
        parent=styles['Normal'],
        leftIndent=20,
        firstLineIndent=-10,
        spaceBefore=3,
        spaceAfter=3
    )

    for line in summary.split('\n'):
        line = line.strip()
        if not line:
            continue

        current_style = styles['Normal']
        is_bullet = False

        # --- Handle Markdown Block Elements ---
        if line.startswith('### '):
            line = line[4:]
            current_style = styles['Heading3']
        elif line.startswith('#### '):
            line = line[5:]
            current_style = styles['Heading4']
        elif line.startswith('* '):
            line = line[2:]
            current_style = bullet_style
            is_bullet = True

        # --- Escape Special Characters (Crucial for SMILES and XML safety) ---
        # This prevents characters like < or > in SMILES from crashing ReportLab
        line = html.escape(line)

        # --- Handle Markdown Inline Elements (Regex) ---
        # Bold: **text** -> <b>text</b>
        line = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', line)
        
        # Italics: *text* -> <i>text</i>
        line = re.sub(r'(?<!\w)\*(.*?)\*(?!\w)', r'<i>\1</i>', line)
        
        # Inline Code/SMILES: `text` -> <font name="Courier">text</font>
        line = re.sub(r'`(.*?)`', r'<font name="Courier" color="blue">\1</font>', line)

        # Add bullet symbol if needed
        if is_bullet:
            line = f"&bull; {line}"

        # Append to story
        story.append(Paragraph(line, current_style))
        
        # Add extra spacing after headings
        if current_style.name.startswith('Heading'):
            story.append(Spacer(1, 4))
    
    story.append(Spacer(1, 15))
    
    # 2D Structures: Before and After
    story.append(Paragraph("Molecular Optimization Profile", heading_style))
    
    mol_images = []
    
    # Initial Molecule (from Pareto seed if exists)
    # The 'feudal_logs/run_{project_id}.json' history[0] contains the seed
    seed_smi = None
    log_path = os.path.join(PROJECTS_DIR, mission_id, "feudal_logs", f"run_{mission_id}.json")
    if os.path.exists(log_path):
        try:
            with open(log_path, 'r') as f:
                history = json.load(f)
                if history: seed_smi = history[0].get('smiles')
        except: pass

    # Lead Candidate
    profile = mission.get("best_candidate_profile", {})
    lead_smi = profile.get('smiles')
    
    def get_mol_img(smi, label):
        if not smi: return Paragraph("N/A", styles['Normal'])
        mol = Chem.MolFromSmiles(smi)
        if not mol: return Paragraph("Invalid SMILES", styles['Normal'])
        img_buf = io.BytesIO()
        canv = Draw.MolToImage(mol, size=(300, 300))
        canv.save(img_buf, format="PNG")
        img_buf.seek(0)
        return [Image(img_buf, width=2.5*inch, height=2.5*inch), Paragraph(label, styles['Normal'])]

    comparison_data = [
        [get_mol_img(seed_smi, "<b>Initial Scaffold</b>"), get_mol_img(lead_smi, "<b>Optimized Lead</b>")]
    ]
    
    comparison_table = Table(comparison_data, colWidths=[3.2*inch, 3.2*inch])
    comparison_table.setStyle(TableStyle([
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ]))
    story.append(comparison_table)
    
    # Property Table
    if profile:
        story.append(Paragraph("Physicochemical Property Analysis", heading_style))
        props = profile.get("properties", {})
        admet = profile.get("admet", {})
        
        table_data = [["Property", "Value", "Metric"]]
        for k, v in props.items():
            val_str = f"{v:.2f}" if isinstance(v, (float, int)) else str(v)
            table_data.append([k.upper(), val_str, "Standard"])
        
        for k, v in admet.items():
            val_str = f"{v:.2f}" if isinstance(v, (float, int)) else str(v)
            table_data.append([k.replace("_", " ").title(), val_str, "ADMET"])
            
        t = Table(table_data, colWidths=[2.5*inch, 1.5*inch, 1.5*inch])
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#7f1d1d")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
        ]))
        story.append(t)

    # Manager Goals Table
    # We extract this from the Feudal Optimization history
    if os.path.exists(log_path):
        try:
            with open(log_path, 'r') as f:
                history = json.load(f)
                if history:
                    story.append(Paragraph("Manager Strategy Alignment", heading_style))
                    goal_data = [["Step", "Manager Strategy", "Score (kcal/mol)", "SMILES (Truncated)"]]

                    prev_score = None
                    for i, entry in enumerate(history[:10]):
                        step = entry.get("step", i)
                        score = entry.get("score", 0)
                        status = entry.get("status", "Confirmed")

                        # Derive human-readable goal from score delta
                        if prev_score is None:
                            goal_label = "Identify seed scaffold"
                        else:
                            delta = score - prev_score
                            if delta < -0.5:
                                goal_label = f"Maximize binding ↑ (Δ {delta:.2f})"
                            elif delta < 0:
                                goal_label = f"Improve affinity (Δ {delta:.2f})"
                            elif delta == 0:
                                goal_label = "Maintain scaffold integrity"
                            else:
                                goal_label = f"Explore scaffold diversity (Δ +{delta:.2f})"
                        prev_score = score

                        goal_data.append([
                            str(step),
                            goal_label,
                            f"{score:.2f}" if score else "N/A",
                            (entry.get("smiles", "")[:25] + "...") if entry.get("smiles") else "N/A"
                        ])

                    gt = Table(goal_data, colWidths=[0.5*inch, 2.5*inch, 1.1*inch, 2.4*inch])
                    gt.setStyle(TableStyle([
                        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#7f1d1d")),
                        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                        ('FONTSIZE', (0, 0), (-1, -1), 8),
                        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor("#f9f9f9")]),
                        ('ALIGN', (2, 0), (2, -1), 'CENTER'),
                    ]))
                    story.append(gt)
        except: pass


    # Retrosynthesis Pathway
    retro_data = mission.get("retro_data", {})
    if retro_data and "routes" in retro_data:
        routes = retro_data["routes"]
        if routes:
            story.append(Paragraph("Retrosynthetic Pathway (AiZynthFinder)", heading_style))
            best_route = routes[0] # Take path 1
            if best_route.get("image"):
                img_p = os.path.join(DATA_DIR, best_route["image"])
                if os.path.exists(img_p):
                    story.append(Image(img_p, width=6.5*inch, height=3*inch, kind='proportional'))
                    # story.append(Paragraph(f"Confidence Score: {best_route.get('score', 0)*100:.1f}%", styles['Italic']))

    # Automated Analysis (Interactions)
    interactions = profile.get("interactions", {})
    if interactions:
        story.append(Paragraph("Automated Structural Interaction Analysis", heading_style))
        plip = interactions.get("plip", [])
        if plip:
            it_data = [["Type", "Residue", "Distance (A)"]]
            for it in plip:
                it_data.append([it.get("type"), it.get("residue"), str(it.get("distance"))])
            
            it_table = Table(it_data, colWidths=[1.5*inch, 2.5*inch, 1.5*inch])
            it_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#555555")),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
            ]))
            story.append(it_table)
            
        rationale = profile.get("rationale")
        if rationale:
            story.append(Spacer(1, 10))
            story.append(Paragraph("<b>Interaction Rationale:</b>", styles['Normal']))
            story.append(Paragraph(rationale, styles['Italic']))

    doc.build(story)
    buffer.seek(0)
    
    return Response(
        content=buffer.getvalue(),
        media_type="application/pdf",
        headers={"Content-Disposition": f"attachment; filename={mission_id}_report.pdf"}
    )

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)