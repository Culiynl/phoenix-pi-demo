import os
import json
import time
import shutil
import threading
from config import PROJECTS_DIR

# Global lock to prevent concurrent file writes to mission_data.json
# on Windows which is very strict about file access.
mission_lock = threading.Lock()

print("[LocalStorage] Module loaded: v2.3 (link_pdb_to_mission check)")

def _get_mission_path(mission_id):
    return os.path.join(PROJECTS_DIR, mission_id, "mission_data.json")

def _ensure_project_dir(mission_id):
    path = os.path.join(PROJECTS_DIR, mission_id)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def get_mission(mission_id):
    path = _get_mission_path(mission_id)
    if not os.path.exists(path):
        return None
    try:
        with mission_lock:
            with open(path, 'r', encoding='utf-8') as f:
                return json.load(f)
    except Exception as e:
        print(f"[LocalStorage] Error reading mission {mission_id}: {e}")
        return None

def create_mission(mission_id, target_name, research_mode="hybrid", description=""):
    _ensure_project_dir(mission_id)
    path = _get_mission_path(mission_id)
    
    # Check if already exists? usually api.py handles generation of ID
    
    data = {
        "id": mission_id,
        "target_name": target_name,
        "description": description,
        "research_mode": research_mode,
        "status": "started",
        "created_at": int(time.time()), # Store as timestamp
        "date": time.strftime("%Y-%m-%d"), # For UI display
        "owner": "Guest Researcher",     # Guest mode default
        "thought_log": [],
        "tool_history": [],
        "current_step": 1,
        "is_thinking": False,
        "awaiting_confirmation": False,
        "research_summary": "",
        "internal_research_summary": "",
        "internal_research_papers": [],
        "best_molecule": None # For the new "Best Molecule" feature
    }
    
    clean_data = _clean_serializable(data)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(clean_data, f, indent=2)
    return clean_data

def update_mission_log(mission_id, log_entry):
    mission = get_mission(mission_id)
    if mission:
        if "thought_log" not in mission:
            mission["thought_log"] = []
        mission["thought_log"].append(log_entry)
        
        # Save back
        _save_mission(mission_id, mission)

def update_mission_field(mission_id, field, value):
    mission = get_mission(mission_id)
    if mission:
        mission[field] = value
        _save_mission(mission_id, mission)

def _clean_serializable(obj):
    """Recursively ensures object is JSON serializable, handling Protobuf/MapComposite."""
    # Basic primitives
    if isinstance(obj, (str, int, float, bool, type(None))):
        return obj
    
    # Common collections
    if isinstance(obj, dict):
        return {str(k): _clean_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_clean_serializable(item) for item in obj]
    
    # Protobuf / Google AI / Special library objects (caught by name)
    obj_type = type(obj).__name__
    if "Composite" in obj_type or "Map" in obj_type or "Repeated" in obj_type:
        # Try map-like extraction
        if hasattr(obj, "items") and callable(getattr(obj, "items")):
            try:
                return {str(k): _clean_serializable(v) for k, v in obj.items()}
            except: pass
        # Try sequence-like extraction
        if hasattr(obj, "__iter__"):
            try:
                return [_clean_serializable(item) for item in obj]
            except: pass
            
    # General map fallback
    if hasattr(obj, "items") and callable(getattr(obj, "items")):
        try:
            return {str(k): _clean_serializable(v) for k, v in obj.items()}
        except: pass
            
    # General iteration fallback
    if hasattr(obj, "__iter__"):
        try:
            return [_clean_serializable(item) for item in obj]
        except: pass

    # NumPy/Pandas fallback
    if hasattr(obj, "tolist") and callable(getattr(obj, "tolist")):
        return obj.tolist()

    # Final Fallback to string
    return str(obj)

def _save_mission(mission_id, data):
    path = _get_mission_path(mission_id)
    try:
        clean_data = _clean_serializable(data)
        with mission_lock:
            # Final safety check before dumping
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(clean_data, f, indent=2)
    except Exception as e:
        import traceback
        print(f"[LocalStorage] CRITICAL Error saving mission {mission_id}: {e}")
        # traceback.print_exc()
        # Fallback to a very minimal save to recover state if possible
        try:
            with mission_lock:
                with open(path, 'w', encoding='utf-8') as f:
                    # Only save a few safe fields
                    minimal = {
                        "id": data.get("id"),
                        "status": "error_recovery",
                        "error": str(e)
                    }
                    json.dump(minimal, f)
        except:
             pass

def get_all_projects():
    projects = []
    if not os.path.exists(PROJECTS_DIR):
        os.makedirs(PROJECTS_DIR)
        
    for dirname in os.listdir(PROJECTS_DIR):
        dirpath = os.path.join(PROJECTS_DIR, dirname)
        if os.path.isdir(dirpath):
            json_path = os.path.join(dirpath, "mission_data.json")
            if os.path.exists(json_path):
                try:
                    with open(json_path, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                        # Minimal data for the list
                        projects.append({
                            "id": dirname,
                            "name": data.get("target_name", "Untitled"),
                            "owner": data.get("owner", "Guest"),
                            "date": data.get("date", "Unknown")
                        })
                except:
                    pass
    return projects

def delete_project(mission_id):
    path = os.path.join(PROJECTS_DIR, mission_id)
    if os.path.exists(path):
        shutil.rmtree(path)
        return True
    return False

def link_pdb_to_mission(mission_id, pdb_id, filename):
    """Updates the mission metadata to include a linked PDB structure."""
    mission = get_mission(mission_id)
    if mission:
        if "linked_pdbs" not in mission:
            mission["linked_pdbs"] = []
        
        # Avoid duplicates
        if not any(p.get("id") == pdb_id for p in mission["linked_pdbs"]):
            mission["linked_pdbs"].append({
                "id": pdb_id,
                "filename": filename,
                "timestamp": int(time.time())
            })
            _save_mission(mission_id, mission)
            return True
    return False
