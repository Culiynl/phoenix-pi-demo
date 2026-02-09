import firebase_admin
from firebase_admin import credentials, firestore
import os

# Global variable to hold the Firestore client
_db = None

def init_firebase():
    global _db
    
    # Check if Firebase is already initialized to prevent "App already exists" errors
    if not firebase_admin._apps:
        # 1. Priority: Check if GOOGLE_APPLICATION_CREDENTIALS is set (Render/Production)
        # firebase_admin.initialize_app() automatically finds this variable.
        if os.getenv("GOOGLE_APPLICATION_CREDENTIALS"):
            print(f"Initializing Firebase with Env Var: {os.getenv('GOOGLE_APPLICATION_CREDENTIALS')}")
            firebase_admin.initialize_app()
            
        # 2. Fallback: Check for local serviceAccountKey.json (Local Development)
        else:
            local_key_path = os.path.join(os.path.dirname(__file__), "..", "serviceAccountKey.json")
            if os.path.exists(local_key_path):
                print("Environment variable not set. Falling back to local serviceAccountKey.json...")
                cred = credentials.Certificate(local_key_path)
                firebase_admin.initialize_app(cred)
            else:
                # 3. Final Fallback: GCloud Default Credentials (if running on Google Cloud directly)
                print("No key found. Attempting Default Credentials...")
                firebase_admin.initialize_app()

    if _db is None:
        _db = firestore.client()
        
    return _db

def get_db():
    """Returns the Firestore client, initializing it if necessary."""
    global _db
    if _db is None:
        init_firebase()
    return _db

def update_mission_log(mission_id, log_entry):
    db = get_db()
    mission_ref = db.collection("missions").document(mission_id)
    # Using array_union ensures we don't overwrite existing logs
    mission_ref.update({
        "thought_log": firestore.ArrayUnion([log_entry])
    })

def create_mission(mission_id, target_name, research_mode="hybrid"):
    db = get_db()
    mission_ref = db.collection("missions").document(mission_id)
    mission_ref.set({
        "target_name": target_name,
        "research_mode": research_mode,
        "status": "started",
        "created_at": firestore.SERVER_TIMESTAMP,
        "thought_log": [],
        "tool_history": [],
        "current_step": 1,
        "is_thinking": False,
        "awaiting_confirmation": False,
        "research_summary": "",
        "internal_research_summary": "", # Added for staging
        "internal_research_papers": []   # Added for staging
    })

def update_mission_field(mission_id, field, value):
    db = get_db()
    mission_ref = db.collection("missions").document(mission_id)
    mission_ref.update({
        field: value
    })