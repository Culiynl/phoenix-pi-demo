import firebase_admin
from firebase_admin import credentials, firestore
import os

# Global variable to hold the Firestore client
_db = None

def init_firebase():
    global _db
    if not firebase_admin._apps:
        # 1. Path for Render (Production)
        render_secret_path = "/etc/secrets/service_account_key.json"
        
        # 2. Path for your Local Machine (Fixed)
        # We go up ONE level (..) from the 'services' folder to hit the 'backend' folder
        local_key_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "serviceAccountKey.json"))

        if os.path.exists(render_secret_path):
            print(f"Using Render Secret Key: {render_secret_path}")
            cred = credentials.Certificate(render_secret_path)
            firebase_admin.initialize_app(cred)
        elif os.path.exists(local_key_path):
            print(f"Using Local Key: {local_key_path}")
            cred = credentials.Certificate(local_key_path)
            firebase_admin.initialize_app(cred)
        else:
            # This helps you debug exactly where the code is looking
            print(f"DEBUG: Looked for local key at: {local_key_path}")
            raise Exception("CRITICAL: No serviceAccountKey.json found! Firestore will not work.")
            
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