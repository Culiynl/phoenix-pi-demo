import firebase_admin
from firebase_admin import credentials, firestore
import os

# Path to the service account key
SERVICE_ACCOUNT_KEY = os.path.join(os.path.dirname(__file__), "..", "serviceAccountKey.json")

_db = None

def init_firebase():
    global _db
    if not firebase_admin._apps:
        if os.path.exists(SERVICE_ACCOUNT_KEY):
            cred = credentials.Certificate(SERVICE_ACCOUNT_KEY)
            firebase_admin.initialize_app(cred)
            print("Initialized Firebase with Service Account Key")
        else:
            # Fallback for Cloud Run (Default Service Account)
            print("Service Account Key not found. Using Application Default Credentials...")
            firebase_admin.initialize_app()
            
    _db = firestore.client()
    return _db

def get_db():
    global _db
    if _db is None:
        init_firebase()
    return _db

def update_mission_log(mission_id, log_entry):
    db = get_db()
    mission_ref = db.collection("missions").document(mission_id)
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
        "research_summary": ""
    })

def update_mission_field(mission_id, field, value):
    db = get_db()
    mission_ref = db.collection("missions").document(mission_id)
    mission_ref.update({
        field: value
    })
