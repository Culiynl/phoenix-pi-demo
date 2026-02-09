import time

# Dictionary to hold mutable state
training_state = {
    "is_running": False,
    "stage": "Idle",
    "progress": 0,
    "status_text": "Idle",
    "logs": [],
    "best_r2": 0.0,
    "current_pocket_center": None,
    "viz_image": None 
}

def log_msg(msg: str):
    """Logs a message to console and the global state."""
    print(msg, flush=True)
    timestamp = time.strftime("%H:%M:%S")
    training_state["logs"].append(f"[{timestamp}] {msg}")

def reset_state():
    training_state["is_running"] = True
    training_state["logs"] = []
    training_state["progress"] = 0
    training_state["best_r2"] = 0.0
    training_state["status_text"] = "Initializing..."