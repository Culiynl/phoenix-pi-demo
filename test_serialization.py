import sys
import os
import json

# Mock projects dir
os.environ["PROJECTS_DIR"] = "./test_projects"
if not os.path.exists("./test_projects"):
    os.makedirs("./test_projects")

# Add backend to path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from services.local_storage_service import _clean_serializable, _save_mission

class MockMapComposite:
    def __init__(self, data):
        self.data = data
    def items(self):
        return self.data.items()

def test_serialization():
    print("Testing serialization...")
    
    complex_obj = {
        "id": "test_mission",
        "nested": MockMapComposite({"key": "value", "list": [1, 2, MockMapComposite({"inner": "inner_val"})]}),
        "proto_like": MockMapComposite({"a": 1})
    }
    
    cleaned = _clean_serializable(complex_obj)
    print(f"Cleaned object: {cleaned}")
    
    # Try saving
    _save_mission("test_mission", complex_obj)
    
    # Verify file
    path = os.path.join("./test_projects", "test_mission", "mission_data.json")
    if os.path.exists(path):
        with open(path, 'r') as f:
            data = json.load(f)
            print("Successfully saved and loaded JSON!")
            print(json.dumps(data, indent=2))
            assert data["nested"]["key"] == "value"
            assert data["nested"]["list"][2]["inner"] == "inner_val"
    else:
        print("File was not saved!")
        sys.exit(1)

if __name__ == "__main__":
    test_serialization()
