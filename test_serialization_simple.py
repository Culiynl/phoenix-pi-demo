
import sys
import os
import json

# Add backend to path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from services.local_storage_service import _clean_serializable

class MockMapComposite:
    def __init__(self, data):
        self.data = data
    def items(self):
        return self.data.items()

def test_serialization():
    print("Testing serialization logic...")
    
    complex_obj = {
        "id": "test_mission",
        "nested": MockMapComposite({"key": "value", "list": [1, 2, MockMapComposite({"inner": "inner_val"})]}),
        "proto_like": MockMapComposite({"a": 1})
    }
    
    cleaned = _clean_serializable(complex_obj)
    print(f"Cleaned object: {json.dumps(cleaned, indent=2)}")
    
    assert cleaned["nested"]["key"] == "value"
    assert cleaned["nested"]["list"][2]["inner"] == "inner_val"
    print("Test PASSED!")

if __name__ == "__main__":
    test_serialization()
