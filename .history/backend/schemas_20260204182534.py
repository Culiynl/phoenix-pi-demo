from pydantic import BaseModel
from typing import Optional, Dict

class DockingRequest(BaseModel):
    project_id: str
    smiles: str
    center_override: Optional[Dict[str, float]] = None
    
class DiseaseQuery(BaseModel):
    disease_name: str
    max_papers: int = 10

class TargetRequest(BaseModel):
    target_name: str

class RetroRequest(BaseModel):
    smiles: str

class FeudalRequest(BaseModel):
    project_id: str
    smiles: str
    steps: Optional[int] = 3
    batch_size: Optional[int] = 8

class ParetoRequest(BaseModel):
    project_id: str
    weights: Dict[str, float]