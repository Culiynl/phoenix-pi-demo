from pydantic import BaseModel
from typing import Optional, Dict

class DockingRequest(BaseModel):
    project_id: str
    smiles: str
    description: Optional[str] = "Manual docking run"
    
class DiseaseQuery(BaseModel):
    disease_name: str
    max_papers: int = 10

class TargetRequest(BaseModel):
    target_name: str

class RetroRequest(BaseModel):
    smiles: str


\class FeudalRequest(BaseModel):
    project_id: str
    smiles: str
    steps: Optional[int] = 3
    batch_size: Optional[int] = 8

class ParetoRequest(BaseModel):
    project_id: str
    # Weights example: {"affinity": 1.0, "qed": 1.0, "sa": 1.0}
    weights: Dict[str, float]