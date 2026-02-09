from pydantic import BaseModel
from pydantic import BaseModel
from typing import Optional

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