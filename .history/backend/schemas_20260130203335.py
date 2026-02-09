from pydantic import BaseModel

class DiseaseQuery(BaseModel):
    disease_name: str
    max_papers: int = 10

class TargetRequest(BaseModel):
    target_name: str

class RetroRequest(BaseModel):
    smiles: str