from pydantic import BaseModel
from typing import Optional, List

class PredictionRequest(BaseModel):
    method: str = "simple"
    num_conformations: int = 5

class PredictionResult(BaseModel):
    job_id: str
    results: Optional[dict]
    status: Optional[str] = "pending" 