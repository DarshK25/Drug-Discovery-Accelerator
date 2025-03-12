from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
import random
import json

app = FastAPI(
    title="Drug Discovery Assistant API",
    description="API for a Drug Discovery Assistant powered by Generative AI",
    version="1.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, replace with specific origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class GenerationRequest(BaseModel):
    target_properties: Dict[str, Any]
    constraints: Optional[Dict[str, Any]] = None
    num_molecules: int = 5

class MoleculeData(BaseModel):
    smiles: str
    properties: Dict[str, Any]
    image_url: Optional[str] = None

class GenerationResponse(BaseModel):
    molecules: List[MoleculeData]
    message: str

class PropertyPredictionRequest(BaseModel):
    smiles: str

class SimilaritySearchRequest(BaseModel):
    smiles: str
    num_results: int = 10

@app.get("/")
async def root():
    return {"message": "Welcome to the Drug Discovery Assistant API"}

@app.post("/generate", response_model=GenerationResponse)
async def generate_molecules(request: GenerationRequest):
    try:
        # Sample molecules for demonstration
        known_drugs = [
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen
            "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=CC=C3",  # Anthraquinone
        ]
        
        # Generate molecules (simplified for demonstration)
        molecules = []
        for i in range(request.num_molecules):
            smiles = random.choice(known_drugs)
            
            # Generate properties
            properties = {
                "molecular_weight": round(random.uniform(200, 500), 2),
                "logP": round(random.uniform(0, 5), 2),
                "num_h_donors": random.randint(0, 5),
                "num_h_acceptors": random.randint(0, 10),
                "rotatable_bonds": random.randint(0, 10),
                "qed": round(random.uniform(0, 1), 3),
                "tpsa": round(random.uniform(50, 150), 2),
                "druglikeness": round(random.uniform(0, 1), 2),
                "solubility": random.choice(["Low", "Moderate", "High"]),
                "bioavailability": random.choice(["Low", "Moderate", "High"]),
                "toxicity_risk": random.choice(["Low", "Moderate", "High"]),
                "synthetic_accessibility": round(random.uniform(1, 10), 1),
                "binding_affinity": round(random.uniform(4, 10), 2),
                "is_valid": True
            }
            
            # Generate image URL (using PubChem for demonstration)
            image_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/PNG"
            
            molecules.append(MoleculeData(
                smiles=smiles,
                properties=properties,
                image_url=image_url
            ))
        
        return GenerationResponse(
            molecules=molecules,
            message=f"Successfully generated {len(molecules)} molecules"
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/predict-properties")
async def predict_properties(request: PropertyPredictionRequest):
    try:
        # Generate properties (simplified for demonstration)
        properties = {
            "molecular_weight": round(random.uniform(200, 500), 2),
            "logP": round(random.uniform(0, 5), 2),
            "num_h_donors": random.randint(0, 5),
            "num_h_acceptors": random.randint(0, 10),
            "rotatable_bonds": random.randint(0, 10),
            "qed": round(random.uniform(0, 1), 3),
            "tpsa": round(random.uniform(50, 150), 2),
            "druglikeness": round(random.uniform(0, 1), 2),
            "solubility": random.choice(["Low", "Moderate", "High"]),
            "bioavailability": random.choice(["Low", "Moderate", "High"]),
            "toxicity_risk": random.choice(["Low", "Moderate", "High"]),
            "synthetic_accessibility": round(random.uniform(1, 10), 1),
            "binding_affinity": round(random.uniform(4, 10), 2),
            "is_valid": True
        }
        
        return {
            "smiles": request.smiles,
            "properties": properties
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/similarity-search")
async def search_similar_molecules(request: SimilaritySearchRequest):
    try:
        # Sample database for demonstration
        database = [
            {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "Aspirin", "target": "COX inhibitor"},
            {"smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "name": "Ibuprofen", "target": "COX inhibitor"},
            {"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "name": "Caffeine", "target": "Adenosine receptor antagonist"},
            {"smiles": "CC(=O)NC1=CC=C(C=C1)O", "name": "Acetaminophen", "target": "COX inhibitor"},
            {"smiles": "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=CC=C3", "name": "Anthraquinone", "target": "DNA intercalator"},
            {"smiles": "COc1ccc2cc(ccc2c1)C(=O)C", "name": "Nabumetone", "target": "COX inhibitor"},
            {"smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "name": "Naproxen", "target": "COX inhibitor"},
            {"smiles": "CCN(CC)CCCC(C)Nc1cccc2ccccc12", "name": "Chloroquine", "target": "Heme polymerase inhibitor"},
            {"smiles": "CCOC(=O)C1=C(C)NC(=C(C1C(=O)OC)C(=O)OC)C", "name": "Amlodipine", "target": "Calcium channel blocker"},
            {"smiles": "CC(CS)C(=O)N1CCCC1C(=O)O", "name": "Captopril", "target": "ACE inhibitor"},
        ]
        
        # Generate random similarity scores
        similar_molecules = []
        for entry in database:
            similarity = round(random.uniform(0, 1), 3)
            entry_with_similarity = entry.copy()
            entry_with_similarity["similarity"] = similarity
            similar_molecules.append(entry_with_similarity)
        
        # Sort by similarity (descending)
        similar_molecules.sort(key=lambda x: x["similarity"], reverse=True)
        
        return {
            "query_smiles": request.smiles,
            "similar_molecules": similar_molecules[:request.num_results]
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/optimize")
async def optimize_molecule(smiles: str, target_properties: str):
    try:
        # Parse target properties
        target_props = json.loads(target_properties)
        
        # Sample molecules for demonstration
        known_drugs = [
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen
            "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=CC=C3",  # Anthraquinone
        ]
        
        # Generate optimized molecules (simplified for demonstration)
        optimized_molecules = []
        for i in range(4):
            # Generate a variation of the input molecule
            optimized_smiles = random.choice(known_drugs)
            
            # Generate properties
            properties = {
                "molecular_weight": round(random.uniform(200, 500), 2),
                "logP": round(random.uniform(0, 5), 2),
                "num_h_donors": random.randint(0, 5),
                "num_h_acceptors": random.randint(0, 10),
                "rotatable_bonds": random.randint(0, 10),
                "qed": round(random.uniform(0, 1), 3),
                "tpsa": round(random.uniform(50, 150), 2),
                "druglikeness": round(random.uniform(0, 1), 2),
                "solubility": random.choice(["Low", "Moderate", "High"]),
                "bioavailability": random.choice(["Low", "Moderate", "High"]),
                "toxicity_risk": random.choice(["Low", "Moderate", "High"]),
                "synthetic_accessibility": round(random.uniform(1, 10), 1),
                "binding_affinity": round(random.uniform(4, 10), 2),
                "is_valid": True
            }
            
            # Generate image URL (using PubChem for demonstration)
            image_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{optimized_smiles}/PNG"
            
            optimized_molecules.append({
                "smiles": optimized_smiles,
                "properties": properties,
                "image_url": image_url
            })
        
        return {
            "original_smiles": smiles,
            "optimized_molecules": optimized_molecules
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("simplified_main:app", host="0.0.0.0", port=8000, reload=True) 