from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import uvicorn
import json
from models.molecule_generator import MoleculeGenerator
from models.property_predictor import PropertyPredictor
from models.similarity_search import SimilaritySearch
from utils.molecule_utils import smiles_to_image, validate_smiles

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

# Initialize models
molecule_generator = MoleculeGenerator()
property_predictor = PropertyPredictor()
similarity_search = SimilaritySearch()

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
        # Generate molecules based on target properties and constraints
        molecules = molecule_generator.generate(
            target_properties=request.target_properties,
            constraints=request.constraints,
            num_molecules=request.num_molecules
        )
        
        # Convert to response format
        molecule_data = []
        for mol in molecules:
            # Predict properties for each molecule
            properties = property_predictor.predict(mol)
            
            # Generate image URL
            image_url = smiles_to_image(mol)
            
            molecule_data.append(MoleculeData(
                smiles=mol,
                properties=properties,
                image_url=image_url
            ))
        
        return GenerationResponse(
            molecules=molecule_data,
            message=f"Successfully generated {len(molecule_data)} molecules"
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/predict-properties")
async def predict_properties(request: PropertyPredictionRequest):
    try:
        # Validate SMILES
        if not validate_smiles(request.smiles):
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # Predict properties
        properties = property_predictor.predict(request.smiles)
        
        return {
            "smiles": request.smiles,
            "properties": properties
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/similarity-search")
async def search_similar_molecules(request: SimilaritySearchRequest):
    try:
        # Validate SMILES
        if not validate_smiles(request.smiles):
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # Search for similar molecules
        similar_molecules = similarity_search.search(
            query_smiles=request.smiles,
            num_results=request.num_results
        )
        
        return {
            "query_smiles": request.smiles,
            "similar_molecules": similar_molecules
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/optimize")
async def optimize_molecule(smiles: str = Form(...), target_properties: str = Form(...)):
    try:
        # Validate SMILES
        if not validate_smiles(smiles):
            raise HTTPException(status_code=400, detail="Invalid SMILES string")
        
        # Parse target properties
        target_props = json.loads(target_properties)
        
        # Optimize molecule
        optimized_molecules = molecule_generator.optimize(
            starting_smiles=smiles,
            target_properties=target_props
        )
        
        # Convert to response format
        molecule_data = []
        for mol in optimized_molecules:
            # Predict properties for each molecule
            properties = property_predictor.predict(mol)
            
            # Generate image URL
            image_url = smiles_to_image(mol)
            
            molecule_data.append({
                "smiles": mol,
                "properties": properties,
                "image_url": image_url
            })
        
        return {
            "original_smiles": smiles,
            "optimized_molecules": molecule_data
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True) 