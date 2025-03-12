# Drug Discovery Assistant Backend

This is the backend for the Drug Discovery Assistant, a tool that uses generative AI to accelerate drug discovery.

## Features

- Generate novel drug candidates based on target properties
- Predict molecular properties for drug candidates
- Search for similar molecules to a query molecule
- Optimize existing molecules to improve properties

## Setup

1. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the server:
```bash
uvicorn main:app --reload
```

The API will be available at http://localhost:8000.

## API Endpoints

- `GET /`: Welcome message
- `POST /generate`: Generate novel drug candidates
- `POST /predict-properties`: Predict properties for a molecule
- `POST /similarity-search`: Search for similar molecules
- `POST /optimize`: Optimize a molecule to improve properties

## Models

The backend includes the following models:

- `MoleculeGenerator`: Generates novel drug candidates
- `PropertyPredictor`: Predicts molecular properties
- `SimilaritySearch`: Searches for similar molecules

## Technologies Used

- FastAPI: Web framework
- RDKit: Cheminformatics toolkit
- PyTorch: Deep learning framework
- Transformers: NLP models for molecular generation

## Future Improvements

- Integration with real generative AI models (VAEs, GANs, transformers)
- Database of molecules for similarity search
- More sophisticated property prediction models
- Integration with external APIs for additional data 