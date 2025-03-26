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

- `MoleculeGenerator`: Generates novel drug candidates using generative AI models
- `PropertyPredictor`: Predicts molecular properties
- `SimilaritySearch`: Searches for similar molecules

### Generative Models

The `MoleculeGenerator` class now implements three types of generative models:

1. **Transformer-based Model**: Uses ChemBERTa (a BERT model fine-tuned on chemical data) to generate novel molecules. This is the default model.

2. **Variational Autoencoder (VAE)**: Encodes molecules into a latent space and generates new molecules by sampling from this space.

3. **Generative Adversarial Network (GAN)**: Currently a placeholder for future implementation.

#### Selecting a Model Type

You can select which generative model to use by modifying the `model_type` parameter in the `MoleculeGenerator` class:

```python
# In backend/models/molecule_generator.py
def __init__(self):
    self.model_type = "transformer"  # Options: "vae", "gan", "transformer"
```

#### Training Your Own Models

To train your own models:

1. **For VAE**:
   - Collect a dataset of SMILES strings
   - Convert them to molecular fingerprints
   - Train the VAE model on these fingerprints
   - Save the model weights to a file
   - Update the VAE model loading code to load your weights

2. **For Transformer**:
   - Fine-tune a pre-trained model like ChemBERTa on your SMILES dataset
   - Save the model and tokenizer
   - Update the transformer model loading code to use your model

## Technologies Used

- FastAPI: Web framework
- RDKit: Cheminformatics toolkit
- PyTorch: Deep learning framework
- Transformers: NLP models for molecular generation
- Accelerate: For faster model inference

## Running with Real Models

To run with real generative models:

1. Make sure you have a GPU available (recommended) or enough CPU resources
2. Install all dependencies including PyTorch, Transformers, and Accelerate
3. If using your own pre-trained models, place them in a `models` directory
4. Update the model paths in the code if necessary
5. Run the server as usual

## Future Improvements

- Full implementation of GAN-based molecule generation
- Fine-tuning models on specific chemical domains
- Conditional generation based on target properties
- Database of molecules for similarity search
- More sophisticated property prediction models
- Integration with external APIs for additional data 