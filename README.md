# Drug Discovery Assistant

A comprehensive platform that leverages generative AI to accelerate drug discovery by enabling researchers to efficiently identify, design, and optimize potential drug candidates.

## Overview

The Drug Discovery Assistant is a full-stack application that combines modern web technologies with advanced AI techniques to streamline the drug discovery process. It provides tools for generating novel molecules, predicting molecular properties, searching for similar compounds, and optimizing existing molecules.

## Features

- **Generate Novel Molecules**: Create new drug candidates based on target properties using generative AI
- **Predict Molecular Properties**: Predict key properties of drug candidates to assess their potential
- **Search Similar Compounds**: Find molecules similar to a query structure in our database
- **Optimize Existing Molecules**: Improve existing drug candidates to enhance specific properties

## Architecture

The application consists of two main components:

### Backend (FastAPI)

- Built with FastAPI, a modern Python web framework
- Uses RDKit for cheminformatics operations
- Implements generative models for molecule generation
- Provides property prediction and similarity search capabilities

### Frontend (React)

- Built with React and Material-UI
- Provides an intuitive user interface for all features
- Visualizes molecular structures and properties
- Offers interactive tools for molecule design and optimization

## Getting Started

### Prerequisites

- Python 3.8+
- Node.js 14+
- npm or yarn

### Backend Setup

1. Navigate to the backend directory:
```bash
cd backend
```

2. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Start the server:
```bash
uvicorn main:app --reload
```

The API will be available at http://localhost:8000.

### Frontend Setup

1. Navigate to the frontend directory:
```bash
cd frontend
```

2. Install dependencies:
```bash
npm install
```

3. Start the development server:
```bash
npm start
```

The application will be available at http://localhost:3000.

## Technologies Used

- **Backend**:
  - FastAPI: Web framework
  - RDKit: Cheminformatics toolkit
  - PyTorch: Deep learning framework
  - Transformers: NLP models for molecular generation

- **Frontend**:
  - React: Frontend library
  - Material-UI: UI component library
  - Axios: HTTP client
  - Chart.js: Data visualization

## Future Improvements

- Integration with real generative AI models (VAEs, GANs, transformers)
- Database of molecules for similarity search
- More sophisticated property prediction models
- User authentication and saved molecules
- Integration with external APIs for additional data
- Batch processing of molecules
- Collaborative features for research teams

## License

This project is licensed under the MIT License - see the LICENSE file for details. # Drug-Discovery-Accelerator
