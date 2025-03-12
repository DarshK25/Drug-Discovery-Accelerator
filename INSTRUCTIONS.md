# Drug Discovery Assistant - Running Instructions

This document provides instructions for running the Drug Discovery Assistant application.

## Prerequisites

- Python 3.8+ installed
- Node.js 14+ installed
- npm installed

## Running the Frontend Only (Recommended)

The frontend includes mock data functionality that allows it to run without the backend. This is the recommended approach for demonstration purposes.

1. Navigate to the frontend directory:
```bash
cd frontend
```

2. Install dependencies:
```bash
npm install
```

3. Start the frontend:
```bash
npm start
```

The application will be available at http://localhost:3000.

## Running Both Frontend and Backend (Optional)

If you want to run both the frontend and backend:

1. Install backend dependencies:
```bash
cd backend
pip install fastapi uvicorn pydantic
```

2. Run the simplified backend:
```bash
cd backend
python simplified_main.py
```

3. In a separate terminal, run the frontend:
```bash
cd frontend
npm install
npm start
```

## Using the Application

1. **Home Page**: Overview of the application and its features
2. **Generate Page**: Generate novel drug candidates based on target properties
3. **Optimize Page**: Optimize existing molecules to improve properties
4. **Search Page**: Search for similar molecules to a query molecule
5. **About Page**: Information about the application and its technologies

## Example Molecules

You can use these example SMILES strings to test the application:

- Aspirin: `CC(=O)Oc1ccccc1C(=O)O`
- Ibuprofen: `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O`
- Caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
- Acetaminophen: `CC(=O)NC1=CC=C(C=C1)O`

## Troubleshooting

- If you encounter issues with the backend, you can still use the frontend with mock data.
- If npm start fails, try clearing npm cache with `npm cache clean --force` and then reinstall dependencies.
- Make sure you have the correct versions of Node.js and Python installed. 