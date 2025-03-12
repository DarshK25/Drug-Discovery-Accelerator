import axios from 'axios';

// Create an axios instance with default config
const api = axios.create({
  baseURL: 'http://localhost:8000',
  headers: {
    'Content-Type': 'application/json',
  },
});

// API functions for drug discovery
const drugDiscoveryApi = {
  // Generate molecules based on target properties
  generateMolecules: async (targetProperties, constraints, numMolecules = 5) => {
    try {
      const response = await api.post('/generate', {
        target_properties: targetProperties,
        constraints: constraints,
        num_molecules: numMolecules,
      });
      return response.data;
    } catch (error) {
      console.error('Error generating molecules:', error);
      // For demo purposes, return mock data if the API is not available
      return mockGenerateMolecules(targetProperties, constraints, numMolecules);
    }
  },

  // Predict properties for a molecule
  predictProperties: async (smiles) => {
    try {
      const response = await api.post('/predict-properties', {
        smiles: smiles,
      });
      return response.data;
    } catch (error) {
      console.error('Error predicting properties:', error);
      // For demo purposes, return mock data if the API is not available
      return mockPredictProperties(smiles);
    }
  },

  // Search for similar molecules
  searchSimilarMolecules: async (smiles, numResults = 10) => {
    try {
      const response = await api.post('/similarity-search', {
        smiles: smiles,
        num_results: numResults,
      });
      return response.data;
    } catch (error) {
      console.error('Error searching similar molecules:', error);
      // For demo purposes, return mock data if the API is not available
      return mockSearchSimilarMolecules(smiles, numResults);
    }
  },

  // Optimize a molecule
  optimizeMolecule: async (smiles, targetProperties) => {
    try {
      const formData = new FormData();
      formData.append('smiles', smiles);
      formData.append('target_properties', JSON.stringify(targetProperties));
      
      const response = await api.post('/optimize', formData, {
        headers: {
          'Content-Type': 'multipart/form-data',
        },
      });
      return response.data;
    } catch (error) {
      console.error('Error optimizing molecule:', error);
      // For demo purposes, return mock data if the API is not available
      return mockOptimizeMolecule(smiles, targetProperties);
    }
  },
};

// Mock functions for demo purposes
const mockGenerateMolecules = (targetProperties, constraints, numMolecules) => {
  const knownDrugs = [
    "CC(=O)Oc1ccccc1C(=O)O",  // Aspirin
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  // Ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  // Caffeine
    "CC(=O)NC1=CC=C(C=C1)O",  // Acetaminophen
    "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=CC=C3",  // Anthraquinone
  ];
  
  const molecules = [];
  for (let i = 0; i < numMolecules; i++) {
    const smiles = knownDrugs[Math.floor(Math.random() * knownDrugs.length)];
    
    molecules.push({
      smiles,
      properties: generateRandomProperties(),
      image_url: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/PNG`
    });
  }
  
  return {
    molecules,
    message: `Successfully generated ${molecules.length} molecules`
  };
};

const mockPredictProperties = (smiles) => {
  return {
    smiles,
    properties: generateRandomProperties()
  };
};

const mockSearchSimilarMolecules = (smiles, numResults) => {
  const database = [
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
  ];
  
  const similarMolecules = database.map(entry => ({
    ...entry,
    similarity: Math.random()
  }));
  
  similarMolecules.sort((a, b) => b.similarity - a.similarity);
  
  return {
    query_smiles: smiles,
    similar_molecules: similarMolecules.slice(0, numResults)
  };
};

const mockOptimizeMolecule = (smiles, targetProperties) => {
  const knownDrugs = [
    "CC(=O)Oc1ccccc1C(=O)O",  // Aspirin
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  // Ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  // Caffeine
    "CC(=O)NC1=CC=C(C=C1)O",  // Acetaminophen
    "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=CC=C3",  // Anthraquinone
  ];
  
  const optimizedMolecules = [];
  for (let i = 0; i < 4; i++) {
    const optimizedSmiles = knownDrugs[Math.floor(Math.random() * knownDrugs.length)];
    
    optimizedMolecules.push({
      smiles: optimizedSmiles,
      properties: generateRandomProperties(),
      image_url: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(optimizedSmiles)}/PNG`
    });
  }
  
  return {
    original_smiles: smiles,
    optimized_molecules: optimizedMolecules
  };
};

// Helper function to generate random properties
const generateRandomProperties = () => {
  return {
    molecular_weight: Math.round((Math.random() * 300 + 200) * 100) / 100,
    logP: Math.round((Math.random() * 5) * 100) / 100,
    num_h_donors: Math.floor(Math.random() * 6),
    num_h_acceptors: Math.floor(Math.random() * 11),
    rotatable_bonds: Math.floor(Math.random() * 11),
    qed: Math.round((Math.random()) * 1000) / 1000,
    tpsa: Math.round((Math.random() * 100 + 50) * 100) / 100,
    druglikeness: Math.round((Math.random()) * 100) / 100,
    solubility: ["Low", "Moderate", "High"][Math.floor(Math.random() * 3)],
    bioavailability: ["Low", "Moderate", "High"][Math.floor(Math.random() * 3)],
    toxicity_risk: ["Low", "Moderate", "High"][Math.floor(Math.random() * 3)],
    synthetic_accessibility: Math.round((Math.random() * 9 + 1) * 10) / 10,
    binding_affinity: Math.round((Math.random() * 6 + 4) * 100) / 100,
    is_valid: true
  };
};

export default drugDiscoveryApi; 