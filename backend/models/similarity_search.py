import numpy as np
from typing import List, Dict, Any
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

class SimilaritySearch:
    """
    A class that searches for molecules similar to a query molecule.
    This is a simplified implementation that uses RDKit fingerprints.
    In a real implementation, this would use a database of molecules.
    """
    
    def __init__(self):
        """Initialize the similarity search model."""
        # In a real implementation, this would load a database of molecules
        self.model_loaded = True
        
        # Create a small database of known drugs for demonstration
        self.database = [
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
            {"smiles": "CC(C)NCC(O)COc1cccc2ccccc12", "name": "Propranolol", "target": "Beta blocker"},
            {"smiles": "CCCC1=NN(C(=O)N1)c1ccc(cc1)N1CCN(CC1)C", "name": "Sildenafil", "target": "PDE5 inhibitor"},
            {"smiles": "Clc1ccccc1C(=O)NCCCn1ccnc1", "name": "Clonazepam", "target": "GABA receptor modulator"},
            {"smiles": "COc1ccc2[nH]cc(CCN(C)C)c2c1", "name": "Dimethyltryptamine", "target": "Serotonin receptor agonist"},
            {"smiles": "CN1C2CCC1CC(C2)OC(=O)C(CO)c1ccccc1", "name": "Cocaine", "target": "Dopamine transporter inhibitor"},
            {"smiles": "CNC(=O)Oc1ccccc1OC(=O)N(C)C", "name": "Physostigmine", "target": "Acetylcholinesterase inhibitor"},
            {"smiles": "CN1c2ccc(cc2C(=NCC1=O)c1ccccc1)Cl", "name": "Diazepam", "target": "GABA receptor modulator"},
            {"smiles": "CC(=O)N(c1ccc(OC(F)(F)F)cc1)C1(CC1)C#N", "name": "Efavirenz", "target": "HIV reverse transcriptase inhibitor"},
            {"smiles": "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", "name": "Salbutamol", "target": "Beta-2 adrenergic receptor agonist"},
            {"smiles": "CC(C)C(=O)OC1(CCN(C)CC1)c1ccccc1", "name": "Atropine", "target": "Muscarinic acetylcholine receptor antagonist"}
        ]
        
        # Compute fingerprints for all molecules in the database
        self.fingerprints = []
        for entry in self.database:
            mol = Chem.MolFromSmiles(entry["smiles"])
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                self.fingerprints.append(fp)
            else:
                self.fingerprints.append(None)
        
        print("Similarity Search initialized")
    
    def search(self, query_smiles: str, num_results: int = 10) -> List[Dict[str, Any]]:
        """
        Search for molecules similar to the query molecule.
        
        Args:
            query_smiles: SMILES string of the query molecule
            num_results: Number of results to return
            
        Returns:
            List of dictionaries containing similar molecules and their similarity scores
        """
        # Convert query SMILES to RDKit molecule
        query_mol = Chem.MolFromSmiles(query_smiles)
        
        if not query_mol:
            # Return empty list if the query molecule is invalid
            return []
        
        # Compute fingerprint for the query molecule
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=1024)
        
        # Compute similarity scores
        similarities = []
        for i, fp in enumerate(self.fingerprints):
            if fp:
                similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
                similarities.append((similarity, i))
        
        # Sort by similarity (descending)
        similarities.sort(reverse=True)
        
        # Return the top results
        results = []
        for similarity, idx in similarities[:num_results]:
            entry = self.database[idx].copy()
            entry["similarity"] = round(similarity, 3)
            results.append(entry)
        
        return results 