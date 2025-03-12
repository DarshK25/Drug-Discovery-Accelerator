import base64
import io
from typing import Optional
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

def validate_smiles(smiles: str) -> bool:
    """
    Validate a SMILES string.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        True if the SMILES string is valid, False otherwise
    """
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def smiles_to_image(smiles: str, width: int = 300, height: int = 200) -> Optional[str]:
    """
    Convert a SMILES string to a base64-encoded PNG image.
    
    Args:
        smiles: SMILES string to convert
        width: Width of the image in pixels
        height: Height of the image in pixels
        
    Returns:
        Base64-encoded PNG image, or None if the SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    # Generate 2D coordinates for the molecule
    AllChem.Compute2DCoords(mol)
    
    # Draw the molecule
    drawer = Draw.MolDraw2DCairo(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    # Get the PNG data
    png_data = drawer.GetDrawingText()
    
    # Convert to base64
    encoded = base64.b64encode(png_data).decode('utf-8')
    
    # Return as a data URL
    return f"data:image/png;base64,{encoded}"

def canonical_smiles(smiles: str) -> Optional[str]:
    """
    Convert a SMILES string to its canonical form.
    
    Args:
        smiles: SMILES string to convert
        
    Returns:
        Canonical SMILES string, or None if the input is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def calculate_molecular_formula(smiles: str) -> Optional[str]:
    """
    Calculate the molecular formula for a SMILES string.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Molecular formula, or None if the SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    # Get the atomic symbols and counts
    atom_dict = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in atom_dict:
            atom_dict[symbol] += 1
        else:
            atom_dict[symbol] = 1
    
    # Order the elements (C and H first, then alphabetically)
    ordered_elements = []
    if 'C' in atom_dict:
        ordered_elements.append(('C', atom_dict['C']))
        del atom_dict['C']
    if 'H' in atom_dict:
        ordered_elements.append(('H', atom_dict['H']))
        del atom_dict['H']
    
    # Add the remaining elements in alphabetical order
    for symbol in sorted(atom_dict.keys()):
        ordered_elements.append((symbol, atom_dict[symbol]))
    
    # Build the formula string
    formula = ''
    for symbol, count in ordered_elements:
        if count == 1:
            formula += symbol
        else:
            formula += f"{symbol}{count}"
    
    return formula

def get_scaffold(smiles: str) -> Optional[str]:
    """
    Get the Murcko scaffold of a molecule.
    
    Args:
        smiles: SMILES string
        
    Returns:
        SMILES string of the scaffold, or None if the input is invalid
    """
    from rdkit.Chem.Scaffolds import MurckoScaffold
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold) if scaffold else None