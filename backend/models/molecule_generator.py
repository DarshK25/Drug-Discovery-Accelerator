import torch
import numpy as np
from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import random
from transformers import AutoTokenizer, AutoModelForCausalLM
import os
from torch import nn
import torch.nn.functional as F

# Import the MoleculeGAN class
try:
    from .molecule_gan import MoleculeGAN
except ImportError:
    # For direct module execution
    try:
        from molecule_gan import MoleculeGAN
    except ImportError:
        print("Warning: Could not import MoleculeGAN")
        MoleculeGAN = None

class MoleculeVAE(nn.Module):
    """
    Variational Autoencoder for molecule generation
    """
    def __init__(self, input_dim=512, hidden_dim=256, latent_dim=64):
        super(MoleculeVAE, self).__init__()
        
        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )
        
        # Latent space
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_var = nn.Linear(hidden_dim, latent_dim)
        
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim),
            nn.Sigmoid()
        )
    
    def encode(self, x):
        h = self.encoder(x)
        mu = self.fc_mu(h)
        log_var = self.fc_var(h)
        return mu, log_var
    
    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        z = mu + eps * std
        return z
    
    def decode(self, z):
        return self.decoder(z)
    
    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        return self.decode(z), mu, log_var

class MoleculeGenerator:
    """
    A class that uses generative AI to create novel drug candidates.
    Implements multiple generative models: VAE, GAN, and Transformer-based approaches.
    """
    
    def __init__(self):
        """Initialize the molecule generator with different model types."""
        self.model_type = "transformer"  # Options: "vae", "gan", "transformer"
        self.model_loaded = False
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # Fallback fragments for when models fail
        self.fragments = [
            "c1ccccc1", "C1CCNCC1", "c1ccncc1", "C(=O)N", "C(=O)O",
            "CC(=O)N", "CN", "CF", "CCl", "CBr", "c1cccnc1", "c1cnccn1"
        ]
        
        # Try to load models
        try:
            self._load_models()
            self.model_loaded = True
            print(f"Molecule Generator initialized with {self.model_type} model")
        except Exception as e:
            print(f"Failed to load generative models: {e}")
            print("Falling back to simplified generation")
    
    def _load_models(self):
        """Load the appropriate generative model based on model_type."""
        if self.model_type == "transformer":
            # Load a pre-trained transformer model for molecule generation
            # In a production environment, you would use a model fine-tuned on SMILES
            try:
                # Try to load a pre-trained model for molecule generation
                # For demonstration, we'll use a small GPT-2 model
                # In production, use a model specifically fine-tuned on SMILES strings
                model_name = "seyonec/ChemBERTa-zinc-base-v1"
                self.tokenizer = AutoTokenizer.from_pretrained(model_name)
                self.model = AutoModelForCausalLM.from_pretrained(model_name)
                self.model.to(self.device)
            except Exception as e:
                print(f"Could not load transformer model: {e}")
                print("Using simplified model instead")
                self.model_type = "simplified"
        
        elif self.model_type == "vae":
            # Initialize and load VAE model
            try:
                self.vae_model = MoleculeVAE()
                # In a real implementation, you would load pre-trained weights
                # self.vae_model.load_state_dict(torch.load("path_to_vae_weights.pt"))
                self.vae_model.to(self.device)
                self.vae_model.eval()
            except Exception as e:
                print(f"Could not load VAE model: {e}")
                print("Using simplified model instead")
                self.model_type = "simplified"
        
        elif self.model_type == "gan":
            # Initialize GAN model
            try:
                if MoleculeGAN is not None:
                    self.gan_model = MoleculeGAN(device=self.device)
                    
                    # Check if pre-trained models exist and load them
                    model_dir = os.path.join(os.path.dirname(__file__), "pretrained_models")
                    generator_path = os.path.join(model_dir, "generator.pt")
                    discriminator_path = os.path.join(model_dir, "discriminator.pt")
                    
                    if os.path.exists(generator_path) and os.path.exists(discriminator_path):
                        self.gan_model.load_models(generator_path, discriminator_path)
                    else:
                        print("Pre-trained GAN models not found. Using untrained model.")
                else:
                    print("MoleculeGAN module not available")
                    self.model_type = "simplified"
            except Exception as e:
                print(f"Could not initialize GAN model: {e}")
                print("Using simplified model instead")
                self.model_type = "simplified"
    
    def _smiles_to_fingerprint(self, smiles):
        """Convert SMILES to molecular fingerprint for VAE input."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)
    
    def _fingerprint_to_tensor(self, fp):
        """Convert fingerprint to PyTorch tensor."""
        return torch.tensor(np.array(fp), dtype=torch.float32).to(self.device)
    
    def _generate_with_transformer(self, num_molecules=5, max_length=100):
        """Generate molecules using a transformer model."""
        generated_smiles = []
        
        try:
            # Use a seed/prompt for generation
            prompt = "C"  # Start with a carbon atom
            
            for _ in range(num_molecules):
                inputs = self.tokenizer(prompt, return_tensors="pt").to(self.device)
                
                # Generate new tokens
                outputs = self.model.generate(
                    inputs.input_ids,
                    max_length=max_length,
                    do_sample=True,
                    top_p=0.9,
                    temperature=1.2,
                    num_return_sequences=1
                )
                
                # Decode the generated tokens to text
                generated_text = self.tokenizer.decode(outputs[0], skip_special_tokens=True)
                
                # Clean up the generated SMILES and validate
                smiles = generated_text.strip()
                mol = Chem.MolFromSmiles(smiles)
                
                if mol:
                    canonical_smiles = Chem.MolToSmiles(mol)
                    generated_smiles.append(canonical_smiles)
            
            return generated_smiles
        except Exception as e:
            print(f"Error in transformer generation: {e}")
            return []
    
    def _generate_with_vae(self, num_molecules=5):
        """Generate molecules using the VAE model."""
        generated_smiles = []
        
        try:
            # Sample from latent space
            for _ in range(num_molecules):
                # Sample from normal distribution
                z = torch.randn(1, 64).to(self.device)
                
                # Decode the latent vector
                with torch.no_grad():
                    decoded = self.vae_model.decode(z)
                
                # Convert back to SMILES (this is a simplified placeholder)
                # In a real implementation, you would convert the decoded tensor
                # back to a valid SMILES string using a proper decoder
                
                # For demonstration, we'll use a random known drug
                known_drugs = [
                    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
                    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
                    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
                    "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen
                ]
                generated_smiles.append(random.choice(known_drugs))
            
            return generated_smiles
        except Exception as e:
            print(f"Error in VAE generation: {e}")
            return []
    
    def _generate_with_gan(self, num_molecules=5):
        """Generate molecules using a GAN model."""
        try:
            if hasattr(self, 'gan_model'):
                return self.gan_model.generate_molecules(num_molecules)
            else:
                print("GAN model not initialized")
                return self._generate_simplified(num_molecules)
        except Exception as e:
            print(f"Error in GAN generation: {e}")
            return self._generate_simplified(num_molecules)
    
    def _generate_simplified(self, num_molecules=5):
        """Fallback method for molecule generation when models fail."""
        generated_molecules = []
        
        for _ in range(num_molecules):
            # Generate a random molecule by combining fragments
            num_fragments = random.randint(2, 5)
            selected_fragments = random.sample(self.fragments, num_fragments)
            
            # Create a molecule from the fragments
            mol = None
            for fragment in selected_fragments:
                frag_mol = Chem.MolFromSmiles(fragment)
                if mol is None:
                    mol = frag_mol
                else:
                    combo = Chem.CombineMols(mol, frag_mol)
                    mol = combo
            
            # Convert to SMILES
            if mol:
                try:
                    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
                    if mol:
                        smiles = Chem.MolToSmiles(mol)
                        generated_molecules.append(smiles)
                except:
                    # If conversion fails, use a default molecule
                    generated_molecules.append("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
        
        # If we couldn't generate enough valid molecules, add some known drugs
        known_drugs = [
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen
            "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=CC=C3",  # Anthraquinone
        ]
        
        while len(generated_molecules) < num_molecules:
            generated_molecules.append(random.choice(known_drugs))
        
        return generated_molecules[:num_molecules]
    
    def generate(self, target_properties: Dict[str, Any], 
                 constraints: Optional[Dict[str, Any]] = None, 
                 num_molecules: int = 5) -> List[str]:
        """
        Generate molecules based on target properties and constraints.
        
        Args:
            target_properties: Dictionary of desired molecular properties
            constraints: Optional constraints on the generation process
            num_molecules: Number of molecules to generate
            
        Returns:
            List of SMILES strings representing the generated molecules
        """
        if not self.model_loaded:
            return self._generate_simplified(num_molecules)
        
        generated_molecules = []
        
        # Generate molecules based on the selected model type
        if self.model_type == "transformer":
            generated_molecules = self._generate_with_transformer(num_molecules)
        elif self.model_type == "vae":
            generated_molecules = self._generate_with_vae(num_molecules)
        elif self.model_type == "gan":
            generated_molecules = self._generate_with_gan(num_molecules)
        else:
            generated_molecules = self._generate_simplified(num_molecules)
        
        # If we couldn't generate enough molecules, fall back to simplified method
        if len(generated_molecules) < num_molecules:
            additional_molecules = self._generate_simplified(num_molecules - len(generated_molecules))
            generated_molecules.extend(additional_molecules)
        
        # Apply property-based filtering (in a real implementation, this would be more sophisticated)
        # For now, we'll just return the generated molecules
        
        return generated_molecules[:num_molecules]
    
    def optimize(self, starting_smiles: str, target_properties: Dict[str, Any]) -> List[str]:
        """
        Optimize a molecule to improve certain properties using generative models.
        
        Args:
            starting_smiles: SMILES string of the starting molecule
            target_properties: Dictionary of desired molecular properties
            
        Returns:
            List of SMILES strings representing optimized molecules
        """
        mol = Chem.MolFromSmiles(starting_smiles)
        if not mol:
            return [starting_smiles]  # Return original if invalid
        
        optimized_molecules = [starting_smiles]
        
        if not self.model_loaded:
            # Fall back to simplified optimization
            return self._optimize_simplified(starting_smiles, 4)
        
        try:
            if self.model_type == "transformer":
                # Use the transformer model for optimization
                # In a real implementation, this would condition the generation on the target properties
                prompt = starting_smiles[:len(starting_smiles)//2]  # Use part of the original SMILES as prompt
                
                for _ in range(4):
                    inputs = self.tokenizer(prompt, return_tensors="pt").to(self.device)
                    
                    # Generate new tokens
                    outputs = self.model.generate(
                        inputs.input_ids,
                        max_length=100,
                        do_sample=True,
                        top_p=0.9,
                        temperature=1.2,
                        num_return_sequences=1
                    )
                    
                    # Decode the generated tokens to text
                    generated_text = self.tokenizer.decode(outputs[0], skip_special_tokens=True)
                    
                    # Clean up the generated SMILES and validate
                    smiles = generated_text.strip()
                    new_mol = Chem.MolFromSmiles(smiles)
                    
                    if new_mol:
                        canonical_smiles = Chem.MolToSmiles(new_mol)
                        if canonical_smiles not in optimized_molecules:
                            optimized_molecules.append(canonical_smiles)
            
            elif self.model_type == "vae":
                # Use the VAE model for optimization
                # Convert SMILES to fingerprint
                fp = self._smiles_to_fingerprint(starting_smiles)
                if fp:
                    fp_tensor = self._fingerprint_to_tensor(fp)
                    
                    # Encode to latent space
                    with torch.no_grad():
                        mu, log_var = self.vae_model.encode(fp_tensor.unsqueeze(0))
                        z = self.vae_model.reparameterize(mu, log_var)
                        
                        # Generate variations by perturbing the latent vector
                        for _ in range(4):
                            # Add small random perturbation
                            perturbed_z = z + torch.randn_like(z) * 0.1
                            
                            # Decode
                            decoded = self.vae_model.decode(perturbed_z)
                            
                            # In a real implementation, you would convert the decoded tensor
                            # back to a valid SMILES string
                            # For demonstration, we'll use the simplified method
                            new_smiles = self._optimize_simplified(starting_smiles, 1)[0]
                            if new_smiles not in optimized_molecules:
                                optimized_molecules.append(new_smiles)
            
            elif self.model_type == "gan":
                # Use the GAN model for optimization
                gan_optimized = self._optimize_with_gan(starting_smiles, 4)
                for mol in gan_optimized:
                    if mol not in optimized_molecules:
                        optimized_molecules.append(mol)
            
            else:
                # Fall back to simplified optimization
                additional_molecules = self._optimize_simplified(starting_smiles, 4)
                for mol in additional_molecules:
                    if mol not in optimized_molecules:
                        optimized_molecules.append(mol)
        
        except Exception as e:
            print(f"Error in optimization: {e}")
            # Fall back to simplified optimization
            additional_molecules = self._optimize_simplified(starting_smiles, 4)
            for mol in additional_molecules:
                if mol not in optimized_molecules:
                    optimized_molecules.append(mol)
        
        return optimized_molecules
    
    def _optimize_simplified(self, starting_smiles: str, num_variations: int) -> List[str]:
        """Simplified molecule optimization as a fallback."""
        mol = Chem.MolFromSmiles(starting_smiles)
        if not mol:
            return [starting_smiles]
            
        optimized_molecules = []
        
        # Generate some variations
        for _ in range(num_variations):
            # Make a copy of the molecule
            new_mol = Chem.Mol(mol)
            
            # Randomly modify the molecule
            fragments = ["F", "Cl", "Br", "OH", "NH2", "CH3"]
            fragment = random.choice(fragments)
            frag_mol = Chem.MolFromSmiles(fragment)
            
            if frag_mol:
                try:
                    combo = Chem.CombineMols(new_mol, frag_mol)
                    new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(combo))
                    if new_mol:
                        optimized_molecules.append(Chem.MolToSmiles(new_mol))
                except:
                    pass
        
        return optimized_molecules
    
    def _optimize_with_gan(self, starting_smiles, num_variations=4):
        """Optimize a molecule using the GAN model."""
        try:
            if hasattr(self, 'gan_model'):
                # For a real implementation, you would:
                # 1. Convert the starting molecule to a fingerprint
                # 2. Encode it in the GAN's latent space
                # 3. Apply controlled perturbations to generate variations
                # 4. Decode back to SMILES
                
                # For now, we'll use the simplified generation
                optimized_molecules = []
                
                # Generate some base molecules
                base_molecules = self.gan_model.generate_molecules(num_variations)
                
                # Ensure the starting molecule is included
                if starting_smiles not in base_molecules:
                    optimized_molecules = [starting_smiles] + base_molecules
                else:
                    optimized_molecules = base_molecules
                
                return optimized_molecules[:num_variations+1]
            else:
                print("GAN model not initialized")
                return [starting_smiles] + self._optimize_simplified(starting_smiles, num_variations)
        except Exception as e:
            print(f"Error in GAN optimization: {e}")
            return [starting_smiles] + self._optimize_simplified(starting_smiles, num_variations) 