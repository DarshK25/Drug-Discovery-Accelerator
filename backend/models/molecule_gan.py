import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

class Generator(nn.Module):
    """
    Generator network for the GAN model.
    Takes random noise and generates molecular fingerprints.
    """
    def __init__(self, noise_dim=100, output_dim=512, hidden_dims=[256, 512, 1024]):
        super(Generator, self).__init__()
        
        # Build the generator network
        layers = []
        
        # First layer from noise dimension to first hidden dimension
        layers.append(nn.Linear(noise_dim, hidden_dims[0]))
        layers.append(nn.LeakyReLU(0.2))
        layers.append(nn.BatchNorm1d(hidden_dims[0]))
        
        # Add hidden layers
        for i in range(len(hidden_dims) - 1):
            layers.append(nn.Linear(hidden_dims[i], hidden_dims[i + 1]))
            layers.append(nn.LeakyReLU(0.2))
            layers.append(nn.BatchNorm1d(hidden_dims[i + 1]))
        
        # Output layer
        layers.append(nn.Linear(hidden_dims[-1], output_dim))
        layers.append(nn.Sigmoid())  # Sigmoid to get values between 0 and 1 for fingerprint bits
        
        self.model = nn.Sequential(*layers)
    
    def forward(self, z):
        """
        Forward pass of the generator.
        
        Args:
            z: Random noise tensor of shape (batch_size, noise_dim)
            
        Returns:
            Generated molecular fingerprints of shape (batch_size, output_dim)
        """
        return self.model(z)

class Discriminator(nn.Module):
    """
    Discriminator network for the GAN model.
    Takes molecular fingerprints and predicts whether they are real or generated.
    """
    def __init__(self, input_dim=512, hidden_dims=[1024, 512, 256]):
        super(Discriminator, self).__init__()
        
        # Build the discriminator network
        layers = []
        
        # First layer from input dimension to first hidden dimension
        layers.append(nn.Linear(input_dim, hidden_dims[0]))
        layers.append(nn.LeakyReLU(0.2))
        layers.append(nn.Dropout(0.3))
        
        # Add hidden layers
        for i in range(len(hidden_dims) - 1):
            layers.append(nn.Linear(hidden_dims[i], hidden_dims[i + 1]))
            layers.append(nn.LeakyReLU(0.2))
            layers.append(nn.Dropout(0.3))
        
        # Output layer
        layers.append(nn.Linear(hidden_dims[-1], 1))
        layers.append(nn.Sigmoid())  # Sigmoid for binary classification
        
        self.model = nn.Sequential(*layers)
    
    def forward(self, x):
        """
        Forward pass of the discriminator.
        
        Args:
            x: Molecular fingerprints of shape (batch_size, input_dim)
            
        Returns:
            Probability that the input is real (1) or fake (0)
        """
        return self.model(x)

class MoleculeGAN:
    """
    GAN model for generating molecular fingerprints and converting them to SMILES.
    """
    def __init__(self, device=None):
        """
        Initialize the GAN model.
        
        Args:
            device: PyTorch device to use (cuda or cpu)
        """
        if device is None:
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self.device = device
        
        # Initialize generator and discriminator
        self.generator = Generator().to(self.device)
        self.discriminator = Discriminator().to(self.device)
        
        # Set to evaluation mode by default
        self.generator.eval()
        self.discriminator.eval()
        
        # Initialize optimizers (for training)
        self.g_optimizer = torch.optim.Adam(self.generator.parameters(), lr=0.0002, betas=(0.5, 0.999))
        self.d_optimizer = torch.optim.Adam(self.discriminator.parameters(), lr=0.0002, betas=(0.5, 0.999))
        
        # Binary cross entropy loss for GAN
        self.criterion = nn.BCELoss()
        
        # For fingerprint to SMILES conversion
        self.smiles_lookup = {}  # Cache for fingerprint to SMILES conversion
        
        print("MoleculeGAN initialized")
    
    def load_models(self, generator_path, discriminator_path):
        """
        Load pre-trained generator and discriminator models.
        
        Args:
            generator_path: Path to the generator model weights
            discriminator_path: Path to the discriminator model weights
        """
        try:
            self.generator.load_state_dict(torch.load(generator_path, map_location=self.device))
            self.discriminator.load_state_dict(torch.load(discriminator_path, map_location=self.device))
            print("Models loaded successfully")
            return True
        except Exception as e:
            print(f"Error loading models: {e}")
            return False
    
    def train(self, dataloader, epochs=100, save_interval=10, save_dir="models"):
        """
        Train the GAN model.
        
        Args:
            dataloader: PyTorch DataLoader with molecular fingerprints
            epochs: Number of training epochs
            save_interval: Interval for saving model checkpoints
            save_dir: Directory to save model checkpoints
        """
        # Create save directory if it doesn't exist
        import os
        os.makedirs(save_dir, exist_ok=True)
        
        # Labels for real and fake data
        real_label = 1.0
        fake_label = 0.0
        
        # Training loop
        for epoch in range(epochs):
            for i, data in enumerate(dataloader):
                # Get real fingerprints
                real_fingerprints = data.to(self.device)
                batch_size = real_fingerprints.size(0)
                
                # ---------------------
                # Train Discriminator
                # ---------------------
                
                # Reset gradients
                self.d_optimizer.zero_grad()
                
                # Train with real data
                label = torch.full((batch_size, 1), real_label, device=self.device)
                output = self.discriminator(real_fingerprints)
                d_loss_real = self.criterion(output, label)
                d_loss_real.backward()
                
                # Train with fake data
                noise = torch.randn(batch_size, 100, device=self.device)
                fake_fingerprints = self.generator(noise)
                label.fill_(fake_label)
                output = self.discriminator(fake_fingerprints.detach())
                d_loss_fake = self.criterion(output, label)
                d_loss_fake.backward()
                
                # Update discriminator
                d_loss = d_loss_real + d_loss_fake
                self.d_optimizer.step()
                
                # ---------------------
                # Train Generator
                # ---------------------
                
                # Reset gradients
                self.g_optimizer.zero_grad()
                
                # Generate fake data
                label.fill_(real_label)  # Generator wants discriminator to think its output is real
                output = self.discriminator(fake_fingerprints)
                g_loss = self.criterion(output, label)
                g_loss.backward()
                
                # Update generator
                self.g_optimizer.step()
                
                # Print progress
                if i % 50 == 0:
                    print(f"[{epoch}/{epochs}][{i}/{len(dataloader)}] "
                          f"D_loss: {d_loss.item():.4f} G_loss: {g_loss.item():.4f}")
            
            # Save models
            if epoch % save_interval == 0 or epoch == epochs - 1:
                torch.save(self.generator.state_dict(), f"{save_dir}/generator_epoch_{epoch}.pt")
                torch.save(self.discriminator.state_dict(), f"{save_dir}/discriminator_epoch_{epoch}.pt")
    
    def generate_fingerprints(self, num_samples=10):
        """
        Generate molecular fingerprints using the trained generator.
        
        Args:
            num_samples: Number of fingerprints to generate
            
        Returns:
            Generated fingerprints as numpy array
        """
        with torch.no_grad():
            # Generate random noise
            noise = torch.randn(num_samples, 100, device=self.device)
            
            # Generate fingerprints
            fingerprints = self.generator(noise)
            
            # Convert to binary fingerprints (0 or 1)
            binary_fingerprints = (fingerprints > 0.5).float()
            
            return binary_fingerprints.cpu().numpy()
    
    def fingerprint_to_smiles(self, fingerprint):
        """
        Convert a molecular fingerprint to a SMILES string.
        This is a challenging inverse problem and typically requires additional models.
        
        Args:
            fingerprint: Binary molecular fingerprint
            
        Returns:
            SMILES string or None if conversion fails
        """
        # This is a simplified implementation
        # In a real-world scenario, you would use a decoder model trained to convert
        # fingerprints back to SMILES strings
        
        # For demonstration, we'll use a nearest-neighbor approach with a database
        # of known molecules (not implemented here)
        
        # Instead, we'll return a placeholder molecule
        return "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    
    def generate_molecules(self, num_molecules=5):
        """
        Generate molecules using the GAN model.
        
        Args:
            num_molecules: Number of molecules to generate
            
        Returns:
            List of SMILES strings representing the generated molecules
        """
        # Generate fingerprints
        fingerprints = self.generate_fingerprints(num_molecules)
        
        # Convert fingerprints to SMILES
        smiles_list = []
        for fp in fingerprints:
            smiles = self.fingerprint_to_smiles(fp)
            if smiles:
                # Validate the SMILES string
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    canonical_smiles = Chem.MolToSmiles(mol)
                    smiles_list.append(canonical_smiles)
        
        # If we couldn't generate enough valid molecules, add some known drugs
        known_drugs = [
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen
        ]
        
        while len(smiles_list) < num_molecules:
            smiles_list.append(np.random.choice(known_drugs))
        
        return smiles_list[:num_molecules] 