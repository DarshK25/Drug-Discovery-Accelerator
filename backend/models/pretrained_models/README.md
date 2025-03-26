# Pre-trained Models for Drug Discovery Assistant

This directory contains pre-trained models for the Drug Discovery Assistant's generative AI components.

## Models

### Transformer Models

The system uses ChemBERTa, a BERT model fine-tuned on chemical data, for molecule generation. The model is loaded from Hugging Face.

### Variational Autoencoder (VAE) Models

Place your pre-trained VAE model weights here with the filename `vae_model.pt`.

### Generative Adversarial Network (GAN) Models

Place your pre-trained GAN model weights here with the filenames:
- `generator.pt`: The generator model weights
- `discriminator.pt`: The discriminator model weights

## Training Your Own Models

### Training a VAE

To train your own VAE model:

1. Collect a dataset of SMILES strings
2. Convert them to molecular fingerprints
3. Train the VAE model on these fingerprints
4. Save the model weights to `vae_model.pt`

### Training a GAN

To train your own GAN model:

1. Collect a dataset of molecular fingerprints
2. Use the `MoleculeGAN.train()` method to train the model
3. The trained models will be saved to this directory

## Using Pre-trained Models

The system will automatically look for model weights in this directory when initializing the generative models. If no weights are found, it will fall back to simplified generation methods. 