import torch
import torch.nn as nn
import torch.nn.functional as F
from .encoder import MLPEncoder
from .decoder import Decoder

class SpaDOT_VAE(nn.Module):
    def __init__(self, model_config, dataloader_dict):
        super(SpaDOT_VAE, self).__init__()
        self.input_dim = model_config['input_dim']
        self.MLP_z_dim = model_config['z_dim']
        self.dtype = model_config['dtype']
        self.device = model_config['device']

        # Build encoders
        self.MLPEncoder = MLPEncoder(
            input_dim=self.input_dim,
            MLP_z_dim=self.MLP_z_dim,
            hidden_dims=model_config['svgp_encoder_layers']
        ).to(dtype=self.dtype)

        # Build decoder
        self.decoder = Decoder(
            input_dim=self.input_dim,
            z_dim=self.MLP_z_dim, # in case input z_dim is odd
            decoder_layers=model_config['decoder_layers']
        ).to(dtype=self.dtype)

    def forward(self, x, y, tp):
        MLP_m, MLP_v = self.MLPEncoder(y)
        MLP_latent_sample = MLP_m + torch.randn_like(MLP_m) * torch.sqrt(MLP_v)
        final_latent = MLP_latent_sample
        # KL loss
        MLP_KL = -0.5 * torch.sum(1 + torch.log(MLP_v) - MLP_m.pow(2) - MLP_v) / self.MLP_z_dim
        # Reconstruction loss
        recon_loss = torch.sum((y-self.decoder(final_latent))**2) / self.input_dim
        return recon_loss, MLP_KL, final_latent
    
    def all_latent_samples(self, X, Y, tp):
        """
        Output latent embedding.

        Parameters:
        -----------
        X: array_like, shape (n_spots, 2)
            Location information.
        Y: array_like, shape (n_spots, n_genes)
            Preprocessed count matrix.
        """ 
        # solve warinings
        X = torch.tensor(X, dtype=self.dtype, device=self.device)
        Y = torch.tensor(Y, dtype=self.dtype, device=self.device)
        MLP_m, _ = self.MLPEncoder(Y)
        latent_samples = MLP_m.data.cpu().detach().numpy()
        return latent_samples
    
    def _gauss_cross_entropy(self, mu1, var1, mu2, var2):
        """
        Computes the element-wise cross entropy
        Given q(z) ~ N(z| mu1, var1)
        returns E_q[ log N(z| mu2, var2) ]
        args:
            mu1:  mean of expectation (batch, tmax, 2) tf variable
            var1: var  of expectation (batch, tmax, 2) tf variable
            mu2:  mean of integrand (batch, tmax, 2) tf variable
            var2: var of integrand (batch, tmax, 2) tf variable
        returns:
            cross_entropy: (batch, tmax, 2) tf variable
        """
        term0 = 1.8378770664093453  # log(2*pi)
        term1 = torch.log(var2)
        term2 = (var1 + mu1 ** 2 - 2 * mu1 * mu2 + mu2 ** 2) / var2
        cross_entropy = -0.5 * (term0 + term1 + term2)
        return cross_entropy
