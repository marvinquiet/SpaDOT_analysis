import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import TransformerConv, GATConv
from torch_geometric.nn import SAGEConv


class SVGPEncoder(nn.Module):
    def __init__(self, input_dim, SVGP_z_dim, hidden_dims):
        super(SVGPEncoder, self).__init__()
        # Create a sequential network with specified hidden dimensions
        layers = [input_dim] + hidden_dims
        SVGP_encoder_net = []
        for i in range(1, len(layers)):
            linear_layer = nn.Linear(layers[i-1], layers[i])
            nn.init.xavier_uniform_(linear_layer.weight)
            SVGP_encoder_net.append(linear_layer)
            # normalization
            SVGP_encoder_net.append(nn.BatchNorm1d(layers[i]))
            # activation
            SVGP_encoder_net.append(nn.LeakyReLU())
        self.SVGP_encoder_net = nn.Sequential(*SVGP_encoder_net)
        self.SVGP_fc = nn.Linear(hidden_dims[-1], SVGP_z_dim*2)
        nn.init.xavier_uniform_(self.SVGP_fc.weight)


    def forward(self, x):
        '''
        x: gene expression
        '''
        h = self.SVGP_encoder_net(x)
        # Compute mean and log variance for the latent space
        SVGP_z = self.SVGP_fc(h)
        SVGP_enc_mu, SVGP_enc_logvar = torch.chunk(SVGP_z, 2, dim=1)
        return SVGP_enc_mu, torch.exp(SVGP_enc_logvar)  # Return mean and variance


class GATEncoder(nn.Module):
    def __init__(self, input_dim, GAT_z_dim, hidden_dim=512, num_heads=4):
        super(GATEncoder, self).__init__()
        # Define multiple GAT layers with specified dimensions and heads
        self.gat1 = TransformerConv(input_dim, hidden_dim, heads=num_heads, concat=False, dropout=0.2)
        # self.norm1 = nn.LayerNorm(hidden_dim * num_heads)
        self.gat2 = TransformerConv(hidden_dim, hidden_dim, heads=num_heads, concat=False, dropout=0.2)
        # self.norm2 = nn.LayerNorm(hidden_dim * num_heads)
        self.gat3 = TransformerConv(hidden_dim, hidden_dim, heads=num_heads, concat=False, dropout=0.2)
        # self.norm3 = nn.LayerNorm(hidden_dim)
        # self.gat1 = TransformerConv(input_dim, hidden_dim, heads=num_heads, concat=False, dropout=0.2)
        # self.gat2 = TransformerConv(hidden_dim, hidden_dim, heads=num_heads, concat=False, dropout=0.2)
        # self.gat3 = TransformerConv(hidden_dim, hidden_dim, heads=num_heads, concat=False, dropout=0.2)
        # self.norm = nn.LayerNorm(hidden_dim)

        # self.gat1 = GATConv(input_dim, hidden_dim, heads=num_heads, concat=True)
        # torch.nn.init.xavier_uniform(self.gat1.lin.weight)
        # self.gat2 = GATConv(hidden_dim * num_heads, hidden_dim, heads=num_heads, concat=True)
        # torch.nn.init.xavier_uniform(self.gat2.lin.weight)
        # self.gat3 = GATConv(hidden_dim * num_heads, hidden_dim, heads=num_heads, concat=False)
        # torch.nn.init.xavier_uniform(self.gat3.lin.weight)
        self.GAT_fc = nn.Linear(hidden_dim, GAT_z_dim*2)
        nn.init.xavier_uniform_(self.GAT_fc.weight)

    def forward(self, x, edge_index):
        '''
        x: gene expression
        edge_index: graph structure
        '''
        # h = F.leaky_relu(self.norm(self.gat1(x, edge_index)))
        # h = F.leaky_relu(self.norm(self.gat2(h, edge_index)))
        # h = self.norm(self.gat3(h, edge_index))
        h = F.leaky_relu(self.gat1(x, edge_index))
        h = F.leaky_relu(self.gat2(h, edge_index))
        h = self.gat3(h, edge_index)
        GAT_z = self.GAT_fc(h)
        GAT_enc_mu, GAT_enc_logvar = torch.chunk(GAT_z, 2, dim=1)
        return GAT_enc_mu, torch.exp(GAT_enc_logvar)


# class GATEncoder_norm(nn.Module):
#     def __init__(self, input_dim, GAT_z_dim, hidden_dim=512, num_heads=4):
#         super(GATEncoder_norm, self).__init__()
#         # Define multiple GAT layers with specified dimensions and heads
#         self.gat1 = TransformerConv(input_dim, hidden_dim, heads=num_heads, concat=True)
#         self.norm1 = nn.LayerNorm(hidden_dim * num_heads)
#         self.gat2 = TransformerConv(hidden_dim * num_heads, hidden_dim, heads=num_heads, concat=True)
#         self.norm2 = nn.LayerNorm(hidden_dim * num_heads)
#         self.gat3 = TransformerConv(hidden_dim * num_heads, hidden_dim, heads=num_heads, concat=False)
#         self.norm3 = nn.LayerNorm(hidden_dim)
#         self.GAT_fc = nn.Linear(hidden_dim, GAT_z_dim*2)
#         nn.init.xavier_uniform_(self.GAT_fc.weight)

#     def forward(self, x, edge_index):
#         '''
#         x: gene expression
#         edge_index: graph structure
#         '''
#         h = F.leaky_relu(self.norm1(self.gat1(x, edge_index)))
#         h = F.leaky_relu(self.norm2(self.gat2(h, edge_index)))
#         h = self.norm3(self.gat3(h, edge_index))
#         GAT_z = self.GAT_fc(h)
#         GAT_enc_mu, GAT_enc_logvar = torch.chunk(GAT_z, 2, dim=1)
#         return GAT_enc_mu, torch.exp(GAT_enc_logvar)


class GCNEncoder(nn.Module):
    def __init__(self, input_dim, GAT_z_dim, hidden_dim=512):
        super(GCNEncoder, self).__init__()
        self.gcn1 = SAGEConv(input_dim, hidden_dim)
        self.gcn2 = SAGEConv(hidden_dim, hidden_dim)
        self.gcn3 = SAGEConv(hidden_dim, hidden_dim)
        self.GCN_fc = nn.Linear(hidden_dim, GAT_z_dim*2)
        nn.init.xavier_uniform_(self.GCN_fc.weight)

    def forward(self, x, edge_index):
        '''
        x: gene expression
        edge_index: graph structure
        '''
        h = F.leaky_relu(self.gcn1(x, edge_index))
        h = F.leaky_relu(self.gcn2(h, edge_index))
        h = self.gcn3(h, edge_index)
        GCN_z = self.GCN_fc(h)
        GCN_enc_mu, GCN_enc_logvar = torch.chunk(GCN_z, 2, dim=1)
        return GCN_enc_mu, torch.exp(GCN_enc_logvar)

# --- define traditional VAE
class MLPEncoder(nn.Module): # same structure with MLP for testing backbone
    def __init__(self, input_dim, MLP_z_dim, hidden_dims):
        super(MLPEncoder, self).__init__()
        # Create a sequential network with specified hidden dimensions
        layers = [input_dim] + hidden_dims
        MLP_encoder_net = []
        for i in range(1, len(layers)):
            linear_layer = nn.Linear(layers[i-1], layers[i])
            nn.init.xavier_uniform_(linear_layer.weight)
            MLP_encoder_net.append(linear_layer)
            # normalization
            MLP_encoder_net.append(nn.BatchNorm1d(layers[i]))
            # activation
            MLP_encoder_net.append(nn.LeakyReLU())
        self.MLP_encoder_net = nn.Sequential(*MLP_encoder_net)
        self.MLP_fc = nn.Linear(hidden_dims[-1], MLP_z_dim*2)
        nn.init.xavier_uniform_(self.MLP_fc.weight)

    def forward(self, x):
        '''
        x: gene expression
        '''
        h = self.MLP_encoder_net(x)
        # Compute mean and log variance for the latent space
        MLP_z = self.MLP_fc(h)
        MLP_enc_mu, MLP_enc_logvar = torch.chunk(MLP_z, 2, dim=1)
        return MLP_enc_mu, torch.exp(MLP_enc_logvar)  # Return mean and variance
