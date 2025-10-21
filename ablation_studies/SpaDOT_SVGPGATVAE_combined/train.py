import os
import anndata
import pandas as pd
import scanpy as sc
import random
from time import time
import torch
from torch import optim
from tqdm.auto import tqdm
from collections import OrderedDict
from torch_geometric.utils import dense_to_sparse

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="torch") # this is for suppressing UserWarning from torch.tensor

import sys
sys.path.append('/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies')
from model import SpaDOT_SVGPGATVAE
from utils import _train_utils, _utils


def train_SpaDOT(dataloader_dict, model_config):
    SpaDOT_model = SpaDOT_SVGPGATVAE.SpaDOT_SVGPGATVAE(model_config, dataloader_dict).to(model_config['device'])
    print(SpaDOT_model)
    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, SpaDOT_model.parameters()), lr=model_config['lr'])

    # --- prepare training
    beta1s = _train_utils._beta_cycle_linear(model_config['maxiter'], stop=model_config['beta1'])
    tp_indexed_list = list(enumerate(model_config['timepoints']))
    
    # --- initialize dictionary to store loss
    loss_dict = OrderedDict()
    loss_names = ['elbo', 'Recon', 'SVGP_KL', 'GAT_KL']
    for epoch in range(model_config['maxiter']):
        loss_dict[epoch] = OrderedDict()
        for name in loss_names:
            loss_dict[epoch][name] = 0

    print("Training SpaDOT model...")
    train_starttime = time()
    for epoch in tqdm(range(model_config['maxiter'])):
        beta1 = beta1s[epoch]
        # --- start training
        SpaDOT_model.train()
        ep_starttime = time()
        tp_loss_dict = OrderedDict()
        # --- random shuffle timepoints
        random.shuffle(tp_indexed_list) 
        for tp_i, tp in tp_indexed_list:
            tp_loss_dict[tp] = OrderedDict()
            for name in loss_names:
                tp_loss_dict[tp][name] = 0
            tp_dataloader = dataloader_dict['dataloaders'][tp]
            for _, batch in enumerate(tp_dataloader):
                # y: gene expression; x: location; tp_ix: timepoint index; edge_index: graph edge index
                y_batch, x_batch, tp_ix, edge_index_batch = batch.x, batch.loc, batch.data_index, batch.edge_index
                x_batch, y_batch, edge_index_batch = x_batch.to(model_config['device']), y_batch.to(model_config['device']), edge_index_batch.to(model_config['device'])
                tp_ix, adj_idx,  = tp_ix[:batch.batch_size], batch.n_id[:batch.batch_size] # subset to seed nodes only
                tp_Recon_val, tp_SVGP_KL_val, tp_GAT_KL_val, tp_p_m = \
                SpaDOT_model.forward(x=x_batch, y=y_batch, edge_index=edge_index_batch, 
                                     tp=tp, batch_size=batch.batch_size)
                tp_elbo_val = model_config['lambda1'] * tp_Recon_val - \
                    beta1 * tp_SVGP_KL_val + \
                    model_config['beta2'] * tp_GAT_KL_val
                # --- backward
                optimizer.zero_grad()
                tp_elbo_val.backward(retain_graph=True)
                torch.nn.utils.clip_grad_norm_(SpaDOT_model.parameters(), 0.3) # stabilize training
                optimizer.step()
                # --- update loss in each timepoint
                for name in loss_names:
                    tp_loss_dict[tp][name] += locals()[f'tp_{name}_val'].detach().cpu().item()
            # --- update loss in loss dictionary
            for name in loss_names:
                tp_loss_dict[tp][name] /= len(tp_dataloader)
                loss_dict[epoch][f'{name}'] += tp_loss_dict[tp][name]
        if epoch % 10 == 0:
            print(f'Epoch {epoch + 1}: Training time: {int(time()-ep_starttime)} seconds, ELBO: {loss_dict[epoch]["elbo"]:.8f}, Recon loss: {loss_dict[epoch]["Recon"]:.8f}, SVGP KL loss: {loss_dict[epoch]["SVGP_KL"]:.8f}, GAT KL loss: {loss_dict[epoch]["GAT_KL"]:.8f}')

    print('Training finished...')
    print('Training time: %d seconds.' % int(time() - train_starttime))
    loss_df = pd.DataFrame.from_dict(loss_dict)
    return SpaDOT_model, loss_df

# --- obtain latent representations
def get_latent(model, model_config, adata, dataloader_dict):
    '''Obtain latent representations
    :param model: SpaDOT model
    :param model_config: model configuration
    :return: latent representations
    '''
    model.eval()
    with torch.no_grad():
        latent_adata_list = list()
        for tp in model_config['timepoints']:
            tp_loc, tp_y, tp_idx = dataloader_dict['datasets'][tp].dataset.tensors
            tp_adj = dataloader_dict['adjacency_matrices'][tp]
            tp_edge_index, _ = dense_to_sparse(tp_adj)
            latent_samples = model.all_latent_samples(tp_loc, tp_y, tp_edge_index, tp)  
            tp_latent_adata = sc.AnnData(latent_samples, obs=adata[tp_idx.numpy()].obs)
            tp_latent_adata.obsm['spatial'] = adata.obsm['spatial'][tp_idx.numpy()]
            latent_adata_list.append(tp_latent_adata)
        latent_adata = anndata.concat(latent_adata_list)
    return latent_adata


def train(args):
    # --- load data ---
    print("Loading data...")
    data_dir = os.path.abspath(args.data)
    if not args.output_dir:
        args.output_dir = os.path.dirname(data_dir)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    adata = anndata.read_h5ad(data_dir)
    model_config = _utils.load_model_config(args)
    # add adata related parameters to model_config
    model_config['input_dim'] = adata.n_vars
    tps = list(adata.obs['timepoint'].unique()) # need to be list
    tps.sort()
    model_config['timepoints'] = tps
    # add device and dtype to model_config
    model_config['device'] = torch.device(args.device)
    model_config['dtype'] = torch.float64

    _utils.set_seed(model_config['seed'])
    # obtain inducing points, graph adjacency matrix, and dataloaders
    print("Preparing data...")
    dataloader_dict = _train_utils.prepare_dataloader(adata, model_config)
    _utils._save_inducing_points(args, dataloader_dict['inducing_points']) # save inducing points

    # train model
    print("Training model...")
    SpaDOT_model, loss_df = train_SpaDOT(dataloader_dict, model_config)
    loss_df.T.to_csv(args.output_dir+os.sep+args.prefix+'loss.csv')
    if args.save_model:
        torch.save(SpaDOT_model.state_dict(), args.output_dir+os.sep+args.prefix+'SpaDOT_model.pth')
        print("Model saved to %s" % (args.output_dir))
    # obtain the latent representation
    latent_adata = get_latent(SpaDOT_model, model_config, adata, dataloader_dict)
    latent_adata.write_h5ad(args.output_dir+os.sep+args.prefix+'latent.h5ad')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="SpaDOT ablation studies.")
    parser.add_argument('--data', dest='data', type=str)
    parser.add_argument('--output_dir', dest='output_dir', type=str)
    parser.add_argument('--prefix', dest='prefix', type=str)
    args = parser.parse_args()
    args.config = '/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies/config.yaml'
    args.device = 'cuda:3'
    args.save_model = True
    print(args)
    # create output directory if not exists
    if 'output_dir' not in args.__dict__:
        args.output_dir = os.path.dirname(args.data)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # train SpaDOT model
    train(args)
