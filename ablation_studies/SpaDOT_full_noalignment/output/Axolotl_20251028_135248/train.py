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
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=UserWarning)

import sys
sys.path.append('/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies')
from model import SpaDOT_SVGPGATVAE_alignment
from utils import _train_utils, _utils, _analyze_utils

def train_SpaDOT(dataloader_dict, model_config):
    SpaDOT_model = SpaDOT_SVGPGATVAE_alignment.SpaDOT_SVGPGATVAE_alignment(model_config, dataloader_dict).to(model_config['device'])
    print(SpaDOT_model)
    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, SpaDOT_model.parameters()), lr=model_config['lr'])

    # --- prepare training
    beta1s = _train_utils._beta_cycle_linear(model_config['maxiter'], stop=model_config['beta1'])
    beta2s = _train_utils._beta_cycle_linear(model_config['maxiter'], stop=model_config['beta2'])
    tp_indexed_list = list(enumerate(model_config['timepoints']))
    # minimum_dataloader_length = min([len(dataloader_dict['dataloaders'][tp]) for tp in model_config['timepoints']])
    
    # --- initialize dictionary to store loss
    loss_dict = OrderedDict()
    loss_names = ['elbo', 'Recon', 'SVGP_KL', 'GAT_KL', 'alignment', 'KMeans', 'OT']
    for epoch in range(model_config['maxiter']):
        loss_dict[epoch] = OrderedDict()
        for name in loss_names:
            loss_dict[epoch][name] = 0

    print("Training SpaDOT model...")
    train_starttime = time()

    # --- add EMA early stopping
    ema_alpha = 0.9
    ema_recon = None
    best_ema = float('inf'); patience = 3; wait = 0

    for epoch in tqdm(range(model_config['maxiter'])):
        beta1 = beta1s[epoch]
        beta2 = beta2s[epoch]
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
            for batch_idx, batch in enumerate(tp_dataloader):
                # y: gene expression; x: location; tp_ix: timepoint index; edge_index: graph edge index
                y_batch, x_batch, tp_ix, edge_index_batch = batch.x, batch.loc, batch.data_index, batch.edge_index
                x_batch, y_batch, edge_index_batch = x_batch.to(model_config['device']), y_batch.to(model_config['device']), edge_index_batch.to(model_config['device'])
                tp_ix, adj_idx  = tp_ix[:batch.batch_size], batch.n_id[:batch.batch_size] # subset to seed nodes only
                tp_Recon_val, tp_SVGP_KL_val, tp_GAT_KL_val, tp_alignment_val, tp_p_m = \
                SpaDOT_model.forward(x=x_batch, y=y_batch, edge_index=edge_index_batch, 
                                     tp=tp, batch_size=batch.batch_size)
                # --- calculate KMeans loss from latent space
                tp_KMeans_val = torch.tensor(0, dtype=model_config['dtype'], device=model_config['device'])
                if epoch >= 1: # add Kmeans loss after the first epoch once the latent space is generated
                    kmeans_loss = _train_utils._compute_kmeans_loss(SpaDOT_model, model_config, tp, tp_ix, tp_p_m)
                    tp_KMeans_val += kmeans_loss
                # --- calculate OT loss from latent space
                tp_OT_val = torch.tensor(0, dtype=model_config['dtype'], device=model_config['device'])
                if epoch >= model_config['ot_epoch']:
                    if tp_i != 0:
                        tp_loss = _train_utils._compute_OT_loss(SpaDOT_model, model_config, tp, tp_ix, tp_p_m, model_config['timepoints'][tp_i-1])
                        tp_OT_val += tp_loss
                tp_elbo_val = model_config['lambda1'] * tp_Recon_val - \
                    beta1 * tp_SVGP_KL_val + \
                    model_config['beta2'] * tp_GAT_KL_val + \
                    model_config['omiga1'] * tp_alignment_val + \
                    model_config['omiga2'] * tp_KMeans_val + \
                    model_config['omiga3'] * tp_OT_val
                # beta2 * tp_GAT_KL_val + \
                # --- backward
                optimizer.zero_grad()
                tp_elbo_val.backward(retain_graph=True)
                torch.nn.utils.clip_grad_norm_(SpaDOT_model.parameters(), 0.5) # stabilize training
                optimizer.step()
                # --- update loss in each timepoint
                for name in loss_names:
                    tp_loss_dict[tp][name] += locals()[f'tp_{name}_val'].detach().cpu().item()
            # --- update loss in loss dictionary
            for name in loss_names:
                tp_loss_dict[tp][name] /= len(tp_dataloader)
                loss_dict[epoch][f'{name}'] += tp_loss_dict[tp][name]
        
        # update ema recon (use epoch recon)
        recon_epoch = loss_dict[epoch]['Recon']
        ema_recon = recon_epoch if ema_recon is None else ema_alpha * ema_recon + (1-ema_alpha) * recon_epoch
        if ema_recon < best_ema - 0.1:
            best_ema = ema_recon
            wait = 0
            best_state = SpaDOT_model.state_dict()
        else:
            wait += 1

        if wait >= patience:
            print(f"Early stopping at epoch {epoch}, best ema recon {best_ema:.6f}")
            SpaDOT_model.load_state_dict(best_state)
            break
        print(f'Best Recon: {best_ema}, Current EMA Recon: {ema_recon}')
        if epoch % 10 == 0:
            print(f'Epoch {epoch + 1}: Training time: {int(time()-ep_starttime)} seconds, ELBO: {loss_dict[epoch]["elbo"]:.8f}, Recon loss: {loss_dict[epoch]["Recon"]:.8f}, SVGP KL loss: {loss_dict[epoch]["SVGP_KL"]:.8f}, GAT KL loss: {loss_dict[epoch]["GAT_KL"]:.8f}, Alignment loss: {loss_dict[epoch]["alignment"]:.8f}, Kmeans loss: {loss_dict[epoch]["KMeans"]:.8f}, OT loss: {loss_dict[epoch]["OT"]:.8f}')
            
        # --- update Kmeans centroids
        _train_utils._update_Kmeans(SpaDOT_model, model_config, dataloader_dict)
        # --- update OT matrix 
        if (epoch + 1) % model_config['ot_config']['ot_epochs'] == 0:
            _train_utils._update_OT_matrix(SpaDOT_model, model_config)

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
    # --- write to record the config yaml
    import yaml
    with open(args.output_dir+os.sep+'config.yaml', 'w') as f:
        yaml.safe_dump(model_config, f)

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
    loss_df.T.to_csv(args.output_dir+os.sep+'loss.csv')
    if args.save_model:
        torch.save(SpaDOT_model.state_dict(), args.output_dir+os.sep+'SpaDOT_model.pth')
        print("Model saved to %s" % (args.output_dir))
    # obtain the latent representation
    latent_adata = get_latent(SpaDOT_model, model_config, adata, dataloader_dict)
    latent_adata.write_h5ad(args.output_dir+os.sep+'latent.h5ad')

def analyze(args):
    data_dir = os.path.abspath(args.data)
    if not args.output_dir:
        args.output_dir = os.path.dirname(data_dir)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if args.n_clusters is None and args.prefix == "":
        args.prefix = 'adaptive_'

    latent = anndata.read_h5ad(data_dir)
    if args.n_clusters is None:
        latent = _analyze_utils.Adaptive_clustering(args, latent)
    else:
        latent = _analyze_utils.KMeans_Clustering(latent, args.n_clusters)
    latent.obs['pixel_x'] = latent.obsm['spatial'][:, 0]
    latent.obs['pixel_y'] = latent.obsm['spatial'][:, 1]
    # draw domains
    _analyze_utils.plot_domains(args, latent)
    # perform optimal transport analysis
    _analyze_utils.OT_analysis(args, latent)
    # plot OT results
    _analyze_utils.plot_OT(args, latent)


if __name__ == "__main__":
    import argparse
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    parser = argparse.ArgumentParser(description="SpaDOT ablation studies.")
    parser.add_argument('--data', dest='data', type=str)
    parser.add_argument('--output_dir', dest='output_dir', type=str)
    parser.add_argument('--prefix', dest='prefix', type=str)
    args = parser.parse_args()
    args.config = '/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies/config.yaml'
    args.device = 'cuda:2'
    args.save_model = True
    print(args)
    # create output directory if not exists
    args.output_dir = os.path.join(args.output_dir, f"{args.prefix}_{timestamp}")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # --- copy folder to output
    import shutil
    import glob
    ablation_dir = '/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies'
    shutil.copytree(os.path.join(ablation_dir+os.sep+'model'), args.output_dir, dirs_exist_ok=True)
    shutil.copytree(os.path.join(ablation_dir+os.sep+'utils'), args.output_dir, dirs_exist_ok=True)
    shutil.copy(os.path.abspath(__file__), args.output_dir)

    # train SpaDOT model
    train(args)
    # analyze results
    args.data = f'{args.output_dir}/latent.h5ad'

    import anndata
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score
    latent_adata = anndata.read_h5ad(args.data)
    tps = list(latent_adata.obs['timepoint'].unique())
    ARI_dict = dict()
    n_clusters = []
    for i, tp in enumerate(tps):
        tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp].copy()
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), 
                           random_state=1993, 
                           n_init=10).fit(tp_adata.X)
        n_clusters.append(len(set(tp_adata.obs['annotation'])))
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        tp_ARI = adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])
        ARI_dict[tp] = tp_ARI
    pd.DataFrame(ARI_dict, index=[0]).to_csv(args.output_dir+os.sep+f'spaDOT_ARI.csv', index=False)
    args.n_clusters = n_clusters
    analyze(args)


