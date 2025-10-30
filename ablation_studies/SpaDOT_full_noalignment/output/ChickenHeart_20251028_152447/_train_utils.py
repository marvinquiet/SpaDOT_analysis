import os, math
from time import time
import random
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

import torch
from torch.utils.data import Dataset, TensorDataset
from torch_geometric.data import Data
from torch_geometric.utils import dense_to_sparse
from torch_geometric.loader import NeighborLoader
from collections import OrderedDict

# --- my package
from utils import _utils
from utils.OT_loss.ot_solvers import compute_transport_map # Waddington, faster C version

class MyDataset(Dataset):
    def __init__(self, pos, counts, tp_ix, dtype=torch.float64):
        self.dataset = TensorDataset(torch.tensor(pos, dtype=dtype), 
                                     torch.tensor(counts, dtype=dtype), 
                                     torch.tensor(tp_ix, dtype=torch.int))
    def __getitem__(self, index):
        pos, counts, tp_ix  = self.dataset[index]
        return pos, counts, tp_ix
    def __len__(self):
        return len(self.dataset)

def prepare_dataloader(adata, model_config):
    # --- obtain location and time point information
    loc = _obtain_tp_loc_info(adata)
    # sample some inducing points
    inducing_idx = random.sample(range(loc.shape[0]), model_config['inducing_point_nums']) 
    inducing_points = loc[inducing_idx, :]
    # --- generate mapper between tp_mat index and real timepoint -> in case the timepoint is not numeric or continous
    idx_to_tp = dict(enumerate(model_config['timepoints']))
    tp_to_idx = {v:k for k,v in idx_to_tp.items()} 
    idx = np.argmax(loc[:, 2:], axis=1)    
    tp_list = [idx_to_tp[_] for _ in idx]
    tp_index_dict = dict() # tp index in timepoint list
    for tp in model_config['timepoints']:
        tp_ix = [_ for _ in range(len(tp_list)) if tp_list[_] == tp]
        tp_index_dict[tp] = tp_ix
    # create inducing point dict / number of training sample dict for separate decoder
    inducing_points_dict, N_train_dict = OrderedDict(), OrderedDict()
    inducing_points_tp = np.argmax(inducing_points[:, 2:], axis=1)
    for _, tp in enumerate(model_config['timepoints']):
        idx = tp_to_idx[tp]
        tp_inducing_idx = np.where(inducing_points_tp == idx)[0]
        inducing_points_dict[tp] = inducing_points[tp_inducing_idx, :2]
        N_train_dict[tp] = np.sum(adata.obs['timepoint'] == tp) # number of training data in each timepoint
    # --- generate timepoint-specific dataloader
    dataloaders_dict, adj_dict, dataset_dict = OrderedDict(), OrderedDict(), OrderedDict()
    for tp in model_config['timepoints']:
        tp_ix = tp_index_dict[tp]
        tp_adata = adata[tp_ix].copy()
        tp_loc = loc[tp_ix, :2]
        tp_dataset = MyDataset(tp_loc, tp_adata.X, tp_ix)
        dataset_dict[tp] = tp_dataset
        # --- Neighborloader for graph data
        k_cutoff_factor =  min(30, model_config['knn_cutoff'] * (1 if tp_adata.n_obs < 3000 else 1 + round((tp_adata.n_obs - 3000) / 1000)))
        # k_cutoff_factor = min(30, model_config['knn_cutoff']*round(1/1000*tp_adata.n_obs))
        _utils._Cal_Spatial_Net(tp_adata, 
                            #    k_cutoff=min(model_config['max_neighbors'], model_config['knn_cutoff']*round(1/1000*tp_adata.n_obs)),
                               k_cutoff=min(model_config['max_neighbors'], k_cutoff_factor),
                               max_neigh=model_config['max_neighbors'])
        tp_adj = torch.tensor(tp_adata.uns['adj'], dtype=model_config['dtype'])
        tp_edge_index, _ = dense_to_sparse(tp_adj)
        tp_graph_data = Data(
            x=torch.tensor(tp_adata.X, dtype=model_config['dtype']),
            # edge_index=tp_edge_index.clone().to(dtype=torch.long), # torch copy warning
            edge_index=torch.tensor(tp_edge_index, dtype=torch.long),
            data_index=torch.tensor(tp_ix, dtype=torch.int),
            loc=torch.tensor(tp_loc, dtype=model_config['dtype'])
        )
        # tp_graph_loader = DataLoader(tp_graph_data, batch_size = model_config['batch_size'], shuffle=True)
        tp_graph_loader = NeighborLoader(
            tp_graph_data, 
            num_neighbors=[max(30, model_config['knn_cutoff']*round(1/1000*tp_adata.n_obs))] * 2, 
            batch_size = model_config['batch_size'],
            subgraph_type="induced",
            worker_init_fn=_utils.seed_worker,
            shuffle=True
        )
        dataloaders_dict[tp] = tp_graph_loader 
        adj_dict[tp] = tp_adj
    return {
        "inducing_points": inducing_points_dict,
        "N_train": N_train_dict,
        "dataloaders": dataloaders_dict,
        "adjacency_matrices": adj_dict,
        "datasets": dataset_dict
    }


def _obtain_tp_loc_info(adata):
    '''
    Obtain the location and time point information from the adata object.
    The location is standardized and concatenated with the time point information.
    The time point information is one-hot encoded.
    '''
    # in case the timepoint is not numeric
    adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
    tp_info = np.array(adata.obs['timepoint_numeric']).astype('int')
    tp_mat = np.zeros((tp_info.size, tp_info.max()+1))
    tp_mat[np.arange(tp_info.size), tp_info] = 1
    n_tp = tp_mat.shape[1]
    # scale locations per time point
    loc = adata.obsm['spatial']
    loc_scaled = np.zeros(loc.shape, dtype=np.float64)
    for i in range(n_tp):
        scaler = StandardScaler()
        tp_loc = loc[tp_mat[:,i]==1, :]
        tp_loc = scaler.fit_transform(tp_loc)
        loc_scaled[tp_mat[:,i]==1, :] = tp_loc
    loc = loc_scaled
    loc = np.concatenate((loc, tp_mat), axis=1)
    return loc


def _beta_cycle_linear(n_iter, start=0.0, stop=1, n_cycle=10, ratio=1):
    L = np.ones(n_iter) * stop
    period = n_iter/n_cycle
    step = (stop-start)/(period*ratio) # linear schedule
    for c in range(n_cycle):
        v, i = start, 0
        while v <= stop and (int(i+c*period) < n_iter):
            L[int(i+c*period)] = v
            v += step
            i += 1
    return L 

# K-Means related
def _compute_kmeans_loss(model, model_config, tp, tp_ix, tp_p_m):
    '''Compute KMeans loss between assignment and centroids
    :param tp: time point
    :param tp_ix: time point specific index
    :param tp_p_m: posterior mean

    :return: K-Means loss
    '''
    cluster_index_dict = model.kmeans_index_dict
    cluster_assignments = [cluster_index_dict[tp][_.item()] for _ in tp_ix]
    cluster_centers = torch.tensor(model.kmeans_center_dict[tp], dtype=model_config['dtype'], device=model_config['device'])
    # Compute K-means loss (mean of squared distances)
    kmeans_loss = torch.sum(torch.norm(tp_p_m - cluster_centers[cluster_assignments]) ** 2 / tp_p_m.shape[1] / len(set(cluster_assignments)) )
    return kmeans_loss

def _update_Kmeans(model, model_config, dataloader_dict):
    model.eval()
    with torch.no_grad():
        for tp in model_config['timepoints']:
            tp_loc, tp_y, tp_idx = dataloader_dict['datasets'][tp].dataset.tensors
            tp_adj = dataloader_dict['adjacency_matrices'][tp]
            tp_edge_index, _ = dense_to_sparse(tp_adj)
            latent_samples = model.all_latent_samples(tp_loc, tp_y, tp_edge_index, tp)  
            # adaptive KMeans
            tp_kmeans = KMeans(n_clusters=model_config['n_clusters'], 
                               random_state=model_config['seed'], 
                               n_init=10).fit(latent_samples)
            model.kmeans_center_dict[tp] = tp_kmeans.cluster_centers_
            model.kmeans_cluster_dict[tp] = tp_kmeans.labels_.tolist()
            model.kmeans_index_dict[tp] = dict(zip(tp_idx.numpy(), tp_kmeans.labels_))

# OT related
def _compute_OT_loss(model, model_config, cur_tp, tp_ix, tp_p_m, prev_tp):
    '''Compute OT loss between timepoints
    :param model:  model
    :param tp: current timepoints
    :param tp_ix: index of current tp to obtain hidden clusters
    :param tp_p_m: posterior mean of current timepoints
    :prev_tp: previous timepoint
    :return: OT transport cost
    '''
    cur_tp_cluster_latent = dict()
    cur_tp_clusters = [model.kmeans_index_dict[cur_tp][_] for _ in tp_ix.numpy()]
    cluster_list = list(set(model.kmeans_cluster_dict[cur_tp]))
    cluster_list.sort()
    for idx, cluster in enumerate(cur_tp_clusters):
        if cluster in cur_tp_cluster_latent:
            cur_tp_cluster_latent[cluster].append(tp_p_m[idx])
        else:
            cur_tp_cluster_latent[cluster] = [tp_p_m[idx]]
    cur_tp_cluster_center = []
    for cluster in cluster_list:
        if cluster not in cur_tp_cluster_latent:
            cur_tp_cluster_center.append(torch.tensor(model.kmeans_center_dict[cur_tp][cluster], 
                                                      dtype=model_config['dtype'], device=model_config['device']))
        else:
            cur_tp_cluster_center.append(torch.mean(torch.stack(cur_tp_cluster_latent[cluster]), dim=0))
    cur_tp_cluster_center = torch.stack(cur_tp_cluster_center)
    gamma = model.gammas[f"{prev_tp}_{cur_tp}"] # obtain transport plan
    gamma = gamma / gamma.sum(axis=1, keepdims=True) # normalize rows
    gamma = np.nan_to_num(gamma, nan=0.0, posinf=0.0, neginf=0.0)  # prune invalid values
    gamma = torch.tensor(gamma, dtype=model_config['dtype']).to(model_config['device'])
    cost_matrix = torch.cdist(
        torch.tensor(model.kmeans_center_dict[prev_tp], dtype=model_config['dtype']).to(model_config['device']), 
        cur_tp_cluster_center, 
        p=2) # (prev_tp_cluster, cur_tp_cluster)
    transport_cost = torch.mean(gamma*cost_matrix)
    return transport_cost

def _update_OT_matrix(model, model_config):
    '''Compute OT matrix
    '''
    model.eval()
    timepoints = model_config['timepoints']
    for tp_i, tp in enumerate(timepoints):
        if tp_i == len(timepoints) - 1: break
        cur_tp = tp
        next_tp = timepoints[tp_i+1]
        gamma = compute_transport_map(model.kmeans_center_dict[cur_tp],
                                       model.kmeans_center_dict[next_tp], 
                                       model_config['ot_config'], G=None)
        model.gammas[f"{timepoints[tp_i]}_{timepoints[tp_i+1]}"] = gamma