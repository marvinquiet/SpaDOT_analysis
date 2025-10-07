# conda activate /media/sda/wenjinma/envs/STAligner_py38/
import warnings
warnings.filterwarnings("ignore")

import STAligner

import os
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import scipy.linalg
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

import torch
import random
from torch.backends import cudnn
seed = 1993
os.environ['PYTHONHASHSEED'] = str(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
random.seed(seed)
cudnn.deterministic = True
cudnn.benchmark = False

used_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_STAligner/MOSTA_brain_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'adata_brain_clean.h5ad')
n_samples = 20000  # Set desired number of cells
selected_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
adata = adata[selected_indices, :].copy()
adata.X = adata.layers['count'].astype(int)
adata.obs['annotation'] = 'Brain'
# --- end loading data

# STAligner analysis
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
# build graph information for each timepoint
Batch_list = []
adj_list = []
for idx, tp in enumerate(timepoints):
    tp_adata = adata[adata.obs['timepoint']==tp].copy()
    tp_adata.X = csr_matrix(tp_adata.X) # covert to csr_matrix, because STAligner.Cal_Spatial_Net only accept sparse matrix
    STAligner.Cal_Spatial_Net(tp_adata, rad_cutoff=1.3)
    # Normalization
    sc.pp.normalize_total(tp_adata, target_sum=1e4)
    sc.pp.log1p(tp_adata)
    sc.pp.highly_variable_genes(tp_adata, flavor="seurat_v3", n_top_genes=5000) #ensure enough common HVGs in the combined matrix
    tp_adata = tp_adata[:, tp_adata.var['highly_variable']]
    adj_list.append(tp_adata.uns['adj'])
    Batch_list.append(tp_adata)
# combine all timepoints
adata_concat = anndata.concat(Batch_list, label="slice_name", keys=timepoints)
adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')
print('adata_concat.shape: ', adata_concat.shape)
adj_concat = np.asarray(adj_list[0].todense())
for batch_id in range(1, len(timepoints)):
    adj_concat = scipy.linalg.block_diag(adj_concat, np.asarray(adj_list[batch_id].todense()))
adata_concat.uns['edgeList'] = np.nonzero(adj_concat)
# use last timepoint as reference
iter_comb = [(0, 4), (1, 4), (2, 4), (3, 4)] ## Fix slice 3 as reference to align
adata_concat = STAligner.train_STAligner(adata_concat, verbose=True, knn_neigh = 100, iter_comb = iter_comb, margin=2.5,  device=used_device)
adata_result = adata_concat.copy()
adata_result.uns.clear()
adata_result.write_h5ad(result_dir+os.sep+'STAligner_adata.h5ad')