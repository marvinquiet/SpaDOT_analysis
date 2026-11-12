# conda activate /media/sda/wenjinma/envs/SpaGCN_py39/
import warnings
warnings.filterwarnings("ignore")

import SpaGCN as spg

import os
from anndata import AnnData
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

import anndata

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

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaGCN/MOSTA_brain_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'adata_brain_clean.h5ad')
n_samples = 20000  # Set desired number of cells
selected_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
adata = adata[selected_indices, :].copy()
adata.X = adata.layers['count'].astype(int).toarray()
adata.obs['annotation'] = 'Brain'
# --- end loading data

# split adata into batches
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
adata_list, l_list, adj_list = [], [], []
for tp in timepoints:
    adata_tp = adata[adata.obs['timepoint']==tp].copy()
    pixel_x = adata_tp.obsm['spatial'][:, 0].tolist()
    pixel_y = adata_tp.obsm['spatial'][:, 1].tolist()
    s=1
    b=49
    adj_tp=spg.calculate_adj_matrix(x=pixel_x,y=pixel_y, beta=b, alpha=s, histology=False)
    adata_tp.raw = adata_tp
    sc.pp.normalize_total(adata_tp)
    sc.pp.log1p(adata_tp)
    p = 0.5
    l_tp = spg.search_l(p, adj_tp, start=0.01, end=1000, tol=0.01, max_run=100)
    adj_list.append(adj_tp)
    adata_list.append(adata_tp)
    l_list.append(l_tp)
# run spaGCN multi-section
res=0.6
#Set seed
r_seed=t_seed=n_seed=100
#Multi-section SpaGCN
clf=spg.multiSpaGCN()
clf.train(adata_list,adj_list,l_list,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
# explore spaGCN embeddings
spaGCN_embed = clf.embed
merge_data = AnnData.concatenate(*adata_list)
merge_data.obsm['spaGCN_embed'] = spaGCN_embed
sc.pp.neighbors(merge_data, use_rep="spaGCN_embed", random_state=42)
sc.tl.umap(merge_data, random_state=42)
sc.pl.umap(merge_data, color=["batch", "annotation"])
plt.savefig(result_dir+os.sep+'spaGCN_UMAP.png')
plt.close()
merge_data.write_h5ad(result_dir+os.sep+'spaGCN.h5ad')