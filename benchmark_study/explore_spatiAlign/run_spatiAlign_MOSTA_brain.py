# conda activate /media/sda/wenjinma/envs/STAligner_py38
import warnings
warnings.filterwarnings("ignore")

from spatialign import Spatialign

# the location of R (used for the mclust clustering)
import os
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
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

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spatiAlign/MOSTA_brain_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'adata_brain_clean.h5ad')
n_samples = 20000  # Set desired number of cells
selected_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
adata = adata[selected_indices, :].copy()
adata.obs['annotation'] = 'Brain'
# --- end loading data

# spatiAlign analysis
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
# split adata into batches
datasets = []
for tp in timepoints:
    adata_tp = adata[adata.obs['timepoint'] == tp].copy()
    adata_tp_file = result_dir+os.sep+str(tp)+'.h5ad'
    adata_tp.write_h5ad(adata_tp_file)
    datasets.append(adata_tp_file)

model = Spatialign(
    *datasets,
    batch_key='batch',
    is_norm_log=True,
    is_scale=False,
    n_neigh=15,
    is_undirected=True,
    latent_dims=100,
    seed=42,
    gpu=0,
    save_path=result_dir+os.sep,
    is_verbose=True
) # use default tau1, tau2 and tau3
model.train()
model.alignment()
correct_adata_list = list()
for i in range(len(timepoints)):
    correct_adata = sc.read_h5ad(result_dir+os.sep+f"res/correct_data{i}.h5ad")
    correct_adata.var['genes'] = correct_adata.var_names.tolist()
    correct_adata_list.append(correct_adata)
merge_data = anndata.AnnData.concatenate(*correct_adata_list)
merge_data.write_h5ad(result_dir+os.sep+'spatiAlign_adata.h5ad')
sc.pp.neighbors(merge_data, use_rep="correct", random_state=42)
sc.tl.umap(merge_data, random_state=42)
sc.pl.umap(merge_data, color=["batch", "annotation"])
plt.savefig(result_dir+os.sep+'spatiAlign_UMAP.png')
plt.close()