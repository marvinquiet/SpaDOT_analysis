import os
import h5py
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

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
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/MOSTA_brain_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data (preprocessed in previous scripts)
adata = anndata.read_h5ad(data_dir+os.sep+'adata_brain_clean.h5ad')
n_samples = 20000  # Set desired number of cells
selected_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
adata = adata[selected_indices, :].copy()
adata.obs['annotation'] = 'Brain'
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
adata.obs.to_csv(result_dir+os.sep+'combined_metadata.csv')
print("Data size:", adata.shape)

loc = np.array(adata.obsm['spatial']).astype('float64')
batch_info = np.array(adata.obs['timepoint_numeric']).astype('int')
batch_mat = np.zeros((batch_info.size, batch_info.max()+1))
batch_mat[np.arange(batch_info.size), batch_info] = 1
 
with h5py.File(result_dir+os.sep+'combined.h5', 'w') as file:
    file.create_dataset('X', data=adata.X.toarray())
    file.create_dataset('pos', data=loc)
    file.create_dataset('batch', data=batch_mat)