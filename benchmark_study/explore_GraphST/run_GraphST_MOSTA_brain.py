# conda activate /media/sda/wenjinma/envs/GraphST/
import os
import torch
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from GraphST import GraphST
from GraphST.utils import clustering

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

device = torch.device('cuda')
os.environ['R_HOME'] = '/media/sda/wenjinma/envs/GraphST/bin/R'
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_GraphST/MOSTA_brain_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'adata_brain_clean.h5ad')
n_samples = 20000  # Set desired number of cells
selected_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
adata = adata[selected_indices, :].copy()
adata.X = adata.layers['count'].astype(int).toarray()
adata.obs['annotation'] = 'Brain'
# --- end loading data

adata.obs['data'] = adata.obs['timepoint']
# plot algined image
plt.rcParams["figure.figsize"] = (3, 3)
ax = sc.pl.embedding(adata, basis='spatial',
                color='data',
                show=False)
ax.set_title('Aligned image')
plt.tight_layout()
plt.savefig(result_dir+os.sep+'aligned_image.png')
plt.close()

# train GraphST model
model = GraphST.GraphST(adata, device=device)
adata = model.train()

# plot images
fig, ax_list = plt.subplots(1, 3, figsize=(12, 3))
### Plotting UMAP before batch effect correction
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color='data', title='Uncorrected',
                ax = ax_list[0],
                show=False)
### Plotting UMAP after batch effect correction
sc.pp.neighbors(adata, use_rep='emb', n_neighbors=10)
sc.tl.umap(adata)
sc.pl.umap(adata,
        color='data',
        ax=ax_list[1],
        title='Batch corrected',
        #legend_loc = 'bottom margin',
        show=False)
### Color by predicted domains
sc.pl.umap(adata, color='annotation', ax=ax_list[2], title='Colored by region', show=False)
plt.tight_layout(w_pad=0.02)
plt.savefig(result_dir+os.sep+'umap.png')
plt.close()
# save adata for future analysis
adata.write(result_dir+os.sep+'GraphST_adata.h5ad')