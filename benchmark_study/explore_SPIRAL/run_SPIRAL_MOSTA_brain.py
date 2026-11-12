# conda activate /media/sda/wenjinma/envs/GraphST
import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import seaborn as sns

from utils import *
import scipy.io

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

## check with our simulation data
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_SPIRAL/MOSTA_brain_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'adata_brain_clean.h5ad')
n_samples = 20000  # Set desired number of cells
selected_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
adata = adata[selected_indices, :].copy()
adata.X = adata.layers['count'].astype(int).toarray()
adata.obs['annotation'] = 'Brain'
# --- end loading data

# preprocess data and train SPIRAL
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
Preprocess(adata, result_dir=result_dir) # default KNN=10, rad=10
train_SPIRAL(result_dir, timepoints) # default KNN=10
spiral_adata = sc.read_h5ad(result_dir+os.sep+"SPIRAL.h5ad")
spiral_adata.obs = adata.obs.loc[spiral_adata.obs_names, :]
spiral_adata.write_h5ad(result_dir+os.sep+"SPIRAL_corrected.h5ad")