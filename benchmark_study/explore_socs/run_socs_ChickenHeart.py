# conda activate /net/zootopia/disk1/wenjinma/envs/socs
import os,sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as io
from scipy import stats
import scanpy as sc
import anndata as ad
import squidpy as sq
import copy
import seaborn as sns
import pickle
import math
import sklearn
import torch

sys.path.append('SOCS')
import socs
from socs.utils import row_normalize
from socs.analysis import map_vector_sampled,hillNumber_1,vec2vec

base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_ChickenHeart"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_socs/ChickenHeart_results"
os.makedirs(result_dir, exist_ok=True)


# --- load data
adata = ad.read_h5ad(data_dir+os.sep+'processed_ChickenHeart.h5ad') # curated region split
adata.obsm['spatial'] = adata.obs[['pixel_x', 'pixel_y']]
adata.obs['annotation'] = adata.obs['region']
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
# --- end loading data


timepoints = list(set(adata.obs['timepoint_numeric']))
timepoints.sort()
for i in range(len(timepoints)-1):
    tp_i = timepoints[i]
    tp_j = timepoints[i+1]
    tp_i_adata = adata[adata.obs['timepoint_numeric'] == tp_i].copy()
    tp_j_adata = adata[adata.obs['timepoint_numeric'] == tp_j].copy()
    tp_i_adata.layers["counts"] = tp_i_adata.X.copy()
    sc.pp.normalize_total(tp_i_adata, inplace=True,target_sum=63)
    sc.pp.log1p(tp_i_adata)
    tp_j_adata.layers["counts"] = tp_j_adata.X.copy()
    sc.pp.normalize_total(tp_j_adata, inplace=True,target_sum=65)
    sc.pp.log1p(tp_j_adata)
    adata_concat= ad.concat([tp_i_adata,tp_j_adata])
    sc.pp.pca(adata_concat, random_state= 0,n_comps=30)
    alpha = 0.1
    eps = 1e-3
    rho1 = 5000.0
    rho2 = 5000.0
    # --- we do not have structure key
    socs_model = socs.ot.SOCSModel(adata_concat,'timepoint_numeric',expr_key='X_pca',struct_key=None,
            alpha=alpha,eps=eps,rho=rho1,rho2=rho2)
    T_socs = socs_model.infer_map(0,1)
    T_socs_n = row_normalize(T_socs)
