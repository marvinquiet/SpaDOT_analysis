# conda activate /net/zootopia/disk1/wenjinma/envs/spateo
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '3'

import torch
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print("Running this notebook on: ", device)

import spateo as st
print("Last run with spateo version:", st.__version__)

# Other imports
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad

import numpy as np
import pandas as pd

base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_Artista"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spateo/Axolotl_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
#adata = anndata.read_h5ad(data_dir+os.sep+'Development_curated.h5ad')
adata = ad.read_h5ad(data_dir+os.sep+'Development_addCP.h5ad')
adata.X = adata.layers['counts'].astype(int)
adata.var['original_genes'] = adata.var_names.tolist()
adata.var_names = adata.var['Gene'].values
adata.obs['status'] = [_.split('_')[0] for _ in  adata.obs['Batch']]
adata.obs['timepoint'] = adata.obs['status'].replace({'Stage44': 't0',
    'Stage54': 't1', 'Stage57': 't2', 'Injury': 't3', 'Batch1': 't4', 'Meta': 't5'})
# focus on first three time points
timepoints = ['t0', 't1', 't2']
adata = adata[adata.obs['timepoint'].isin(timepoints)].copy()
adata.obs['pixel_x'] = adata.obsm['spatial'][:, 0]
adata.obs['pixel_y'] = -adata.obsm['spatial'][:, 1]
adata.obs_names = adata.obs_names+'_'+adata.obs['timepoint']
adata.obs['annotation'] = adata.obs['region_annotation']
# --- end loading data

slice1 = adata[adata.obs['timepoint'] == 't1']
slice2 = adata[adata.obs['timepoint'] == 't2']
# preprocess slice1
sc.pp.filter_genes(slice1, min_cells=3)
slice1.layers["counts"] = slice1.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(slice1)
# Logarithmize the data
sc.pp.log1p(slice1)
# annotates highly variable genes
sc.pp.highly_variable_genes(slice1, n_top_genes=2000)
# preprocess slice2
sc.pp.filter_genes(slice2, min_cells=3)
# Saving count data
slice2.layers["counts"] = slice2.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(slice2)
# Logarithmize the data
sc.pp.log1p(slice2)
# annotates highly variable genes
sc.pp.highly_variable_genes(slice2, n_top_genes=2000)

spatial_key = 'spatial'
cluster_key = 'annotation'
st.pl.slices_2d(
    slices = [slice1, slice2],
    label_key = cluster_key,
    spatial_key = spatial_key,
    height=4,
    center_coordinate=True,
    show_legend=True,
    legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, 0) ,'ncol': 5, 'borderaxespad': -4, 'frameon': False},
)
plt.savefig(result_dir+os.sep+'slices_1_vs_slices_2_spatial.png', bbox_inches='tight')

st.align.group_pca([slice1,slice2], pca_key='X_pca')
key_added = 'align_spatial'
# spateo return aligned slices as well as the mapping matrix
aligned_slices, pis = st.align.morpho_align(
    models=[slice1, slice2],
    ## Uncomment this if use highly variable genes
    # models=[slice1[:, slice1.var.highly_variable], slice2[:, slice2.var.highly_variable]],
    ## Uncomment the following if use pca embeddings
    # rep_layer='X_pca',
    # rep_field='obsm',
    # dissimilarity='cos',
    verbose=False,
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
)

st.pl.overlay_slices_2d(slices = aligned_slices, spatial_key = key_added, height=5, overlay_type='backward')
plt.savefig(result_dir+os.sep+'aligned_slices_1_vs_slices_2_spatial.png', bbox_inches='tight')


# --- perform nonrigid alignment
slice1 = adata[adata.obs['timepoint'] == 't0']
slice2 = adata[adata.obs['timepoint'] == 't1']
# preprocess slice1
sc.pp.filter_genes(slice1, min_cells=3)
slice1.layers["counts"] = slice1.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(slice1)
# Logarithmize the data
sc.pp.log1p(slice1)
# annotates highly variable genes
sc.pp.highly_variable_genes(slice1, n_top_genes=2000)
# preprocess slice2
sc.pp.filter_genes(slice2, min_cells=3)
# Saving count data
slice2.layers["counts"] = slice2.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(slice2)
# Logarithmize the data
sc.pp.log1p(slice2)
# annotates highly variable genes
sc.pp.highly_variable_genes(slice2, n_top_genes=2000)

spatial_key = 'spatial'
key_added = 'align_spatial'
st.align.group_pca([slice1, slice2], pca_key='X_pca')
aligned_slices, pis = st.align.morpho_align(
    models=[slice1, slice2],
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
    verbose=True,
    rep_layer='X_pca',
    rep_field='obsm',
    dissimilarity='cos',
    beta=1,
    lambdaVF=1,
    max_iter=500,
    K=100,
)
st.pl.overlay_slices_2d(slices = aligned_slices, spatial_key = key_added, height=3, overlay_type='backward', legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, 0),'borderaxespad': -2, 'frameon': False},)
plt.savefig(result_dir+os.sep+'rigid_aligned_slices_1_vs_slices_2_spatial.png', bbox_inches='tight')

st.pl.overlay_slices_2d(slices = aligned_slices, spatial_key = key_added+'_nonrigid', height=3, overlay_type='backward', legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, 0),'borderaxespad': -2, 'frameon': False},)
plt.savefig(result_dir+os.sep+'nonrigid_aligned_slices_1_vs_slices_2_spatial.png', bbox_inches='tight')

st.pl.overlay_slices_2d(slices = aligned_slices, spatial_key = key_added+'_nonrigid', label_key="annotation", height=3, overlay_type='backward', legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, 0),'borderaxespad': -2, 'frameon': False},)
plt.savefig(result_dir+os.sep+'nonrigid_aligned_slices_1_vs_slices_2_spatial_annotated.png', bbox_inches='tight')