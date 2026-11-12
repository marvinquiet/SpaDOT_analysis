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
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spateo/ATAC-Visium_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
section_dict = {1: 't2', 2: 't1', 3: 't0'}
data_dict = {1: 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz', # E15.5
             2: 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz', # E13.5
             3: 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz' # E12.5
            }
adata = ad.read_h5ad(data_dir+os.sep+'exprs.h5ad')
adata.obs['x'] = None
adata.obs['y'] = None
adata.obs['pixel_x'] = None
adata.obs['pixel_y'] = None
loc_metadata_list = []
for key, value in data_dict.items():
    metadata = pd.read_csv(data_dir+os.sep+value, sep=',', index_col=0, header=None)
    metadata.columns = ['tissue', 'y', 'x', 'pixel_y', 'pixel_x']
    metadata['barcodes'] = metadata.index
    metadata.index = metadata['barcodes']+'_'+str(key)
    key_barcodes = adata[adata.obs['section'] == key].obs_names.tolist()
    common_barcodes = list(set(metadata.index.tolist()).intersection(set(key_barcodes)))
    section_metadata = metadata.loc[common_barcodes, :]
    loc_metadata_list.append(section_metadata)
loc_metadata = pd.concat(loc_metadata_list, axis=0)
common_barcodes = list(set(loc_metadata.index.tolist()).intersection(set(adata.obs_names.tolist())))
adata = adata[common_barcodes].copy()
adata.obs.loc[common_barcodes, 'x'] = loc_metadata.loc[common_barcodes, 'x'].values
adata.obs.loc[common_barcodes, 'y'] = loc_metadata.loc[common_barcodes, 'y'].values
adata.obs.loc[common_barcodes, 'pixel_x'] = loc_metadata.loc[common_barcodes, 'pixel_x'].values
adata.obs.loc[common_barcodes, 'pixel_y'] = loc_metadata.loc[common_barcodes, 'pixel_y'].values
adata.obs['annotation'] = adata.obs['seurat_clusters_harmony']
adata.obs['annotation'].replace({3: 'CNS', 2: 'CNS', 4: 'CNS',
                                 1: 'Men/PNS', 9: 'Men/PNS', 11: 'Men/PNS',
                                 0: 'Mesenchyme', 5: 'Mesenchyme', 6: 'Mesenchyme', 8: 'Mesenchyme',
                                 7: 'Liver', 10: 'Heart'}, inplace=True)
adata.obs['annotation'] = adata.obs['annotation'].astype('category')
adata.obs['timepoint'] = [section_dict[_] for _ in adata.obs['section']]
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
# remove those cells that are not in the image
with open(data_dir+os.sep+'t2_removed_cells.txt', 'r') as f:
    t2_removed_cells = f.readlines()
t2_removed_cells = [_.strip() for _ in t2_removed_cells]
adata = adata[~adata.obs_names.isin(t2_removed_cells)]
adata.obs.to_csv(result_dir+os.sep+'combined_ATAC-Visium_metadata.csv')
print("Data size:", adata.shape)
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']
adata.layers['count'] = adata.X.astype('float32')
# --- end loading data

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