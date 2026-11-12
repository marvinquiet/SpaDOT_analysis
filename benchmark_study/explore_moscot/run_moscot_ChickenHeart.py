# conda activate /net/zootopia/disk1/wenjinma/envs/moscot
import os
import warnings
import anndata

import moscot as mt
import moscot.plotting as mtp
from moscot import datasets
from moscot.problems.spatiotemporal import SpatioTemporalProblem

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import scanpy as sc
import squidpy as sq

warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", FutureWarning)

# --- example data
example_adata = datasets.mosta() # this used 2000 HVGs
example_adata.obs["time"] = example_adata.obs["time"].astype("category")

base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_ChickenHeart"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_moscot/ChickenHeart_results"
os.makedirs(result_dir, exist_ok=True)


# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'processed_ChickenHeart.h5ad') # curated region split
adata.obsm['spatial'] = adata.obs[['pixel_x', 'pixel_y']]
adata.obs['annotation'] = adata.obs['region']
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
# --- end loading data

adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3') # the example data contains 2,000 HVGs, so we do the same
adata = adata[:, adata.var['highly_variable']].copy()

library_ids = list(adata.obs['timepoint_numeric'].unique())
sq.pl.spatial_scatter(
    adata,
    color="annotation",
    frameon=False,
    shape=None,
    library_key="timepoint_numeric",
    library_id=library_ids,
    title=library_ids,
    figsize=(10, 10),
    size=[5, 5, 5, 5],
    legend_fontsize=15,
)
plt.savefig(result_dir+os.sep+"spatial_scatter.png", bbox_inches="tight")

# --- run moscot
stp = SpatioTemporalProblem(adata)
stp = stp.score_genes_for_marginals(
    gene_set_proliferation="human", gene_set_apoptosis="human"
) # chicken is not implemented, only mouse/human available

stp = stp.prepare(
    time_key="timepoint_numeric",
    spatial_key="spatial",
    joint_attr=None,
    # callback="local-pca", # would produce error
)

stp = stp.solve(alpha=0.5, epsilon=0, rank=200)

timepoints = list(set(adata.obs['timepoint_numeric']))
timepoints.sort()
for i in range(len(timepoints)-1):
    tp_i = timepoints[i]
    tp_j = timepoints[i+1]
    tp_i_adata = adata[adata.obs['timepoint_numeric'] == tp_i].copy()
    tp_j_adata = adata[adata.obs['timepoint_numeric'] == tp_j].copy()
    celltypes_source_1 = tp_i_adata.obs["annotation"].cat.categories
    celltypes_target_1 = tp_j_adata.obs["annotation"].cat.categories
    celltypes_intersection_1 = list(
        set(celltypes_source_1).intersection(celltypes_target_1)
    )
    source_ordered_1 = celltypes_intersection_1 + list(
        set(celltypes_source_1) - set(celltypes_target_1)
    )
    target_ordered_1 = celltypes_intersection_1 + list(
        set(celltypes_target_1) - set(celltypes_source_1)
    )
    moscot_transition_df = stp.cell_transition(
        source=tp_i,
        target=tp_j,
        source_groups={"annotation": source_ordered_1},
        target_groups={"annotation": target_ordered_1},
        forward=True,
    )
    mtp.cell_transition(stp, fontsize=12, figsize=(7, 7))
    plt.savefig(result_dir+os.sep+f'{tp_i}_{tp_j}_transition_row_normalized.png', bbox_inches="tight")
    moscot_transition_df.to_csv(result_dir+os.sep+f'{tp_i}_{tp_j}_row_normalized.csv')
    moscot_transition_df = stp.cell_transition(
        source=tp_i,
        target=tp_j,
        source_groups={"annotation": source_ordered_1},
        target_groups={"annotation": target_ordered_1},
        forward=False
    ) 
    moscot_transition_df.to_csv(result_dir+os.sep+f'{tp_i}_{tp_j}_column_normalized.csv')