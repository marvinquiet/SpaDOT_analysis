# conda activate /net/zootopia/disk1/wenjinma/envs/moscot
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "3" 

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

base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_Artista"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_moscot/Axolotl_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
#adata = anndata.read_h5ad(data_dir+os.sep+'Development_curated.h5ad')
adata = anndata.read_h5ad(data_dir+os.sep+'Development_addCP.h5ad')
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

adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
# --- obtain converted genes
convert_gene_df = pd.read_csv(data_dir+os.sep+'Axolotl_human_orthologs.csv', sep=',', header=0, index_col=0)
adata.var = adata.var.merge(convert_gene_df, left_on='original_genes', right_index=True, suffixes=('', '_human'))
# remove those NaN in hs_genes
adata_filtered = adata[:, ~adata.var['hs_genes'].isna()].copy()
# keep first gene if duplicate
adata_filtered = adata_filtered[:, ~adata_filtered.var['hs_genes'].duplicated(keep='first')].copy()
sc.pp.highly_variable_genes(adata_filtered, n_top_genes=2000, flavor='seurat_v3') # the example data contains 2,000 HVGs, so we do the same
adata_filtered = adata_filtered[:, adata_filtered.var['highly_variable']].copy()
adata_filtered.var_names = adata_filtered.var['hs_genes'].values

library_ids = list(adata_filtered.obs['timepoint_numeric'].unique())
sq.pl.spatial_scatter(
    adata_filtered,
    color="annotation",
    frameon=False,
    shape=None,
    library_key="timepoint_numeric",
    library_id=library_ids,
    title=library_ids,
    figsize=(10, 10),
    size=[5]*len(library_ids),
    legend_fontsize=15,
)
plt.savefig(result_dir+os.sep+"spatial_scatter.png", bbox_inches="tight")

# --- run moscot
stp = SpatioTemporalProblem(adata_filtered)
stp = stp.score_genes_for_marginals(
    gene_set_proliferation="human", gene_set_apoptosis="human"
) # we use Axolotl human orthologs

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
