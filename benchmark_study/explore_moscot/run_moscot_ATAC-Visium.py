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

base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_moscot/ATAC-Visium_noremoval_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
section_dict = {1: 't2', 2: 't1', 3: 't0'}
data_dict = {1: 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz', # E15.5
             2: 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz', # E13.5
             3: 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz' # E12.5
            }
adata = anndata.read_h5ad(data_dir+os.sep+'exprs.h5ad')
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
    size=[5]*len(library_ids),
    legend_fontsize=15,
)
plt.savefig(result_dir+os.sep+"spatial_scatter.png", bbox_inches="tight")

# --- run moscot
stp = SpatioTemporalProblem(adata)
stp = stp.score_genes_for_marginals(
    gene_set_proliferation="mouse", gene_set_apoptosis="mouse"
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
    # cell_transition_df = stp.cell_transition(
    #     source=tp_i,
    #     target=tp_j,
    #     source_groups='annotation', # does not matter what source_groups/target_groups are
    #     target_groups='annotation',
    #     aggregation_mode='cell', 
    #     forward=False,
    # )
    # cell_transition_df.to_csv(result_dir+os.sep+f'{tp_i}_{tp_j}_cell_transition.csv')
    # tp_j_adata.obs['transitioned_annotation'] = cell_transition_df.idxmax(axis=0).values
    # from sklearn.metrics import adjusted_rand_score
    # ari = adjusted_rand_score(tp_j_adata.obs['annotation'], tp_j_adata.obs['transitioned_annotation'])