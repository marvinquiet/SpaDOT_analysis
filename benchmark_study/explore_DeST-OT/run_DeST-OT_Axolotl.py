# conda activate ~/envs/spatialATAC
import os,sys
sys.path.append('/net/mulan/home/wenjinma/projects/spatialATAC-explore/explore_DeST-OT/DeST_OT')

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from sklearn.metrics import adjusted_rand_score

from src.destot.DESTOT import align, xi_to_growth_rate
from src.destot.metrics import growth_distortion_metric, migration_metric
from src.destot.plotting import plot_slice_value
from src.destot.plotting import get_transition

base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_Artista"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_DeST-OT/Axolotl_results"
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
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
# --- end loading data

# adata.write_h5ad("/net/mulan/home/wenjinma/projects/SpaDOT/examples/Axolotl.h5ad")
timepoints = list(set(adata.obs['timepoint_numeric']))
timepoints.sort()
for i in range(len(timepoints)-1):
    tp_i = timepoints[i]
    tp_j = timepoints[i+1]
    tp_i_adata = adata[adata.obs['timepoint_numeric'] == tp_i].copy()
    tp_j_adata = adata[adata.obs['timepoint_numeric'] == tp_j].copy()
    Pi, xi = align(tp_i_adata, tp_j_adata, alpha=0.2, gamma=50, epsilon=0.1, max_iter=100, 
                       balanced=False, use_gpu=False, normalize_xi=True, check_convergence=False)
    plot_slice_value(tp_i_adata, xi)
    plt.savefig(result_dir+os.sep+f'{tp_i}_growth_value.png')
    print(f'The cell migration distance between {tp_i} and {tp_j} is: \
        {growth_distortion_metric(tp_i_adata, tp_j_adata, Pi, annotation_key="annotation")}')
    # --- column normalized
    T, labelX, labelY = get_transition(Pi, tp_i_adata, tp_j_adata, 
                                        annotation_key = 'annotation', 
                                        plotting=True, column = True)
    plt.tight_layout()
    plt.savefig(result_dir+os.sep+f'{tp_i}_{tp_j}_column_normalized.png')
    transition_df = pd.DataFrame(T, index=labelX, columns=labelY)
    transition_df.to_csv(result_dir+os.sep+f'{tp_i}_{tp_j}_column_normalized.csv')
    # --- row normalized
    T, labelX, labelY = get_transition(Pi, tp_i_adata, tp_j_adata, 
                                        annotation_key = 'annotation', 
                                        plotting=True, column = False)
    plt.tight_layout()
    plt.savefig(result_dir+os.sep+f'{tp_i}_{tp_j}_row_normalized.png')
    # --- write transition table
    transition_df = pd.DataFrame(T, index=labelX, columns=labelY)
    transition_df.to_csv(result_dir+os.sep+f'{tp_i}_{tp_j}_row_normalized.csv')
