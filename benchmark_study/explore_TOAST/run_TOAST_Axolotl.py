# conda activate /net/zootopia/disk1/wenjinma/envs/TOAST
import os,sys
import warnings
warnings.filterwarnings('ignore')
sys.path.append('TOAST/')
from TOAST_utils import *

import anndata

base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_Artista"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_TOAST/Axolotl_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
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

# --- from TOAST
n_comps = 30  
k = 5 
alpha = 0.5  
epsilon = 0.1  

# adata.write_h5ad("/net/mulan/home/wenjinma/projects/SpaDOT/examples/Axolotl.h5ad")
timepoints = list(set(adata.obs['timepoint_numeric']))
timepoints.sort()
for i in range(len(timepoints)-1):
    tp_i = timepoints[i]
    tp_j = timepoints[i+1]
    tp_i_adata = adata[adata.obs['timepoint_numeric'] == tp_i].copy()
    tp_j_adata = adata[adata.obs['timepoint_numeric'] == tp_j].copy()
    fgw_sn = compute_spatial_fgw_alignment(tp_i_adata, tp_j_adata, alpha, epsilon, n_comps, k, annotation_col='annotation')
    # --- aggregate by annotation
    # make DataFrame with proper row/col labels
    transition_df = pd.DataFrame(fgw_sn,
                      index=tp_i_adata.obs_names,
                      columns=tp_j_adata.obs_names)
    # series of cluster labels aligned to rows and columns
    row_clusters = tp_i_adata.obs['annotation']
    col_clusters = tp_j_adata.obs['annotation']
    # sum rows by row-cluster, then sum columns by col-cluster => cluster x cluster matrix
    cluster_transition = transition_df.groupby(row_clusters).sum().groupby(by=col_clusters, axis=1).sum()
    transition_prob_col_norm_df = cluster_transition.div(cluster_transition.sum(axis=0), axis=1)
    transition_prob_col_norm_df.to_csv(result_dir+os.sep+f'{tp_i}_{tp_j}_column_normalized.csv')
    transition_prob_row_norm_df = cluster_transition.div(cluster_transition.sum(axis=1), axis=0)
    transition_prob_row_norm_df.to_csv(result_dir+os.sep+f'{tp_i}_{tp_j}_row_normalized.csv')