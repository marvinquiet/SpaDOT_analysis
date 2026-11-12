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
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_DeST-OT/ATAC-Visium_noremoval_results"
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
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes

adata.obs.to_csv(result_dir+os.sep+'combined_ATAC-Visium_metadata.csv')
print("Data size:", adata.shape)
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']

# write anndata to data_dir
# adata.obs['timepoint'] = adata.obs['timepoint'].astype(str)
# adata.write_h5ad('/net/mulan/home/wenjinma/projects/SpaDOT/examples/MouseOrganogenesis.h5ad')

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
