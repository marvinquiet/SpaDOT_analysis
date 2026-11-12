import os
import h5py
import anndata
import numpy as np
import pandas as pd
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_Artista"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/Axolotl_results"
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
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
adata.obs['pixel_x'] = adata.obsm['spatial'][:, 0]
adata.obs['pixel_y'] = -adata.obsm['spatial'][:, 1]
adata.obs_names = adata.obs_names+'_'+adata.obs['timepoint']
adata.obs['annotation'] = adata.obs['region_annotation']
adata.obs.to_csv(result_dir+os.sep+'combined_metadata.csv')
print("Data size:", adata.shape)

loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
batch_info = np.array(adata.obs['timepoint_numeric']).astype('int')
batch_mat = np.zeros((batch_info.size, batch_info.max()+1))
batch_mat[np.arange(batch_info.size), batch_info] = 1
 
with h5py.File(result_dir+os.sep+'combined.h5', 'w') as file:
    file.create_dataset('X', data=adata.X.toarray())
    file.create_dataset('pos', data=loc)
    file.create_dataset('batch', data=batch_mat)
