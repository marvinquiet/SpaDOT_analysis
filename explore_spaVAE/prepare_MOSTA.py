import os
import h5py
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/MOSTA_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'adata_clean.h5ad')
adata.X = adata.layers['count'].astype(int)
adata.obs['timepoint'] = adata.obs['time'].replace({9.5: 't0',
    10.5: 't1', 11.5: 't2', 12.5: 't3', 13.5: 't4', 14.5: 't5', 15.5: 't6', 16.5: 't7'})
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
# count each time points
# t7    121724
# t6    113218
# t5    102492
# t4     77126
# t3     51315
# t2     30051
# t1      8494
# t0      5892
# downsample data
adata = sc.pp.subsample(adata, n_obs = 100000, copy = True)
# count each time points after subsampling
# t7    23834
# t6    22289
# t5    20049
# t4    15102
# t3    10003
# t2     5931
# t1     1658
# t0     1134
# focus on first three time points
adata.obs['pixel_x'] = adata.obsm['spatial'][:, 0]
adata.obs['pixel_y'] = adata.obsm['spatial'][:, 1]
adata.obs_names = adata.obs_names+'_'+adata.obs['timepoint']
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