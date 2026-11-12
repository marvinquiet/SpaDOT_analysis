import os
import h5py
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_ChickenHeart"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/ChickenHeart_results"
os.makedirs(result_dir, exist_ok=True)


# --- load data (preprocessed in previous scripts)
adata = sc.read_h5ad(data_dir+os.sep+'processed_ChickenHeart.h5ad')
adata.obsm['spatial'] = adata.obs[['pixel_x', 'pixel_y']]
adata.obs['annotation'] = adata.obs['region']
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
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