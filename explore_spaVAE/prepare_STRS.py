import os
import h5py
import anndata
import numpy as np
import pandas as pd
import scipy

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_STRS"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/STRS_results"
os.makedirs(result_dir, exist_ok=True)

# load data
exprs = scipy.io.mmread(data_dir+os.sep+'STRS.mtx.gz').toarray()
genes = np.loadtxt(data_dir+os.sep+'STRS_genes.txt', dtype='str')
metadata = pd.read_csv(data_dir+os.sep+'STRS_metadata.csv', header=0, index_col=0)
adata = anndata.AnnData(X=exprs.T, obs=metadata, var=pd.DataFrame(index=genes, data={'gene_id':genes}))
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
batch_info = np.array(adata.obs['timepoint_numeric']).astype('int')
batch_mat = np.zeros((batch_info.size, batch_info.max()+1))
batch_mat[np.arange(batch_info.size), batch_info] = 1
adata.obs.to_csv(result_dir+os.sep+'combined_STRS_metadata.csv')

with h5py.File(result_dir+os.sep+'combined_STRS.h5', 'w') as file:
    file.create_dataset('X', data=adata.X)
    file.create_dataset('pos', data=loc)
    file.create_dataset('batch', data=batch_mat)