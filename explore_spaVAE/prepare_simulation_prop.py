import os, argparse
import random
import h5py
import anndata
import numpy as np

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_simulation_brain"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/simulation_prop_results"
os.makedirs(result_dir, exist_ok=True)

simulation_data_dir = "/media/sda/wenjinma/projects/spatialATAC-explore/results/spatialPCA_simulation_brain_proportion"
# load simulated data
adata = anndata.read_h5ad(simulation_data_dir+os.sep+'adata.h5ad')
random.seed(1234)
sampled_idx = random.sample(range(adata.shape[0]), 3000*3)
adata = adata[sampled_idx].copy()
adata.obs.rename(columns={'x_coor': 'pixel_x', 'y_coor': 'pixel_y', 'label': 'annotation'}, inplace=True)
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
batch_info = np.array(adata.obs['timepoint_numeric']).astype('int')
batch_mat = np.zeros((batch_info.size, batch_info.max()+1))
batch_mat[np.arange(batch_info.size), batch_info] = 1
adata.obs.to_csv(result_dir+os.sep+'combined_simulation_metadata.csv')

with h5py.File(result_dir+os.sep+'combined_simulation.h5', 'w') as file:
    file.create_dataset('X', data=adata.X.A)
    file.create_dataset('pos', data=loc)
    file.create_dataset('batch', data=batch_mat)