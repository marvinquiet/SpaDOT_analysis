# conda activate /media/sda/wenjinma/envs/spatialATAC
import os, argparse
import random
import h5py
import anndata
import numpy as np

data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop"
result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_spaVAE/simulation_prop_results"
os.makedirs(result_dir, exist_ok=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Run spaVAE on simulation data")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    simulation_data_dir = data_dir+os.sep+str(args.seed)
    simulation_result_dir = result_dir+os.sep+str(args.seed)
    os.makedirs(simulation_result_dir, exist_ok=True)
    # load simulated data
    adata = anndata.read_h5ad(simulation_data_dir+os.sep+'adata_downsampled.h5ad')
    adata.obs.rename(columns={'x_coor': 'pixel_x', 'y_coor': 'pixel_y', 'label': 'annotation'}, inplace=True)
    loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
    adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
    batch_info = np.array(adata.obs['timepoint_numeric']).astype('int')
    batch_mat = np.zeros((batch_info.size, batch_info.max()+1))
    batch_mat[np.arange(batch_info.size), batch_info] = 1
    adata.obs.to_csv(simulation_result_dir+os.sep+'combined_simulation_metadata.csv')
    
    with h5py.File(simulation_result_dir+os.sep+'combined_simulation.h5', 'w') as file:
        file.create_dataset('X', data=adata.X)
        file.create_dataset('pos', data=loc)
        file.create_dataset('batch', data=batch_mat)