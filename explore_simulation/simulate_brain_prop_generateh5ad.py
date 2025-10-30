import os, argparse
import anndata
from scipy.io import mmread
import numpy as np
import pandas as pd

data_dir = '/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop'

def parse_args():
    parser = argparse.ArgumentParser(description="Generate simulation h5ad")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    parser.add_argument('--prop', type=float, help='Major cell type proportion')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    simulation_data_dir = data_dir+os.sep+str(args.seed)+os.sep+'prop_'+str(args.prop)
    # load data from three time points
    t0_mtx = mmread(simulation_data_dir+os.sep+'t0.mtx')
    gene_names = [f'Feature{i}' for i in range(1, t0_mtx.shape[0]+1)]
    t0_metadata = pd.read_csv(simulation_data_dir+os.sep+'t0_metadata.csv', index_col=0)
    t0_metadata['timepoint'] = 't0'
    t0_metadata.index = t0_metadata.index.astype(str) + '_t0'
    t0_adata = anndata.AnnData(X=t0_mtx.T, obs=t0_metadata, var=pd.DataFrame(index=gene_names))
    t1_mtx = mmread(simulation_data_dir+os.sep+'t1.mtx')
    t1_metadata = pd.read_csv(simulation_data_dir+os.sep+'t1_metadata.csv', index_col=0)
    t1_metadata['timepoint'] = 't1'
    t1_metadata.index = t1_metadata.index.astype(str) + '_t1'
    t1_adata = anndata.AnnData(X=t1_mtx.T, obs=t1_metadata, var=pd.DataFrame(index=gene_names))
    t2_mtx = mmread(simulation_data_dir+os.sep+'t2.mtx')
    t2_metadata = pd.read_csv(simulation_data_dir+os.sep+'t2_metadata.csv', index_col=0)
    t2_metadata['timepoint'] = 't2'
    t2_metadata.index = t2_metadata.index.astype(str) + '_t2'
    t2_adata = anndata.AnnData(X=t2_mtx.T, obs=t2_metadata, var=pd.DataFrame(index=gene_names))
    adata = anndata.concat([t0_adata, t1_adata, t2_adata])
    adata.X = adata.X.toarray()

    # --- downsample adata
    n_cells = 6000
    np.random.seed(args.seed)
    idx = np.random.choice(adata.n_obs, n_cells, replace=False)
    adata_sub = adata[idx, :]
    adata_sub.write_h5ad(simulation_data_dir+os.sep+'adata_downsampled.h5ad')