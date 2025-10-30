import os, argparse
import anndata
from scipy.io import mmread
import pandas as pd

base_dir = "/net/mulan/home/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_simulation_brain"

def parse_args():
    parser = argparse.ArgumentParser(description="Run spaGCN on simulation data")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    parser.add_argument('--DEG_x', type=int, help='DEG between region split')
    parser.add_argument('--DEG_y', type=int, help='DEG betweeen different regions')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    simulation_data_dir = data_dir+os.sep+str(args.seed)+os.sep+str(args.DEG_x)+'_'+str(args.DEG_y)
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
    adata.write_h5ad(simulation_data_dir+os.sep+'adata.h5ad')