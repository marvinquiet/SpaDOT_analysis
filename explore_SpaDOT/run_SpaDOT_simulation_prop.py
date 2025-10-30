# conda activate ~/envs/SpaDOT
import os
from argparse import Namespace
import random
import argparse
import torch
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import time

# --- for neat output
import warnings
# Suppress all UserWarnings
warnings.filterwarnings("ignore", category=UserWarning)
# Suppress all FutureWarnings
warnings.filterwarnings("ignore", category=FutureWarning)

import SpaDOT

device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')
base_dir = "/media/sda/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop"
result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_SpaDOT/simulation_prop_results"
os.makedirs(result_dir, exist_ok=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Run SEDR on simulation data")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    parser.add_argument('--prop', type=float, help='Simulated proportion of major cell type in each region')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    simulation_data_dir = data_dir+os.sep+str(args.seed)+os.sep+'prop_'+str(args.prop)
    simulation_result_dir = result_dir+os.sep+str(args.seed)+os.sep+'prop_'+str(args.prop)
    os.makedirs(simulation_result_dir, exist_ok=True)
    adata = sc.read_h5ad(simulation_data_dir+os.sep+'adata_downsampled.h5ad')
    adata.obs.rename(columns={'x_coor': 'pixel_x', 'y_coor': 'pixel_y', 'label': 'annotation'}, inplace=True)
    loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
    adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
    adata.obsm['spatial'] = loc
    adata.obs['data'] = adata.obs['timepoint']
    adata.var['gene_names'] = adata.var_names
    if not sp.issparse(adata.X):
        adata.X = sp.csc_matrix(adata.X)
    adata.write_h5ad(simulation_result_dir+os.sep+'adata.h5ad')

    start_time = time.time()
    # --- create args
    preprocess_args = Namespace(
        data=simulation_result_dir+os.sep+'adata.h5ad',
        prefix='preprocessed_',
        feature_selection=True,
        output_dir=simulation_result_dir
    )
    SpaDOT.preprocess(preprocess_args)
    train_args = Namespace(
        data=simulation_result_dir+os.sep+'preprocessed_adata.h5ad',
        output_dir=simulation_result_dir,
        prefix="",
        config = None, # use default configuration
        save_model = True,
        device = 'cuda:2' # if GPU is not available, you can also use cpu, but not recommended
    )
    SpaDOT.train(train_args)
    # end timer
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'SpaDOT_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')

    # explore clustering results revised embeddings
    latent_adata = anndata.read_h5ad(simulation_result_dir+os.sep+'latent.h5ad')
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in set(latent_adata.obs['timepoint']):
        tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.X
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), random_state=0).fit(tp_emb)
        tp_kmeans_dict[tp] = tp_kmeans
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        tp_ARI = adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])
        ARI_dict[tp] = tp_ARI
        print('%s Kmeans ARI: %.3f' % (tp, tp_ARI))
        plt.figure(figsize=(5,5))
        ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='kmeans', s=10)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(simulation_result_dir+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
        plt.close()
    pd.DataFrame(ARI_dict, index=[0]).to_csv(simulation_result_dir+os.sep+f'ARI.csv', index=False)

