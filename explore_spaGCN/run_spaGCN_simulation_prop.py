# conda activate /media/sda/wenjinma/envs/SpaGCN_py39/
import warnings
warnings.filterwarnings("ignore")

import SpaGCN as spg

# the location of R (used for the mclust clustering)
import os, argparse
from anndata import AnnData
import scanpy as sc
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import time

import torch

data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop"
result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_spaGCN/simulation_prop_results"
os.makedirs(result_dir, exist_ok=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Run spaGCN on simulation data")
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
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()
    # split adata into batches
    start_time = time.time()
    adata_list, l_list, adj_list = [], [], []
    for tp in timepoints:
        adata_tp = adata[adata.obs['timepoint']==tp].copy()
        pixel_x = adata_tp.obs['pixel_x'].tolist()
        pixel_y = adata_tp.obs['pixel_y'].tolist()
        s=1
        b=49
        adj_tp=spg.calculate_adj_matrix(x=pixel_x,y=pixel_y, beta=b, alpha=s, histology=False)
        # np.savetxt(simulation_result_dir+os.sep+str(tp)+'_adj.csv', adj_tp, delimiter=',')
        adata_tp.raw = adata_tp
        sc.pp.normalize_total(adata_tp)
        sc.pp.log1p(adata_tp)
        p = 0.5
        l_tp = spg.search_l(p, adj_tp, start=0.01, end=1000, tol=0.01, max_run=100)
        adj_list.append(adj_tp)
        adata_list.append(adata_tp)
        l_list.append(l_tp)
    # run spaGCN multi-section
    res=0.6
    #Set seed
    r_seed=t_seed=n_seed=100
    #Multi-section SpaGCN
    clf=spg.multiSpaGCN()
    clf.train(adata_list,adj_list,l_list,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'spaGCN_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')
    
    # explore spaGCN embeddings
    spaGCN_embed = clf.embed
    merge_data = AnnData.concatenate(*adata_list)
    merge_data.obsm['spaGCN_embed'] = spaGCN_embed
    sc.pp.neighbors(merge_data, use_rep="spaGCN_embed", random_state=42)
    sc.tl.umap(merge_data, random_state=42)
    merge_data.obs['tp_annotation'] = merge_data.obs['timepoint'].astype(str) + '_' + merge_data.obs['annotation'].astype(str)
    sc.pl.umap(merge_data, color=["batch", "celltype", "tp_annotation"])
    plt.savefig(simulation_result_dir+os.sep+'spaGCN_UMAP.png')
    plt.close()
    merge_data.write_h5ad(simulation_result_dir+os.sep+'spaGCN.h5ad')
    # --- calculate ARI per timepoint
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in timepoints:
        tp_adata = merge_data[merge_data.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['spaGCN_embed']
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), random_state=0).fit(tp_emb)
        tp_kmeans_dict[tp] = tp_kmeans
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        tp_ARI = adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])
        print('%s Kmeans ARI: %.3f' % (tp, tp_ARI))
        ARI_dict[tp] = tp_ARI
        plt.figure(figsize=(5,5))
        ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='kmeans', s=10)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(simulation_result_dir+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
        plt.close()
    pd.DataFrame(ARI_dict, index=[0]).to_csv(simulation_result_dir+os.sep+'ARI.csv', index=False)
