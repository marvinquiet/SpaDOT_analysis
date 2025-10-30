# conda activate /media/sda/wenjinma/envs/GraphST/
# conda activate /net/zootopia/disk1/wenjinma/envs/GraphST/
import os
import random
import argparse
import torch
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import time

from GraphST import GraphST
from GraphST.utils import clustering

device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')
os.environ['R_HOME'] = '/usr/bin/R'
# base_dir = "/media/sda/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop"
result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_GraphST/simulation_prop_results"
os.makedirs(result_dir, exist_ok=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Run GraphST on simulation data")
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
    # plot algined image
    plt.rcParams["figure.figsize"] = (3, 3)
    ax = sc.pl.embedding(adata, basis='spatial',
                    color='data',
                    show=False)
    ax.set_title('Aligned image')
    plt.tight_layout()
    plt.savefig(simulation_result_dir+os.sep+'aligned_image.png')
    plt.close()
    # train GraphST model
    start_time = time.time()
    model = GraphST.GraphST(adata, device=device)
    adata = model.train()
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'GraphST_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')

    # # clustering(adata, n_clusters=10, method='louvain', start=0.1, end=2.0, increment=0.01)
    # # plot images
    # fig, ax_list = plt.subplots(1, 4, figsize=(15, 3))
    # ### Plotting UMAP before batch effect correction
    # sc.pp.normalize_total(adata)
    # sc.pp.log1p(adata)
    # sc.pp.pca(adata)
    # sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=10, n_pcs=40)
    # sc.tl.umap(adata)
    # sc.pl.umap(adata, color='data', title='Uncorrected',
    #                 ax = ax_list[0],
    #                 show=False)
    # ### Plotting UMAP after batch effect correction
    # sc.pp.neighbors(adata, use_rep='emb', n_neighbors=10)
    # sc.tl.umap(adata)
    # sc.pl.umap(adata,
    #         color='data',
    #         ax=ax_list[1],
    #         title='Batch corrected',
    #         #legend_loc = 'bottom margin',
    #         show=False)
    # ### Color by predicted domains
    # adata.obs['tp_annotation'] = adata.obs['timepoint'].astype(str)+'-'+adata.obs['annotation'].astype(str)
    # sc.pl.umap(adata, color='tp_annotation', ax=ax_list[2], title='Colored by region', show=False)
    # sc.pl.umap(adata, color='celltype', ax=ax_list[3], title='Colored by celltype', show=False)
    # plt.tight_layout(w_pad=0.02)
    # plt.savefig(simulation_result_dir+os.sep+'umap.png')
    # plt.close()
    # save adata for future analysis
    adata.write(simulation_result_dir+os.sep+'GraphST_adata.h5ad')
    # calculate clustering results per timepoint
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in set(adata.obs['timepoint']):
        tp_adata = adata[adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['emb']
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
    pd.DataFrame(ARI_dict, index=[0]).to_csv(simulation_result_dir+os.sep+'ARI.csv', index=False)
