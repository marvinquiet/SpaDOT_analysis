# conda activate /media/sda/wenjinma/envs/SEDR/
import os
import random
import argparse
import torch
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import time

import SEDR
import harmonypy as hm

device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain"
result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_SEDR/simulation_results"
os.makedirs(result_dir, exist_ok=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Run GraphST on simulation data")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    parser.add_argument('--DEG_x', type=int, help='DEG between region split')
    parser.add_argument('--DEG_y', type=int, help='DEG betweeen different regions')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    simulation_data_dir = data_dir+os.sep+str(args.seed)+os.sep+str(args.DEG_x)+'_'+str(args.DEG_y)
    simulation_result_dir = result_dir+os.sep+str(args.seed)+os.sep+str(args.DEG_x)+'_'+str(args.DEG_y)
    os.makedirs(simulation_result_dir, exist_ok=True)
    # adata = sc.read_h5ad(simulation_data_dir+os.sep+'adata.h5ad')
    # random.seed(1234)
    # sampled_idx = random.sample(range(adata.shape[0]), 3000*3)
    # adata = adata[sampled_idx].copy()
    adata = sc.read_h5ad(simulation_data_dir+os.sep+'adata_downsampled.h5ad')
    adata.obs.rename(columns={'x_coor': 'pixel_x', 'y_coor': 'pixel_y', 'label': 'annotation'}, inplace=True)
    loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
    adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
    adata.obsm['spatial'] = loc
    adata.obs['data'] = adata.obs['timepoint']

    # --- SEDR
    start_time = time.time()
    for idx, tp in enumerate(list(set(adata.obs['timepoint']))):
        tp_adata = adata[adata.obs['timepoint']==tp].copy()
        tp_adata.var_names_make_unique()
        tp_adata.obs['batch_name'] = tp
        graph_dict_tmp = SEDR.graph_construction(tp_adata, 12)
        if idx == 0:
            SEDR_adata = tp_adata
            graph_dict = graph_dict_tmp
            SEDR_adata.obs['proj_name'] = tp
            tp_adata.obs['proj_name'] = tp_adata.obs['proj_name'].astype(str)
        else:
            var_names = SEDR_adata.var_names.intersection(tp_adata.var_names)
            SEDR_adata = SEDR_adata[:, var_names]
            tp_adata = tp_adata[:, var_names]
            tp_adata.obs['proj_name'] = tp
            tp_adata.obs['proj_name'] = tp_adata.obs['proj_name'].astype(str)
            SEDR_adata = anndata.concat((SEDR_adata, tp_adata))
            graph_dict = SEDR.combine_graph_dict(graph_dict, graph_dict_tmp)
    SEDR_adata.layers['count'] = SEDR_adata.X
    sc.pp.normalize_total(SEDR_adata, target_sum=1e6)
    sc.pp.highly_variable_genes(SEDR_adata, flavor="seurat_v3", layer='count', n_top_genes=2000)
    SEDR_adata = SEDR_adata[:, SEDR_adata.var['highly_variable'] == True]
    sc.pp.scale(SEDR_adata)
    # Run PCA on embeddings
    from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable.
    adata_X = PCA(n_components=200, random_state=42).fit_transform(SEDR_adata.X)
    SEDR_adata.obsm['X_pca'] = adata_X
    sedr_net = SEDR.Sedr(SEDR_adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
    using_dec = False
    if using_dec:
        sedr_net.train_with_dec()
    else:
        sedr_net.train_without_dec()
    sedr_feat, _, _, _ = sedr_net.process()
    SEDR_adata.obsm['SEDR'] = sedr_feat
    # Use harmony to integrate embeddings from different batches and calculate revised embeddings
    meta_data = SEDR_adata.obs[['batch_name']]
    data_mat = SEDR_adata.obsm['SEDR']
    vars_use = ['batch_name']
    ho = hm.run_harmony(data_mat, meta_data, vars_use)
    res = pd.DataFrame(ho.Z_corr).T
    res_df = pd.DataFrame(data=res.values, columns=['X{}'.format(i+1) for i in range(res.shape[1])], index=SEDR_adata.obs.index)
    SEDR_adata.obsm[f'SEDR.Harmony'] = res_df
    # end timer
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'SEDR_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')

    # # plot embeddings from SEDR
    sc.pp.neighbors(SEDR_adata, use_rep='SEDR.Harmony', metric='cosine') # In the tutorial, cosine is used
    sc.tl.umap(SEDR_adata)
    SEDR_adata.obs['tp_annotation'] = SEDR_adata.obs['timepoint'].astype(str) + '_' + SEDR_adata.obs['annotation'].astype(str)
    sc.pl.umap(SEDR_adata, color=['tp_annotation', 'celltype', 'batch_name'], show=False)
    plt.savefig(simulation_result_dir+os.sep+'SEDR_umap.png')
    plt.close()
    SEDR_adata.write(simulation_result_dir+os.sep+'SEDR_adata.h5ad')
    # explore clustering results revised embeddings
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in set(SEDR_adata.obs['timepoint']):
        tp_adata = SEDR_adata[SEDR_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['SEDR.Harmony']
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


