# conda activate /media/sda/wenjinma/envs/STAligner_py38
# conda activate /net/zootopia/disk1/wenjinma/envs/spatialign_py38
import warnings
warnings.filterwarnings("ignore")

from spatialign import Spatialign

# the location of R (used for the mclust clustering)
import os, argparse
import anndata
import torch
import scanpy as sc
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import time


def parse_args():
    parser = argparse.ArgumentParser(description="Run spatiAlign on simulation data")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    parser.add_argument('--prop', type=float, help='Simulated proportion of major cell type in each region')
    parser.add_argument('--DEG', type=int, help='Number of DEGs per cell type')
    parser.add_argument('--mixedNPCs', action='store_true', help='Whether the NPCs is a mixture')
    parser.add_argument('--regionNPCs', action='store_true', help='Whether the NPCs is region specific')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    if args.mixedNPCs:
        data_dir = f"/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_mixedNPCs_DEG{args.DEG}"
        result_dir = f"/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_spatiAlign/simulation_prop_results_mixedNPCs_DEG{args.DEG}"
    elif args.regionNPCs:
        data_dir = f"/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_regionNPCs_DEG{args.DEG}"
        result_dir = f"/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_spatiAlign/simulation_prop_results_regionNPCs_DEG{args.DEG}"
    else:
        data_dir = f"/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_DEG{args.DEG}"
        result_dir = f"/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_spatiAlign/simulation_prop_results_DEG{args.DEG}"
    os.makedirs(result_dir, exist_ok=True)
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
    datasets = []
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()
    for tp in timepoints:
        adata_tp = adata[adata.obs['timepoint'] == tp].copy()
        adata_tp_file = simulation_result_dir+os.sep+str(tp)+'.h5ad'
        adata_tp.write_h5ad(adata_tp_file)
        datasets.append(adata_tp_file)

    start_time = time.time()
    model = Spatialign(
        *datasets,
        min_genes=0,
        min_cells=0,
        batch_key='batch',
        is_norm_log=True,
        is_scale=False,
        n_neigh=15,
        is_undirected=True,
        latent_dims=100,
        seed=42,
        gpu=1 if torch.cuda.is_available() else None,
        save_path=simulation_result_dir+os.sep,
        is_verbose=True,
        tau1=0.05, tau2=0.01, tau3=0.1) # use recommended tau1, tau2 and tau3
    model.train()
    model.alignment()
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'spatiAlign_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')
    correct_adata_list = list()
    for i in range(len(timepoints)):
        correct_adata = sc.read_h5ad(simulation_result_dir+os.sep+f"res/correct_data{i}.h5ad")
        correct_adata.var['genes'] = correct_adata.var_names.tolist()
        correct_adata_list.append(correct_adata)
    merge_data = anndata.AnnData.concatenate(*correct_adata_list)
    sc.pp.neighbors(merge_data, use_rep="correct", random_state=42)
    sc.tl.umap(merge_data, random_state=42)
    merge_data.obs['tp_annotation'] = merge_data.obs['timepoint'].astype(str) + '_' + merge_data.obs['annotation'].astype(str)
    sc.pl.umap(merge_data, color=["batch", "celltype", "tp_annotation"])
    plt.savefig(simulation_result_dir+os.sep+'spatiAlign_UMAP.png')
    plt.close()
    # --- calculate ARI per timepoint
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in timepoints:
        tp_adata = merge_data[merge_data.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['correct']
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
