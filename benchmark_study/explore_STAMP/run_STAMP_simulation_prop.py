# conda activate /net/zootopia/disk1/wenjinma/envs/stamp
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import numpy as np
import sctm
import squidpy as sq
import pandas as pd
import anndata

# the location of R (used for the mclust clustering)
import os
import argparse
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import time



def parse_args():
    parser = argparse.ArgumentParser(description="Run STAMP on simulation data")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    parser.add_argument('--prop', type=float, help='Simulated proportion of major cell type in each region')
    parser.add_argument('--DEG', type=int, default=500, help='Number of DEGs per cell type')
    parser.add_argument('--mixedNPCs', action='store_true', help='Whether the NPCs is a mixture')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    if args.mixedNPCs:
        data_dir = f"/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_mixedNPCs_DEG{args.DEG}"
        result_dir = f"/net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_STAMP/simulation_prop_results_mixedNPCs_DEG{args.DEG}"
    else:
        data_dir = f"/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_DEG{args.DEG}"
        result_dir = f"/net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_STAMP/simulation_prop_results_DEG{args.DEG}"
    os.makedirs(result_dir, exist_ok=True)
    simulation_data_dir = data_dir+os.sep+str(args.seed)+os.sep+'prop_'+str(args.prop)
    simulation_result_dir = result_dir+os.sep+str(args.seed)+os.sep+'prop_'+str(args.prop)
    os.makedirs(simulation_result_dir, exist_ok=True)
    adata = sc.read_h5ad(simulation_data_dir+os.sep+'adata_downsampled.h5ad')
    adata.obs.rename(columns={'x_coor': 'pixel_x', 'y_coor': 'pixel_y', 'label': 'annotation'}, inplace=True)
    loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
    adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
    adata.obs['timepoint_numeric'] = adata.obs['timepoint_numeric'].astype('float32')
    adata.obsm['spatial'] = loc
    adata.obs['data'] = adata.obs['timepoint']
    adata.layers['count'] = adata.X.astype('float32')
    
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()

    # --- run STAMP
    start_time = time.time()
    sq.gr.spatial_neighbors(adata, n_neighs=round(1 / 1000 * adata.n_obs))
    sctm.seed.seed_everything(0)
    print(adata)
    sctm.pp.filter_genes(adata, 0.05)
    print(adata)

    model = sctm.stamp.STAMP(
        adata,
        n_topics=20,
        layer="count",
        time_covariate_keys="timepoint_numeric",
        # gene_likelihood="poisson",
    )
    # We used a weighted sampler here as that the first timepoint is extremly small compared to the last timepoint.
    model.train(device=device)
    # end timer
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'STAMP_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')

    topic_prop = model.get_cell_by_topic()
    adata.obsm['X_stamp'] = topic_prop.values
    for i in topic_prop.columns:
        adata.obs[i] = topic_prop[i]
    adata.obs.to_csv(simulation_result_dir+os.sep+'STAMP_topic_proportions.csv')

    for tp in timepoints:
        tp_adata = adata[adata.obs['timepoint'] == tp]
        for i in topic_prop.columns:
            sctm.pl.spatial(tp_adata, color=i, size=5, vmax='p99')
            plt.savefig(simulation_result_dir+os.sep+f'{i}_{tp}_spatial.png')

    # --- use topic proportions to perform KMeans clustering
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in set(adata.obs['timepoint']):
        tp_adata = adata[adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['X_stamp']
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
    pd.DataFrame(ARI_dict, index=[0]).to_csv(simulation_result_dir+os.sep+'STAMP_ARI.csv', index=False)


