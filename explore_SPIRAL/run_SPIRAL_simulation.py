import os, argparse
import sys
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import random
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import seaborn as sns
import time

from utils import *

## check with our simulation data
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_simulation_brain"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_SPIRAL/simulation_results"
os.makedirs(result_dir, exist_ok=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Run SPIRAL on simulation data")
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
    if os.path.exists(simulation_result_dir+os.sep+'SPIRAL_running_time.txt'):
        print('SPIRAL has been run on this simulation data')
        print(simulation_result_dir)
        sys.exit(0)

    adata = sc.read_h5ad(simulation_data_dir+os.sep+'adata.h5ad')
    random.seed(1234)
    sampled_idx = random.sample(range(adata.shape[0]), 3000*3)
    adata = adata[sampled_idx].copy()
    adata.obs.rename(columns={'x_coor': 'pixel_x', 'y_coor': 'pixel_y', 'label': 'annotation'}, inplace=True)
    loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
    adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
    adata.obsm['spatial'] = loc
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()
    # preprocess data and train SPIRAL
    start_time = time.time()
    Preprocess(adata, result_dir=simulation_result_dir) # default KNN=10, rad=10
    train_SPIRAL(simulation_result_dir, timepoints) # default KNN=10
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'SPIRAL_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')

    # spiral_adata = sc.read_h5ad(simulation_result_dir+os.sep+"SPIRAL.h5ad")
    # spiral_adata.obs = adata.obs.loc[spiral_adata.obs_names, :]
    # spiral_adata.write_h5ad(simulation_result_dir+os.sep+"SPIRAL_corrected.h5ad")
    # # --- explore kmeans clustering on the STAligner embedding
    # tp_kmeans_dict = dict()
    # ARI_dict = dict()
    # for tp in timepoints:
    #     tp_adata = spiral_adata[spiral_adata.obs['timepoint'] == tp]
    #     tp_emb = tp_adata.obsm['spiral']
    #     tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), random_state=0).fit(tp_emb)
    #     tp_kmeans_dict[tp] = tp_kmeans
    #     tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
    #     tp_ARI = adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])
    #     print('%s Kmeans ARI: %.3f' % (tp, tp_ARI))
    #     ARI_dict[tp] = tp_ARI
    #     plt.figure(figsize=(5,5))
    #     ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='kmeans', s=10)
    #     plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #     plt.tight_layout()
    #     plt.savefig(simulation_result_dir+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
    #     plt.close()
    # pd.DataFrame(ARI_dict, index=[0]).to_csv(simulation_result_dir+os.sep+'SPIRAL_ARI.csv', index=False)
    # sc.tl.umap(spiral_adata)
    # spiral_adata.obs['tp_annotation'] = spiral_adata.obs['timepoint'].astype(str) + '_' + spiral_adata.obs['annotation'].astype(str)
    # sc.pl.umap(spiral_adata, color=["timepoint", "celltype", "tp_annotation"])
    # plt.savefig(simulation_result_dir+os.sep+"SPIRAL_umap.png")
    # plt.close()