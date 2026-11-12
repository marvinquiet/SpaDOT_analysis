# conda activate /media/sda/wenjinma/envs/spatialATAC
import os, sys, argparse
os.environ["CUDA_VISIBLE_DEVICES"] = "3"
import math
import h5py
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import time

import matplotlib.pyplot as plt
from matplotlib import transforms
import seaborn as sns
import colorcet as cc
sns.set_palette(cc.glasbey)
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


sys.path.append('/net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_spaVAE/spaVAE/src/spaVAE_Batch/')

def parse_args():
    parser = argparse.ArgumentParser(description="Run spaVAE on simulation data")
    parser.add_argument('--seed', type=int, help='Simulation seed')
    parser.add_argument('--prop', type=float, help='Simulated proportion of major cell type in each region')
    parser.add_argument('--DEG', type=int, default=500, help='Number of DEGs per cell type')
    parser.add_argument('--mixedNPCs', action='store_true', help='Whether the NPCs is a mixture')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    if args.mixedNPCs:
        data_dir = f"/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_mixedNPCs_DEG{args.DEG}"
        result_dir = f"/net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_spaVAE/simulation_prop_results_mixedNPCs_DEG{args.DEG}"
    else:
        data_dir = f"/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_DEG{args.DEG}"
        result_dir = f"/net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_spaVAE/simulation_prop_results_DEG{args.DEG}"
    os.makedirs(result_dir, exist_ok=True)
    simulation_data_dir = data_dir+os.sep+str(args.seed)+os.sep+'prop_'+str(args.prop)
    simulation_result_dir = result_dir+os.sep+str(args.seed)+os.sep+'prop_'+str(args.prop)

    # run spaVAE Batch version
    start_time = time.time()
    os.system('python /net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_spaVAE/spaVAE/src/spaVAE_Batch/run_spaVAE_Batch.py --data_file {0}/combined_simulation.h5 --model_file {0}/model.pt --final_latent_file {0}/final_latent.txt --denoised_counts_file {0}/denoised_counts.txt > {0}/log.out 2>&1'.format(simulation_result_dir))  
    end_time = time.time()
    with open(simulation_result_dir+os.sep+'spaVAE_running_time.txt', 'w') as f:
        f.write(str(end_time-start_time))
        f.write('\n')

    # analyze results
    data_mat = h5py.File(simulation_result_dir+os.sep+'combined_simulation.h5', 'r')
    loc = np.array(data_mat['pos']).astype('float64')
    data_mat.close()
    
    metadata = pd.read_csv(simulation_result_dir+os.sep+'combined_simulation_metadata.csv', index_col=0)
    latent = np.loadtxt(simulation_result_dir+os.sep+'final_latent.txt', delimiter=',')
    latent_adata = anndata.AnnData(X=latent, obs=metadata)
    latent_adata.obs['annotation'] = latent_adata.obs['annotation'].astype('category')
    sc.pp.neighbors(latent_adata, n_neighbors=round(math.sqrt(latent_adata.shape[0])))
    sc.tl.umap(latent_adata)
    sc.pl.umap(latent_adata, color=['annotation', 'timepoint'],
            ncols=3)
    plt.savefig(simulation_result_dir+os.sep+'spaVAE_UMAP.png')
    plt.close()

    # --- calculate ARI per timepoint
    timepoints = list(set(latent_adata.obs['timepoint']))
    timepoints.sort()
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in timepoints:
        tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp]
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), random_state=0).fit(tp_adata.X)
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
    pd.DataFrame(ARI_dict, index=[0]).to_csv(simulation_result_dir+os.sep+'spaVAE_ARI.csv', index=False)
