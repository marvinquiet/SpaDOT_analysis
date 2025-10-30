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
result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies/test_SpaDOT_package"
os.makedirs(result_dir, exist_ok=True)

if __name__ == "__main__":

    import argparse
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    parser = argparse.ArgumentParser(description="SpaDOT ablation studies.")
    parser.add_argument('--data', dest='data', type=str)
    parser.add_argument('--prefix', dest='prefix', type=str)
    args = parser.parse_args()
    args.device = 'cuda:2'
    args.save_model = True
    print(args)
    # create output directory if not exists
    args.output_dir = os.path.join(result_dir, f"{args.prefix}_{timestamp}")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    train_args = Namespace(
        data=args.data,
        output_dir=args.output_dir,
        prefix="",
        config = None, # use default configuration
        save_model = True,
        device = 'cuda:3' # if GPU is not available, you can also use cpu, but not recommended
    )
    SpaDOT.train(train_args)
    # analyze results
    args.data = f'{args.output_dir}/latent.h5ad'

    latent_adata = anndata.read_h5ad(args.data)
    tp_kmeans_dict = dict()
    ARI_dict = dict()
    for tp in set(latent_adata.obs['timepoint']):
        tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp].copy()
        tp_adata.obs['pixel_x'] = tp_adata.obsm['spatial'][:, 0]
        tp_adata.obs['pixel_y'] = tp_adata.obsm['spatial'][:, 1]
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), random_state=0).fit(tp_adata.X)
        tp_kmeans_dict[tp] = tp_kmeans
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        tp_ARI = adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])
        ARI_dict[tp] = tp_ARI
        print('%s Kmeans ARI: %.3f' % (tp, tp_ARI))
        plt.figure(figsize=(5,5))
        ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='kmeans', s=10)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(args.output_dir+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
        plt.close()
    pd.DataFrame(ARI_dict, index=[0]).to_csv(args.output_dir+os.sep+f'spaDOT_ARI.csv', index=False)