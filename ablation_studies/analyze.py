import os
import anndata

import sys
sys.path.append('/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies')
from utils import _analyze_utils

def analyze(args):
    data_dir = os.path.abspath(args.data)
    if not args.output_dir:
        args.output_dir = os.path.dirname(data_dir)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if args.n_clusters is None and args.prefix == "":
        args.prefix = 'adaptive_'

    latent = anndata.read_h5ad(data_dir)
    if args.n_clusters is None:
        latent = _analyze_utils.Adaptive_clustering(args, latent)
    else:
        latent = _analyze_utils.KMeans_Clustering(latent, args.n_clusters)
    latent.obs['pixel_x'] = latent.obsm['spatial'][:, 0]
    latent.obs['pixel_y'] = latent.obsm['spatial'][:, 1]
    # draw domains
    _analyze_utils.plot_domains(args, latent)
    # perform optimal transport analysis
    _analyze_utils.OT_analysis(args, latent)
    # plot OT results
    _analyze_utils.plot_OT(args, latent)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="SpaDOT ablation studies analysis.")
    parser.add_argument('--output_dir', dest='output_dir', type=str)
    parser.add_argument('--prefix', dest='prefix', type=str)
    args = parser.parse_args()
    args.data = f'{args.output_dir}/{args.prefix}latent.h5ad'

    # create arguments for testing
    # class Args:
    #     prefix = ""
    # args = Args()
    # args.prefix = 'MouseOrganogenesis_'
    # args.data = f'SpaDOT/output/{args.prefix}latent.h5ad'
    # args.n_clusters = [5, 5, 5]
    # create output directory if not exists
    if 'output_dir' not in args.__dict__:
        args.output_dir = os.path.dirname(args.data)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # analyze latent representations

    import anndata
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score
    latent_adata = anndata.read_h5ad(args.data)
    tps = list(latent_adata.obs['timepoint'].unique())
    ARI_dict = dict()
    n_clusters = []
    for i, tp in enumerate(tps):
        tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp].copy()
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), 
                           random_state=1993, 
                           n_init=10).fit(tp_adata.X)
        n_clusters.append(len(set(tp_adata.obs['annotation'])))
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        tp_ARI = adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])
        ARI_dict[tp] = tp_ARI
    pd.DataFrame(ARI_dict, index=[0]).to_csv(args.output_dir+os.sep+f'{args.prefix}spaDOT_ARI.csv', index=False)
    args.n_clusters = n_clusters
    analyze(args)
