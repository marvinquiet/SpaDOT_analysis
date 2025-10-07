import os, math
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
sns.set_palette(cc.glasbey)

from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

#result_dir = "/net/mulan/home/wenjinma/projects/spatialATAC-explore/results/SpatialPCA_GSE214989_separatesample"
#result_dir = "/net/mulan/home/wenjinma/projects/spatialATAC-explore/results/SpatialPCA_GSE214989_multiplesample_combinedlabels"
result_dir = "/net/mulan/home/wenjinma/projects/spatialATAC-explore/results/SpatialPCA_Artista_development_multiplesample_combinedlabels"
pc_df = pd.read_csv(result_dir+os.sep+'spatialPCA_pcs.csv', index_col=0)
final_latent = pc_df.iloc[:, 0:20]

latent_adata = sc.AnnData(final_latent, dtype="float64")
latent_adata.obs = pc_df
latent_adata.obs['annotation'] = latent_adata.obs['annotation'].astype('category')
sc.pp.neighbors(latent_adata, n_neighbors=round(math.sqrt(latent_adata.shape[0])))
sc.tl.umap(latent_adata)
plt.clf()
sc.pl.umap(latent_adata, color=['annotation', 'timepoint'],
        ncols=3)
plt.savefig(result_dir+os.sep+'latent_umap.png')

# run KMeans separately
tp_kmeans_dict = dict()
for tp in set(latent_adata.obs['timepoint']):
    tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp]
    #tp_adata.obs['pixel_x_bak'] = tp_adata.obs['pixel_x']
    #tp_adata.obs['pixel_y_bak'] = tp_adata.obs['pixel_y']
    #if tp == 't1':
    #    tp_adata.obs['pixel_x'] = tp_adata.obs['pixel_y_bak']
    #    tp_adata.obs['pixel_y'] = -tp_adata.obs['pixel_x_bak']
    #elif tp == 't0' or tp == 't2':
    #    tp_adata.obs['pixel_x'] = -tp_adata.obs['pixel_y_bak']
    #    tp_adata.obs['pixel_y'] = -tp_adata.obs['pixel_x_bak']
    tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), random_state=0).fit(tp_adata.X)
    tp_kmeans_dict[tp] = tp_kmeans
    tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
    print('Timepoint %s Seprate KMeans ARI: %.3f' % (tp, adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])))
    plt.figure(figsize=(4, 4))
    ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='kmeans', s=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(result_dir+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
    # separate louvain clustering
    #import igraph as ig
    #from sklearn.neighbors import NearestNeighbors
    #nn_model = NearestNeighbors(n_neighbors=round(math.sqrt(tp_adata.shape[0])) + 1)
    #nn_model.fit(tp_adata.X)
    #_, indices = nn_model.kneighbors(tp_adata.X)
    #snn_graph = ig.Graph()
    #snn_graph.add_vertices(len(tp_adata.X))
    ## Add edges based on shared neighbors
    #for i, neighbors in enumerate(indices):
    #    for neighbor in neighbors[1:]:  # Exclude the point itself
    #        shared_neighbors = set(indices[neighbor][1:])
    #        snn_graph.add_edge(i, neighbor, shared_neighbors=shared_neighbors)
    #g_walk = ig.Graph.community_walktrap(snn_graph)
    #cluster_results = g_walk.as_clustering(len(set(tp_adata.obs['annotation']))).membership
    #tp_adata.obs['walktrap_cluster'] = cluster_results
    #print('Timepoint %s Separate Walktrap ARI: %.3f' % (tp, adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['walktrap_cluster'])))
    #plt.clf()
    #ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='walktrap_cluster', s=10)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.tight_layout()
    #plt.savefig(result_dir+os.sep+str(tp)+'_separate_walktrap_clustering.png')

from OT_loss.ot_solvers import compute_transport_map
ot_config = {
    "growth_iters": 3,
    "ot_epochs": 10,
    "epsilon": 0.05, # larger epsilon can result in a blurred transport plan
    "lambda1": 1, # row imbalance
    "lambda2": 1, # column imbalance
    "epsilon0": 1,
    "tau": 1000,
    "scaling_iter": 3000,
    "inner_iter_max": 50,
    "tolerance": 1e-8,
    "max_iter": 1e7,
    "batch_size": 5,
    "extra_iter": 1000,
    "numItermax": 1000000,
    "use_Py": False,
    "use_C": True,
    "profiling": False,
    "method": "waddington"
}  # specified in tutorial
 
timepoints = list(set(latent_adata.obs['timepoint']))
timepoints.sort()
for tp_i in range(len(timepoints)-1):
    cur_tp = timepoints[tp_i]
    next_tp = timepoints[tp_i+1]
    gamma = compute_transport_map(tp_kmeans_dict[cur_tp].cluster_centers_, tp_kmeans_dict[next_tp].cluster_centers_, ot_config, G=None)  ## do not set growth rate for now
    ## write OT matrix
    np.save(result_dir+os.sep+'gamma_tp{}_tp{}_final_OT.npy'.format(cur_tp, next_tp), gamma)
    plt.clf()
    plt.figure(figsize=(4, 4))
    ax = sns.heatmap(gamma, linewidth=0.5)
    plt.savefig(result_dir+os.sep+'gamma_tp{}_tp{}_final_heatmap.png'.format(cur_tp, next_tp))
 
