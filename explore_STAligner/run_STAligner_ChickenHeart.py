# conda activate /media/sda/wenjinma/envs/STAligner_py38/
import warnings
warnings.filterwarnings("ignore")

import STAligner

import os
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import scipy.linalg
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

import torch
used_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_ChickenHeart"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_STAligner/ChickenHeart_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'processed_ChickenHeart.h5ad') # curated region split
adata.obsm['spatial'] = adata.obs[['pixel_x', 'pixel_y']]
adata.obs['annotation'] = adata.obs['region']
# --- end loading data

# STAligner analysis
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
# build graph information for each timepoint
Batch_list = []
adj_list = []
for idx, tp in enumerate(timepoints):
    tp_adata = adata[adata.obs['timepoint']==tp].copy()
    tp_adata.X = csr_matrix(tp_adata.X) # covert to csr_matrix, because STAligner.Cal_Spatial_Net only accept sparse matrix
    STAligner.Cal_Spatial_Net(tp_adata, rad_cutoff=1.3)
    # Normalization
    sc.pp.normalize_total(tp_adata, target_sum=1e4)
    sc.pp.log1p(tp_adata)
    sc.pp.highly_variable_genes(tp_adata, flavor="seurat_v3", n_top_genes=5000) #ensure enough common HVGs in the combined matrix
    tp_adata = tp_adata[:, tp_adata.var['highly_variable']]
    adj_list.append(tp_adata.uns['adj'])
    Batch_list.append(tp_adata)
# combine all timepoints
adata_concat = anndata.concat(Batch_list, label="slice_name", keys=timepoints)
adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')
print('adata_concat.shape: ', adata_concat.shape)
adj_concat = np.asarray(adj_list[0].todense())
for batch_id in range(1, len(timepoints)):
    adj_concat = scipy.linalg.block_diag(adj_concat, np.asarray(adj_list[batch_id].todense()))
adata_concat.uns['edgeList'] = np.nonzero(adj_concat)
# use last timepoint as reference
iter_comb = [(0, 2), (1, 2)] ## Fix slice 3 as reference to align
adata_concat = STAligner.train_STAligner(adata_concat, verbose=True, knn_neigh = 100, iter_comb = iter_comb, margin=2.5,  device=used_device)
adata_result = adata_concat.copy()
adata_result.uns.clear()
adata_result.write_h5ad(result_dir+os.sep+'STAligner_adata.h5ad')
# --- explore kmeans clustering on the STAligner embedding
sc.pp.neighbors(adata_concat, use_rep='STAligner', random_state=666)
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in timepoints:
    tp_adata = adata_concat[adata_concat.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['STAligner']
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
    plt.savefig(result_dir+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'STAligner_ARI.csv', index=False)
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['annotation', 'batch_name'], show=False)
plt.savefig(result_dir+os.sep+'STAligner_UMAP.png')
plt.close()
# --- explore kmeans clustering on the STAGATE embedding
sc.pp.neighbors(adata_concat, use_rep='STAGATE', random_state=666)
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in timepoints:
    tp_adata = adata_concat[adata_concat.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['STAGATE']
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
    plt.savefig(result_dir+os.sep+str(tp)+'_Kmeans_separate_clustering_STAGATE.png')
    plt.close()
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'STAGATE_ARI.csv', index=False)
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['annotation', 'batch_name'], show=False)
plt.savefig(result_dir+os.sep+'STAGATE_UMAP.png')
plt.close()