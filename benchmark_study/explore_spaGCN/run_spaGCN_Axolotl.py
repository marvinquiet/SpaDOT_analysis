# conda activate /media/sda/wenjinma/envs/SpaGCN_py39/
import warnings
warnings.filterwarnings("ignore")

import SpaGCN as spg

import os
from anndata import AnnData
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

import anndata

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_Artista"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaGCN/Axolotl_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
#adata = anndata.read_h5ad(data_dir+os.sep+'Development_curated.h5ad') # curated region split
adata = anndata.read_h5ad(data_dir+os.sep+'Development_addCP.h5ad')
adata.X = adata.layers['counts'].astype(int)
adata.var['original_genes'] = adata.var_names.tolist()
adata.var_names = adata.var['Gene'].values
adata.obs['status'] = [_.split('_')[0] for _ in  adata.obs['Batch']]
adata.obs['timepoint'] = adata.obs['status'].replace({'Stage44': 't0',
    'Stage54': 't1', 'Stage57': 't2', 'Injury': 't3', 'Batch1': 't4', 'Meta': 't5'})
# focus on first three time points
timepoints = ['t0', 't1', 't2']
adata = adata[adata.obs['timepoint'].isin(timepoints)].copy()
adata.obs['pixel_x'] = adata.obsm['spatial'][:, 0]
adata.obs['pixel_y'] = -adata.obsm['spatial'][:, 1]
adata.obs_names = adata.obs_names+'_'+adata.obs['timepoint']
adata.obs['annotation'] = adata.obs['region_annotation']
# --- end loading data

# split adata into batches
adata_list, l_list, adj_list = [], [], []
for tp in timepoints:
    adata_tp = adata[adata.obs['timepoint']==tp].copy()
    pixel_x = adata_tp.obs['pixel_x'].tolist()
    pixel_y = adata_tp.obs['pixel_y'].tolist()
    s=1
    b=49
    adj_tp=spg.calculate_adj_matrix(x=pixel_x,y=pixel_y, beta=b, alpha=s, histology=False)
    adata_tp.raw = adata_tp
    sc.pp.normalize_total(adata_tp)
    sc.pp.log1p(adata_tp)
    p = 0.5
    l_tp = spg.search_l(p, adj_tp, start=0.01, end=1000, tol=0.01, max_run=100)
    adj_list.append(adj_tp)
    adata_list.append(adata_tp)
    l_list.append(l_tp)
# run spaGCN multi-section
res=0.6
#Set seed
r_seed=t_seed=n_seed=100
#Multi-section SpaGCN
clf=spg.multiSpaGCN()
clf.train(adata_list,adj_list,l_list,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
# explore spaGCN embeddings
spaGCN_embed = clf.embed
merge_data = AnnData.concatenate(*adata_list)
merge_data.obsm['spaGCN_embed'] = spaGCN_embed
sc.pp.neighbors(merge_data, use_rep="spaGCN_embed", random_state=42)
sc.tl.umap(merge_data, random_state=42)
sc.pl.umap(merge_data, color=["batch", "annotation"])
plt.savefig(result_dir+os.sep+'spaGCN_UMAP.png')
plt.close()
merge_data.write_h5ad(result_dir+os.sep+'spaGCN.h5ad')
# --- calculate ARI per timepoint
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in timepoints:
    tp_adata = merge_data[merge_data.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['spaGCN_embed']
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
    plt.close()
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'spaGCN_ARI.csv', index=False)

