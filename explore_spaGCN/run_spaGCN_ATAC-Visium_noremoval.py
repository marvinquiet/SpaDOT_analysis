# conda activate /media/sda/wenjinma/envs/SpaGCN_py39/
import warnings
warnings.filterwarnings("ignore")

import SpaGCN as spg

# the location of R (used for the mclust clustering)
import os, argparse
from anndata import AnnData
import scanpy as sc
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

import anndata

import torch
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaGCN/ATAC-Visium_noremoval_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
section_dict = {1: 't2', 2: 't1', 3: 't0'}
data_dict = {1: 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz', # E15.5
             2: 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz', # E13.5
             3: 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz' # E12.5
            }
adata = anndata.read_h5ad(data_dir+os.sep+'exprs.h5ad')
adata.obs['x'] = None
adata.obs['y'] = None
adata.obs['pixel_x'] = None
adata.obs['pixel_y'] = None
loc_metadata_list = []
for key, value in data_dict.items():
    metadata = pd.read_csv(data_dir+os.sep+value, sep=',', index_col=0, header=None)
    metadata.columns = ['tissue', 'y', 'x', 'pixel_y', 'pixel_x']
    metadata['barcodes'] = metadata.index
    metadata.index = metadata['barcodes']+'_'+str(key)
    key_barcodes = adata[adata.obs['section'] == key].obs_names.tolist()
    common_barcodes = list(set(metadata.index.tolist()).intersection(set(key_barcodes)))
    section_metadata = metadata.loc[common_barcodes, :]
    loc_metadata_list.append(section_metadata)
loc_metadata = pd.concat(loc_metadata_list, axis=0)
common_barcodes = list(set(loc_metadata.index.tolist()).intersection(set(adata.obs_names.tolist())))
adata = adata[common_barcodes].copy()
adata.obs.loc[common_barcodes, 'x'] = loc_metadata.loc[common_barcodes, 'x'].values
adata.obs.loc[common_barcodes, 'y'] = loc_metadata.loc[common_barcodes, 'y'].values
adata.obs.loc[common_barcodes, 'pixel_x'] = loc_metadata.loc[common_barcodes, 'pixel_x'].values
adata.obs.loc[common_barcodes, 'pixel_y'] = loc_metadata.loc[common_barcodes, 'pixel_y'].values
adata.obs['annotation'] = adata.obs['seurat_clusters_harmony']
adata.obs['annotation'].replace({3: 'CNS', 2: 'CNS', 4: 'CNS',
                                 1: 'Men/PNS', 9: 'Men/PNS', 11: 'Men/PNS',
                                 0: 'Mesenchyme', 5: 'Mesenchyme', 6: 'Mesenchyme', 8: 'Mesenchyme',
                                 7: 'Liver', 10: 'Heart'}, inplace=True)
adata.obs['annotation'] = adata.obs['annotation'].astype('category')
adata.obs['timepoint'] = [section_dict[_] for _ in adata.obs['section']]
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
# remove those cells that are not in the image
adata.obs.to_csv(result_dir+os.sep+'combined_ATAC-Visium_metadata.csv')
print("Data size:", adata.shape)
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']
# --- end loading data

timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
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
