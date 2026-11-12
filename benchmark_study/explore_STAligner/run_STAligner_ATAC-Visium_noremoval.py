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
data_dir = base_dir+os.sep+"data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_STAligner/ATAC-Visium_noremoval_results"
os.makedirs(result_dir, exist_ok=True)


# --- load data
section_dict = {1: 't2', 2: 't1', 3: 't0'}
data_dict = {1: 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz', # E15.5
             2: 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz', # E13.5
             3: 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz' # E12.5
            }
adata = sc.read_h5ad(data_dir+os.sep+'exprs.h5ad')
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
adata.obs.to_csv(result_dir+os.sep+'combined_ATAC-Visium_metadata.csv')
print("Data size:", adata.shape)
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']

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
adata_result.obs['x'] = adata_result.obs['x'].astype(int)
adata_result.obs['y'] = adata_result.obs['y'].astype(int)
adata_result.obs['pixel_x'] = adata_result.obs['pixel_x'].astype(int)
adata_result.obs['pixel_y'] = adata_result.obs['pixel_y'].astype(int)
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
adata_concat.obs['tp_annotation'] = adata_concat.obs['timepoint'].astype(str) + '_' + adata_concat.obs['annotation'].astype(str)
sc.pl.umap(adata_concat, color=['tp_annotation', 'batch_name'], show=False)
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
adata_concat.obs['tp_annotation'] = adata_concat.obs['timepoint'].astype(str) + '_' + adata_concat.obs['annotation'].astype(str)
sc.pl.umap(adata_concat, color=['tp_annotation', 'batch_name'], show=False)
plt.savefig(result_dir+os.sep+'STAGATE_UMAP.png')
plt.close()
