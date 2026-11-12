# conda activate /media/sda/wenjinma/envs/SEDR/
import os
import anndata
import torch
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

import random
from torch.backends import cudnn
seed = 1993
os.environ['PYTHONHASHSEED'] = str(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
random.seed(seed)
cudnn.deterministic = True
cudnn.benchmark = False

import SEDR
import harmonypy as hm

device = torch.device('cuda')
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_MOSTA"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_SEDR/MOSTA_brain_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'adata_brain_clean.h5ad')
n_samples = 20000  # Set desired number of cells
selected_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
adata = adata[selected_indices, :].copy()
adata.X = adata.layers['count'].astype(int)
# focus on first three time points
adata.obs['pixel_x'] = adata.obsm['spatial'][:, 0]
adata.obs['pixel_y'] = adata.obsm['spatial'][:, 1]
adata.obs['annotation'] = 'Brain'
# --- end loading data

# --- SEDR analysis
for idx, tp in enumerate(list(set(adata.obs['timepoint']))):
    tp_adata = adata[adata.obs['timepoint']==tp].copy()
    tp_adata.var_names_make_unique()
    tp_adata.obs['batch_name'] = tp
    graph_dict_tmp = SEDR.graph_construction(tp_adata, 12)
    if idx == 0:
        SEDR_adata = tp_adata
        graph_dict = graph_dict_tmp
        SEDR_adata.obs['proj_name'] = tp
    else:
        var_names = SEDR_adata.var_names.intersection(tp_adata.var_names)
        SEDR_adata = SEDR_adata[:, var_names]
        tp_adata = tp_adata[:, var_names]
        tp_adata.obs['proj_name'] = tp
        SEDR_adata = SEDR_adata.concatenate(tp_adata)
        graph_dict = SEDR.combine_graph_dict(graph_dict, graph_dict_tmp)
SEDR_adata.layers['count'] = SEDR_adata.X
sc.pp.normalize_total(SEDR_adata, target_sum=1e6)
sc.pp.highly_variable_genes(SEDR_adata, flavor="seurat_v3", layer='count', n_top_genes=2000)
SEDR_adata = SEDR_adata[:, SEDR_adata.var['highly_variable'] == True]
sc.pp.scale(SEDR_adata)
# Run PCA on embeddings
from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable.
adata_X = PCA(n_components=200, random_state=42).fit_transform(SEDR_adata.X)
SEDR_adata.obsm['X_pca'] = adata_X
sedr_net = SEDR.Sedr(SEDR_adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
using_dec = False
if using_dec:
    sedr_net.train_with_dec()
else:
    sedr_net.train_without_dec()
sedr_feat, _, _, _ = sedr_net.process()
SEDR_adata.obsm['SEDR'] = sedr_feat
# Use harmony to integrate embeddings from different batches and calculate revised embeddings
meta_data = SEDR_adata.obs[['batch']]
data_mat = SEDR_adata.obsm['SEDR']
vars_use = ['batch']
ho = hm.run_harmony(data_mat, meta_data, vars_use)
res = pd.DataFrame(ho.Z_corr).T
res_df = pd.DataFrame(data=res.values, columns=['X{}'.format(i+1) for i in range(res.shape[1])], index=SEDR_adata.obs.index)
SEDR_adata.obsm[f'SEDR.Harmony'] = res_df
# plot embeddings from SEDR
sc.pp.neighbors(SEDR_adata, use_rep='SEDR.Harmony', metric='cosine') # In the tutorial, cosine is used
sc.tl.umap(SEDR_adata)
sc.pl.umap(SEDR_adata, color=['annotation', 'batch_name'], show=False)
plt.savefig(result_dir+os.sep+'SEDR_umap.png')
plt.close()
SEDR_adata.write(result_dir+os.sep+'SEDR_adata.h5ad')
# explore clustering results revised embeddings
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in set(SEDR_adata.obs['timepoint']):
    tp_adata = SEDR_adata[SEDR_adata.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['SEDR.Harmony']
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
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'SEDR_ARI.csv', index=False)
