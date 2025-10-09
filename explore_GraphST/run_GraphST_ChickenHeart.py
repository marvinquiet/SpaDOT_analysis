# conda activate /media/sda/wenjinma/envs/GraphST/
import os
import torch
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from GraphST import GraphST
from GraphST.utils import clustering

device = torch.device('cuda')
os.environ['R_HOME'] = '/media/sda/wenjinma/envs/GraphST/bin/R'
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_ChickenHeart"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_GraphST/ChickenHeart_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'processed_ChickenHeart.h5ad') # curated region split
adata.obsm['spatial'] = adata.obs[['pixel_x', 'pixel_y']]
adata.obs['annotation'] = adata.obs['region']
# --- end loading data

loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']
# plot algined image
plt.rcParams["figure.figsize"] = (3, 3)
ax = sc.pl.embedding(adata, basis='spatial',
                color='data',
                show=False)
ax.set_title('Aligned image')
plt.tight_layout()
plt.savefig(result_dir+os.sep+'aligned_image.png')
plt.close()

# train GraphST model
model = GraphST.GraphST(adata, device=device)
adata = model.train()

# plot images
fig, ax_list = plt.subplots(1, 3, figsize=(12, 3))
### Plotting UMAP before batch effect correction
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color='data', title='Uncorrected',
                ax = ax_list[0],
                show=False)
### Plotting UMAP after batch effect correction
sc.pp.neighbors(adata, use_rep='emb', n_neighbors=10)
sc.tl.umap(adata)
sc.pl.umap(adata,
        color='data',
        ax=ax_list[1],
        title='Batch corrected',
        #legend_loc = 'bottom margin',
        show=False)
### Color by predicted domains
sc.pl.umap(adata, color='annotation', ax=ax_list[2], title='Colored by region', show=False)
plt.tight_layout(w_pad=0.02)
plt.savefig(result_dir+os.sep+'umap.png')
plt.close()
# save adata for future analysis
adata.write(result_dir+os.sep+'GraphST_adata.h5ad')
# calculate clustering results per timepoint
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in set(adata.obs['timepoint']):
    tp_adata = adata[adata.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['emb']
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
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'GraphST_ARI.csv', index=False)

