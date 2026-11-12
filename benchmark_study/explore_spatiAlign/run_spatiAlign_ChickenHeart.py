# conda activate /media/sda/wenjinma/envs/STAligner_py38
import warnings
warnings.filterwarnings("ignore")

from spatialign import Spatialign

# the location of R (used for the mclust clustering)
import os
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_ChickenHeart"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spatiAlign/ChickenHeart_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
adata = anndata.read_h5ad(data_dir+os.sep+'processed_ChickenHeart.h5ad') # curated region split
adata.obsm['spatial'] = np.array(adata.obs[['pixel_x', 'pixel_y']])
adata.obs['annotation'] = adata.obs['region']
# --- end loading data

# spatiAlign analysis
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
# split adata into batches
datasets = []
for tp in timepoints:
    adata_tp = adata[adata.obs['timepoint'] == tp].copy()
    adata_tp_file = result_dir+os.sep+str(tp)+'.h5ad'
    adata_tp.write_h5ad(adata_tp_file)
    datasets.append(adata_tp_file)

model = Spatialign(
    *datasets,
    batch_key='batch',
    is_norm_log=True,
    is_scale=False,
    n_neigh=15,
    is_undirected=True,
    latent_dims=100,
    seed=42,
    gpu=0,
    save_path=result_dir+os.sep,
    is_verbose=True
) # use default tau1, tau2 and tau3
model.train()
model.alignment()
correct_adata_list = list()
for i in range(len(timepoints)):
    correct_adata = sc.read_h5ad(result_dir+os.sep+f"res/correct_data{i}.h5ad")
    correct_adata.var['genes'] = correct_adata.var_names.tolist()
    correct_adata_list.append(correct_adata)
merge_data = anndata.AnnData.concatenate(*correct_adata_list)
# merge_data = correct_t1.concatenate(correct_t2, correct_t3)
sc.pp.neighbors(merge_data, use_rep="correct", random_state=42)
sc.tl.umap(merge_data, random_state=42)
sc.pl.umap(merge_data, color=["batch", "annotation"])
plt.savefig(result_dir+os.sep+'spatiAlign_UMAP.png')
plt.close()
# --- calculate ARI per timepoint
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in timepoints:
    tp_adata = merge_data[merge_data.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['correct']
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
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'spatiAlign_ARI.csv', index=False)
