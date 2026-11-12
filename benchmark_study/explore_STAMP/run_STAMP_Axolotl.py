# conda activate /net/zootopia/disk1/wenjinma/envs/stamp
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import numpy as np
import sctm
import squidpy as sq
import pandas as pd
import anndata

# the location of R (used for the mclust clustering)
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_Artista"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_STAMP/Axolotl_results"
os.makedirs(result_dir, exist_ok=True)


# --- load data
#adata = anndata.read_h5ad(data_dir+os.sep+'Development_curated.h5ad')
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

adata.layers['count'] = adata.X.astype('float32')
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
adata.obs['timepoint_numeric'] = adata.obs['timepoint_numeric'].astype('float32')
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()

# --- run STAMP
sq.gr.spatial_neighbors(adata, n_neighs=round(1 / 1000 * adata.n_obs))
sctm.seed.seed_everything(0)
print(adata)
sctm.pp.filter_genes(adata, 0.05)
print(adata)

model = sctm.stamp.STAMP(
    adata,
    n_topics=20,
    layer="count",
    time_covariate_keys="timepoint_numeric",
    # gene_likelihood="poisson",
)
# We used a weighted sampler here as that the first timepoint is extremly small compared to the last timepoint.
model.train(device="cuda:0")
topic_prop = model.get_cell_by_topic()
adata.obsm['X_stamp'] = topic_prop.values
for i in topic_prop.columns:
    adata.obs[i] = topic_prop[i]
adata.obs.to_csv(result_dir+os.sep+'STAMP_topic_proportions.csv')

for tp in timepoints:
    tp_adata = adata[adata.obs['timepoint'] == tp]
    for i in topic_prop.columns:
        sctm.pl.spatial(tp_adata, color=i, size=5, vmax='p99')
        plt.savefig(result_dir+os.sep+f'{i}_{tp}_spatial.png')

# --- use topic proportions to perform KMeans clustering
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in set(adata.obs['timepoint']):
    tp_adata = adata[adata.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['X_stamp']
    tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['region_annotation'])), random_state=0).fit(tp_emb)
    tp_kmeans_dict[tp] = tp_kmeans
    tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
    tp_ARI = adjusted_rand_score(tp_adata.obs['region_annotation'], tp_adata.obs['kmeans'])
    print('%s Kmeans ARI: %.3f' % (tp, tp_ARI))
    ARI_dict[tp] = tp_ARI
    plt.figure(figsize=(5,5))
    ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='kmeans', s=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(result_dir+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
    plt.close()
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'STAMP_ARI.csv', index=False)
