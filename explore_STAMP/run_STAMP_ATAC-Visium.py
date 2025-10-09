# conda activate /net/zootopia/disk1/wenjinma/envs/stamp
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import numpy as np
import sctm
import squidpy as sq
import pandas as pd
import seaborn as sns
import anndata as ad

# the location of R (used for the mclust clustering)
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


base_dir = "/net/mulan/home/wenjinma"
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_STAMP/ATAC-Visium_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
section_dict = {1: 't2', 2: 't1', 3: 't0'}
data_dict = {1: 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz', # E15.5
             2: 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz', # E13.5
             3: 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz' # E12.5
            }
adata = ad.read_h5ad(data_dir+os.sep+'exprs.h5ad')
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
with open(data_dir+os.sep+'t2_removed_cells.txt', 'r') as f:
    t2_removed_cells = f.readlines()
t2_removed_cells = [_.strip() for _ in t2_removed_cells]
adata = adata[~adata.obs_names.isin(t2_removed_cells)]
adata.obs.to_csv(result_dir+os.sep+'combined_ATAC-Visium_metadata.csv')
print("Data size:", adata.shape)
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']
adata.layers['count'] = adata.X.astype('float32')
# --- end loading data

timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
# --- run STAMP
sq.gr.spatial_neighbors(adata, n_neighs=round(1 / 1000 * adata.n_obs))
sctm.seed.seed_everything(0)
print(adata)

# --- select top 6,000 highly variable genes as suggested in the paper -> does not work, would generate problem
# sc.pp.highly_variable_genes(adata, n_top_genes=6000, flavor="seurat_v3", layer="count")
# adata = adata[:, adata.var.highly_variable].copy()
sctm.pp.filter_genes(adata, 0.05)
adata.obs['timepoint_numeric'] = adata.obs['timepoint_numeric'].astype('float32')

# --- for this dataset, each slice contains 5 clusters, so we set n_topics=5
model = sctm.stamp.STAMP(
    adata,
    n_topics=5,
    layer="count",
    time_covariate_keys="timepoint_numeric",
    # gene_likelihood="poisson",
)

# We used a weighted sampler here as that the first timepoint is extremly small compared to the last timepoint.
model.train(device="cuda:0") # here I did not use the sampler because the size is the same, I will use on other datasets
topic_prop = model.get_cell_by_topic()
beta = model.get_feature_by_topic()
for i in topic_prop.columns:
    adata.obs[i] = topic_prop[i]


for tp in timepoints:
    tp_adata = adata[adata.obs['timepoint'] == tp]
    for i in topic_prop.columns:
        sctm.pl.spatial(tp_adata, color=i, size=5, vmax='p99')
        plt.savefig(result_dir+os.sep+f'{i}_{tp}_spatial.png')


for topic in topic_prop.columns:
    genes = []
    for i in range(len(timepoints)):
        genes += beta.xs(topic, level=1).transpose().nlargest(10, i).index.tolist()
    genes = np.unique(genes)
    beta.xs(topic, level=1).rank(ascending=False, axis=1).transpose().loc[genes]
    sns.clustermap(
        beta.xs(topic, level=1).transpose().loc[genes],
        col_cluster=False,
        cmap="RdBu_r",
        vmin=-3,
        xticklabels=True,
        yticklabels=True,
        annot=True,
    )
    plt.savefig(result_dir+os.sep+f'{topic}_heatmap.png')


# ---  I will use topic = 20 and the topic proportion take as latent space
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

# --- use topic proportions to perform KMeans clustering
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in set(adata.obs['timepoint']):
    tp_adata = adata[adata.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['X_stamp']
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
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'STAMP_ARI.csv', index=False)
