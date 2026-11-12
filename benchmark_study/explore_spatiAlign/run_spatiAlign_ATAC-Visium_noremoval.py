# conda activate /media/sda/wenjinma/envs/STAligner_py38
import warnings
warnings.filterwarnings("ignore")

from spatialign import Spatialign

# the location of R (used for the mclust clustering)
import os
import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spatiAlign/ATAC-Visium_noremoval_results"
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
adata.obs['batch'] = adata.obs['timepoint']
adata.obs['x'] = adata.obs['x'].astype(np.float)
adata.obs['y'] = adata.obs['y'].astype(np.float)
adata.obs['pixel_x'] = adata.obs['pixel_x'].astype(np.float)
adata.obs['pixel_y'] = adata.obs['pixel_y'].astype(np.float)

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
correct_t1 = sc.read_h5ad(result_dir+os.sep+"res/correct_data0.h5ad")
correct_t2 = sc.read_h5ad(result_dir+os.sep+"res/correct_data1.h5ad")
correct_t3 = sc.read_h5ad(result_dir+os.sep+"res/correct_data2.h5ad")
# merge_data = correct_t1.concatenate(correct_t2, correct_t3)
merge_data = ad.concat([correct_t1, correct_t2, correct_t3])
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
