import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import seaborn as sns

from utils import *
import scipy.io

## check with our simulation data
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_STRS"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_SPIRAL/STRS_results"
os.makedirs(result_dir, exist_ok=True)

# --- load data
exprs = scipy.io.mmread(data_dir+os.sep+'STRS.mtx.gz').toarray()
genes = np.loadtxt(data_dir+os.sep+'STRS_genes.txt', dtype='str')
metadata = pd.read_csv(data_dir+os.sep+'STRS_metadata.csv', header=0, index_col=0)
adata = anndata.AnnData(X=exprs.T, obs=metadata, var=pd.DataFrame(index=genes, data={'gene_id':genes}))
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
adata.obsm['spatial'] = loc
adata.obs['annotation'] = adata.obs['injury.zones']
timepoints = list(set(adata.obs['timepoint']))
timepoints.sort()
# --- end loading data

# preprocess data and train SPIRAL
Preprocess(adata, result_dir=result_dir) # default KNN=10, rad=10
train_SPIRAL(result_dir, timepoints) # default KNN=10
spiral_adata = sc.read_h5ad(result_dir+os.sep+"SPIRAL.h5ad")
spiral_adata.obs = adata.obs.loc[spiral_adata.obs_names, :]
spiral_adata.write_h5ad(result_dir+os.sep+"SPIRAL_corrected.h5ad")
# --- explore kmeans clustering on the STAligner embedding
tp_kmeans_dict = dict()
ARI_dict = dict()
for tp in timepoints:
    tp_adata = spiral_adata[spiral_adata.obs['timepoint'] == tp]
    tp_emb = tp_adata.obsm['spiral']
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
pd.DataFrame(ARI_dict, index=[0]).to_csv(result_dir+os.sep+'SPIRAL_ARI.csv', index=False)
sc.tl.umap(spiral_adata)
sc.pl.umap(spiral_adata, color=["timepoint", "annotation"])
plt.savefig(result_dir+os.sep+"SPIRAL_umap.png")
plt.close()