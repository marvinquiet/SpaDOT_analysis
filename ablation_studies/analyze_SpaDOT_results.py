# conda activate ~/envs/SpaDOT
import os
import anndata
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/ablation_studies/test_SpaDOT_package"

for folder in os.listdir(result_dir):
    subfolder_path = os.path.join(result_dir, folder)
    if os.path.isdir(subfolder_path):  # ensure it's a folder
        target_latent_path = os.path.join(subfolder_path, 'latent.h5ad')
        if os.path.exists(target_latent_path):
            latent_adata = anndata.read_h5ad(target_latent_path)
            tp_kmeans_dict = dict()
            ARI_dict = dict()
            for tp in set(latent_adata.obs['timepoint']):
                tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp].copy()
                tp_adata.obs['pixel_x'] = tp_adata.obsm['spatial'][:, 0]
                tp_adata.obs['pixel_y'] = tp_adata.obsm['spatial'][:, 1]
                tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs['annotation'])), random_state=0).fit(tp_adata.X)
                tp_kmeans_dict[tp] = tp_kmeans
                tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
                tp_ARI = adjusted_rand_score(tp_adata.obs['annotation'], tp_adata.obs['kmeans'])
                ARI_dict[tp] = tp_ARI
                print('%s Kmeans ARI: %.3f' % (tp, tp_ARI))
                plt.figure(figsize=(5,5))
                ax = sns.scatterplot(data=tp_adata.obs, x='pixel_x', y='pixel_y', hue='kmeans', s=10)
                plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                plt.tight_layout()
                plt.savefig(subfolder_path+os.sep+str(tp)+'_Kmeans_separate_clustering.png')
                plt.close()
            pd.DataFrame(ARI_dict, index=[0]).to_csv(subfolder_path+os.sep+f'spaDOT_ARI.csv', index=False)
