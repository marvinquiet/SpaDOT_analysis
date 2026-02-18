import os
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import wot
import harmonypy as hm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_rel

result_dir = "/net/zootopia/disk1/wenjinma/minicluster/minicluster_backup_06102025/wenjinma/projects/spatialATAC-explore"

def run_OT(latent_adata, prefix=''):
    ot_model = wot.ot.OTModel(latent_adata, epsilon = 0.01, lambda1 = 0.1,lambda2 = 5, growth_iters=3) # imbalanced OT with imbalanced row and comparatiely balanced target
    ot_model.compute_all_transport_maps(tmap_out=prefix) # compute transport maps
    tmap_model = wot.tmap.TransportMapModel.from_directory(prefix)
    # calculate transition matrix
    latent_cell_sets = latent_adata.obs.groupby('true_labels').apply(lambda x: x.index.tolist()).to_dict()
    days = latent_adata.obs['day'].unique()
    days.sort()
    for tp_i in range(len(days)-1):
        prev_day = days[tp_i]
        next_day = days[tp_i+1]
        prev_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=prev_day)
        next_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=next_day)
        transition_table = tmap_model.transition_table(prev_day_populations, next_day_populations) # aggregated OT matrix
        transition_table.write_h5ad(f'{prefix}_transition_table_{prev_day}_{next_day}.h5ad')
        # --- normalize by column sum
        transition_prob = transition_table.X
        transition_prob_col = transition_prob/transition_prob.sum(axis=0, keepdims=True)
        transition_prob_col_df = pd.DataFrame(transition_prob_col, index=transition_table.obs_names, columns=transition_table.var_names)
        transition_prob_col_df.to_csv(f'{prefix}_transition_table_{prev_day}_{next_day}_col_norm.csv')
        # --- normalize by row sum
        transition_prob = transition_table.X
        transition_prob_row = transition_prob/transition_prob.sum(axis=1, keepdims=True)
        transition_prob_row_df = pd.DataFrame(transition_prob_row, index=transition_table.obs_names, columns=transition_table.var_names)
        transition_prob_row_df.to_csv(f'{prefix}_transition_table_{prev_day}_{next_day}_row_norm.csv')
        OT_mat = np.minimum(transition_prob_col_df, transition_prob_row_df)
        OT_mat.to_csv(f'{prefix}_transition_prob_{prev_day}_{next_day}.csv')

# --- analyze simulation results
for dataset in ['ATAC-Visium_noremoval_results', 'Axolotl_results', 'ChickenHeart_results']:
    dataset_dir = f"{result_dir}/explore_GATVAE-SVGP-VAE/{dataset}"
    latent_adata = anndata.read_h5ad(f"{dataset_dir}/prediction_adata.h5ad")
    latent_adata.obsm['X_latent'] = latent_adata.X
    sc.pp.neighbors(
        latent_adata,
        use_rep="X_latent",
    )
    sc.tl.umap(latent_adata)
    sc.pl.umap(
        latent_adata,
        color="true_labels",
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"{dataset_dir}/latent_UMAP.pdf")
    ho = hm.run_harmony(latent_adata.obsm['X_latent'], latent_adata.obs, vars_use=['timepoint'])
    res = pd.DataFrame(ho.Z_corr.T)
    res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
    res.to_csv(f"{dataset_dir}/harmony_latent.csv", sep = "\t", index = False)
    latent_adata.obsm['X_latent_harmony']= ho.Z_corr.T
    sc.pp.neighbors(
        latent_adata,
        use_rep="X_latent_harmony",
    )
    sc.tl.umap(latent_adata)
    sc.pl.umap(
        latent_adata,
        color="true_labels",
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"{dataset_dir}/latent_UMAP_harmony.pdf")
    latent_adata.obs['day'] = latent_adata.obs['timepoint_numeric']
    latent_adata.obs['cell_growth_rate'] = 1 # initialize as 1
    latent_adata.X = latent_adata.obsm['X_latent']
    run_OT(latent_adata, prefix=dataset_dir+os.sep+"original")
    latent_adata.X = latent_adata.obsm['X_latent_harmony']
    run_OT(latent_adata, prefix=dataset_dir+os.sep+"harmony")
    days = latent_adata.obs['day'].unique()
    days.sort()
    for tp_i in range(len(days)-1):
        prev_day = days[tp_i]
        next_day = days[tp_i+1]
        harmony_OT_mat = pd.read_csv(dataset_dir+os.sep+f"harmony_transition_prob_{prev_day}_{next_day}.csv", index_col=0)
        original_OT_mat = pd.read_csv(dataset_dir+os.sep+f"original_transition_prob_{prev_day}_{next_day}.csv", index_col=0)
        # flatten the matrices
        mat1_flat = original_OT_mat.values.flatten()
        mat2_flat = harmony_OT_mat.values.flatten()
        stat, pvalue = ttest_rel(mat1_flat, mat2_flat)
        print("p-value:", pvalue)
        diff_sum = np.sum(np.abs(mat1_flat - mat2_flat))
        print("Abs diff: ", diff_sum)
        print("Effect size:", diff_sum / np.sum(mat1_flat))
        plt.figure(figsize=(8, 8))
        sns.heatmap(original_OT_mat, linewidth=0.5, annot=True, cmap="Reds", fmt=".1e", vmin=0, vmax=1)
        plt.savefig(dataset_dir+os.sep+f"original_gamma_{prev_day}_{next_day}_final_heatmap.png")
        plt.tight_layout()
        plt.close()
        plt.figure(figsize=(8, 8))
        sns.heatmap(harmony_OT_mat, linewidth=0.5, annot=True, cmap="Reds", fmt=".1e", vmin=0, vmax=1)
        plt.savefig(dataset_dir+os.sep+f"harmony_gamma_{prev_day}_{next_day}_final_heatmap.png")
        plt.tight_layout()
        plt.close()

