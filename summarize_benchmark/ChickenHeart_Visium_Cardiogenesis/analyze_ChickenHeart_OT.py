import os
import anndata
import numpy as np
import pandas as pd
import wot
import matplotlib.pyplot as plt

base_dir = "/media/sda/wenjinma"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_GATVAE-SVGP-VAE/ChickenHeart_results"
latent_adata = anndata.read_h5ad(result_dir+os.sep+'prediction_adata.h5ad')
latent_adata.obs['day'] = latent_adata.obs['timepoint_numeric']
latent_adata.obs['cell_growth_rate'] = 1 # initialize as 1
ot_model = wot.ot.OTModel(latent_adata, epsilon = 0.01, lambda1 = 0.1,lambda2 = 5, growth_iters=3) # imbalanced OT with imbalanced row and comparatiely balanced target
ot_model.compute_all_transport_maps(tmap_out=result_dir+os.sep+'ChickenHeart') # compute transport maps
tmap_model = wot.tmap.TransportMapModel.from_directory(result_dir+os.sep+'ChickenHeart')
# generate region dict
latent_adata.obs['spaDOT_pred_labels'] = latent_adata.obs['timepoint'].astype('str')+'_'+latent_adata.obs['pred_labels'].astype('str')
latent_cell_sets = latent_adata.obs.groupby('spaDOT_pred_labels').apply(lambda x: x.index.tolist()).to_dict()
days = latent_adata.obs['day'].unique()
days.sort()
for tp_i in range(len(days)-1):
    prev_day = days[tp_i]
    next_day = days[tp_i+1]
    # obtain 
    prob_adata = anndata.read_h5ad(result_dir+os.sep+'ChickenHeart_'+str(prev_day)+'_'+str(next_day)+'.h5ad')
    prev_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=prev_day)
    next_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=next_day)
    # for i in range(len(prev_day_populations)):
    #     plt.figure()
    #     for j in range(len(next_day_populations)):
    #         prob_i_j = prob_adata.X[prev_day_populations[i].p.astype(bool)][:, next_day_populations[j].p.astype(bool)]
    #         prob_i_j_flattened = np.log(prob_i_j.flatten())
    #         prob_i_j_flattened[np.isneginf(prob_i_j_flattened)] = np.nan
    #         plt.hist(prob_i_j_flattened, bins=100, alpha=0.5, label='from '+str(i)+' to '+str(j))
    #     plt.legend()
    #     plt.savefig(result_dir+os.sep+'transition_prob_'+str(prev_day)+'_'+str(next_day)+'_'+str(i)+'.png')
    #     plt.close()

    transition_table = tmap_model.transition_table(prev_day_populations, next_day_populations) # aggregated OT matrix
    transition_table.write_h5ad(result_dir+os.sep+'transition_table_'+str(prev_day)+'_'+str(next_day)+'.h5ad')
    # analyze transition table
    # --- normalize by column sum
    transition_prob = transition_table.X
    transition_prob = transition_prob/transition_prob.sum(axis=0, keepdims=True)
    transition_prob_df = pd.DataFrame(transition_prob, index=transition_table.obs_names, columns=transition_table.var_names)
    transition_prob_df.to_csv(result_dir+os.sep+'transition_prob_'+str(prev_day)+'_'+str(next_day)+'_col_norm.csv')
    # --- normalize by row sum
    transition_prob = transition_table.X
    transition_prob = transition_prob/transition_prob.sum(axis=1, keepdims=True)
    transition_prob_df = pd.DataFrame(transition_prob, index=transition_table.obs_names, columns=transition_table.var_names)
    transition_prob_df.to_csv(result_dir+os.sep+'transition_prob_'+str(prev_day)+'_'+str(next_day)+'_row_norm.csv')
