import os
import anndata
import numpy as np
import pandas as pd
import wot

base_dir = "/media/sda/wenjinma"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_GATVAE-SVGP-VAE"

# analyze results from seed 1993
sim_combinations = ['1993/20_30', '1993/20_50', '1993/30_30', '1993/30_50']
for sim_comb in sim_combinations:
    result_prefix = result_dir+os.sep+"simulation_results"+os.sep+sim_comb
    spaDOT_result_prefix = result_prefix+os.sep+'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1'
    latent_adata = anndata.read_h5ad(spaDOT_result_prefix+os.sep+'prediction_adata.h5ad')
    latent_adata.obs['day'] = latent_adata.obs['timepoint_numeric']
    latent_adata.obs['cell_growth_rate'] = 1 # initialize as 1
    ot_model = wot.ot.OTModel(latent_adata, epsilon = 0.01, lambda1 = 0.1,lambda2 = 5, growth_iters=3) # imbalanced OT with imbalanced row and comparatiely balanced target
    ot_model.compute_all_transport_maps(tmap_out=spaDOT_result_prefix+os.sep+'simulation') # compute transport maps
    tmap_model = wot.tmap.TransportMapModel.from_directory(spaDOT_result_prefix+os.sep+'simulation')
    # generate region dict
    latent_adata.obs['spaDOT_pred_labels'] = latent_adata.obs['timepoint'].astype('str')+'_'+latent_adata.obs['pred_labels'].astype('str')
    latent_cell_sets = latent_adata.obs.groupby('spaDOT_pred_labels').apply(lambda x: x.index.tolist()).to_dict()
    days = latent_adata.obs['day'].unique()
    days.sort()
    for tp_i in range(len(days)-1):
        prev_day = days[tp_i]
        next_day = days[tp_i+1]
        prev_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=prev_day)
        next_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=next_day)
        transition_table = tmap_model.transition_table(prev_day_populations, next_day_populations) # aggregated OT matrix
        transition_table.write_h5ad(spaDOT_result_prefix+os.sep+'transition_table_'+str(prev_day)+'_'+str(next_day)+'.h5ad')
        # analyze transition table
        transition_prob = transition_table.X
        transition_prob = transition_prob/transition_prob.sum(axis=0)
        transition_prob_df = pd.DataFrame(transition_prob, index=transition_table.obs_names, columns=transition_table.var_names)
        transition_prob_df.to_csv(spaDOT_result_prefix+os.sep+'transition_prob_'+str(prev_day)+'_'+str(next_day)+'.csv')