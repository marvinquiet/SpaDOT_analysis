# conda activate /media/sda/wenjinma/envs/spatialATAC
import os
import pandas as pd

from method_utils import *
import utils

base_dir = "/media/sda/wenjinma"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore"

methods = ["GATVAE-SVGP-VAE", "spatialPCA", "GraphST", "SEDR", "spaGCN", "spaVAE", "STAligner", "PRECAST", "spatiAlign", "SPIRAL", "STADIA"]
annotation_key = "annotation"
seeds = [1993, 2024, 1105, 327, 113, 905, 789, 1314, 888, 666]
combinations = ["20_30", "30_30", "20_50", "30_50"]
simulation_combinations = [f"{seed}/{comb}" for seed in seeds for comb in combinations]
#simulation_combinations = ["666/10_30", "666/20_30", "666/30_30", "666/10_50", "666/20_50", "666/30_50"]
for simulation_comb in simulation_combinations:
    seed, DEG_settings = simulation_comb.split("/")
    DEG_x, DEG_y = DEG_settings.split("_")
    # create figure directory
    figure_dir = base_dir+os.sep+f"projects/spatialATAC-explore/summarize_benchmark/simulation_figures/{seed}-{DEG_x}_{DEG_y}"
    os.makedirs(figure_dir, exist_ok=True)
    spatial_dir = figure_dir+os.sep+"spatial_layout"
    os.makedirs(spatial_dir, exist_ok=True)
    latent_dir = figure_dir+os.sep+"latent_layout"
    os.makedirs(latent_dir, exist_ok=True)
    result_prefix = "simulation_results"+os.sep+simulation_comb

    spaDOT_result_prefix = result_prefix+os.sep+'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1'
    # --- split spaDOT results into two parts
    spaDOT_SVGP_adata = analyze_spaDOT_domains(result_dir, spaDOT_result_prefix, method="GATVAE-SVGP-VAE",
                                        label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"},
                                        latent_index=list(range(10)))
    utils.visualize_spatial(spaDOT_SVGP_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaDOT_SVGP", height=6)
    #utils.visualize_latent(spaDOT_SVGP_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_SVGPVAE")
    spaDOT_SVGP_adata.obs.to_csv(result_dir+os.sep+"explore_GATVAE-SVGP-VAE"+os.sep+spaDOT_result_prefix+os.sep+"SVGPVAE_predictions.csv")

    spaDOT_GAT_adata = analyze_spaDOT_domains(result_dir, spaDOT_result_prefix, method="GATVAE-SVGP-VAE",
                                        label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"},
                                        latent_index=list(range(11, 20)))
    utils.visualize_spatial(spaDOT_GAT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaDOT_GAT", height=6)
    #utils.visualize_latent(spaDOT_GAT_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_GATVAE")
    spaDOT_GAT_adata.obs.to_csv(result_dir+os.sep+"explore_GATVAE-SVGP-VAE"+os.sep+spaDOT_result_prefix+os.sep+"GATVAE_predictions.csv")

    # obtain results
    spaDOT_adata = analyze_spaDOT_domains(result_dir, spaDOT_result_prefix, method="GATVAE-SVGP-VAE",
                                    label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    spaDOT_adata.obs['pred_labels'] = spaDOT_adata.obs['pred_labels'].astype(str)
    spaDOT_adata.obs['spaDOT_pred_labels'] = spaDOT_adata.obs['timepoint'].astype(str) + '_' + spaDOT_adata.obs['pred_labels'].astype(str)
    custom_palette = {
        't0_0': '#FEFC53',
        't0_1': '#EA3423',
        't0_2': '#76FC4C',
        't1_0': '#FEFC53',
        't1_1': '#891AF5',
        't1_2': '#EA33F4',
        't1_3': '#76FC4C',
        't2_0': '#716A1E',
        't2_4': '#F5C242',
        't2_2': '#76FC4C',
        't2_1': '#891AF5',
        't2_3': '#EA33F4'
    }
    utils.visualize_spatial(spaDOT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="spaDOT_pred_labels", result_dir=spatial_dir, method="spaDOT", color_palette=custom_palette, height=4)

    spaDOT_adata.obs['true_labels'] = spaDOT_adata.obs.apply(utils.rename_GT_labels, axis=1)
    utils.visualize_spatial(spaDOT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="true_labels", result_dir=spatial_dir, method="GT", height=4, color_palette=utils.GT_custom_palette)
    #utils.visualize_latent(spaDOT_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT",
    #color_key=['timepoint', 'celltype', 'true_labels', 'spaDOT_pred_labels'])

    spatialPCA_adata = analyze_spatialPCA_domains(result_dir, result_prefix,
                           label_key='label', 
                           rename_dict={'label': "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
    utils.visualize_spatial(spatialPCA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spatialPCA")
    #utils.visualize_latent(spatialPCA_adata, use_rep_key='X', result_dir=latent_dir, method="spatialPCA")

    GraphST_adata = analyze_GraphST_domains(result_dir, result_prefix,
                                        label_key=annotation_key, 
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    utils.visualize_spatial(GraphST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="GraphST")
    #utils.visualize_latent(GraphST_adata, use_rep_key='emb', result_dir=latent_dir, method="GraphST")

    SEDR_adata = analyze_SEDR_domains(result_dir, result_prefix, label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    utils.visualize_spatial(SEDR_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="SEDR")
    #utils.visualize_latent(SEDR_adata, use_rep_key='SEDR.Harmony', result_dir=latent_dir, method="SEDR")

    spaGCN_adata = analyze_spaGCN_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    utils.visualize_spatial(spaGCN_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaGCN")
    #utils.visualize_latent(spaGCN_adata, use_rep_key='spaGCN_embed', result_dir=latent_dir, method="spaGCN")

    # --- integration methods
    spaVAE_adata = analyze_spaVAE_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    utils.visualize_spatial(spaVAE_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaVAE")
    #utils.visualize_latent(spaVAE_adata, use_rep_key='X', result_dir=latent_dir, method="spaVAE")

    STAligner_adata = analyze_STAligner_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    utils.visualize_spatial(STAligner_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="STAligner")
    #utils.visualize_latent(STAligner_adata, use_rep_key='STAligner', result_dir=latent_dir, method="STAligner")

    PRECAST_adata = analyze_PRECAST_domains(result_dir, result_prefix, label_key='label',
                            rename_dict={'label': "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
    utils.visualize_spatial(PRECAST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="PRECAST")
    #utils.visualize_latent(PRECAST_adata, use_rep_key='X', result_dir=latent_dir, method="PRECAST")

    spatiAlign_adata = analyze_spatiAlign_domains(result_dir, result_prefix)
    utils.visualize_spatial(spatiAlign_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spatiAlign")
    spatiAlign_adata.obs['spatiAlign_pred_labels'] = spatiAlign_adata.obs['timepoint'].astype(str) + '_' + spatiAlign_adata.obs['pred_labels'].astype(str)
    custom_palette = {
        't0_0': '#EA3423',
        't0_1': '#FEFC53',
        't0_2': '#76FC4C',
        't1_0': '#76FC4C',
        't1_1': '#EA3423',
        't1_2': '#FEFC53',
        't2_0': '#FEFC53',
        't2_2': '#EA3423',
        't2_3': '#76FC4C',
        't1_3': 'black',
        't2_1': 'black',
        't2_4': 'blue'
    }
    utils.visualize_spatial(spatiAlign_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="spatiAlign_pred_labels", result_dir=spatial_dir, method="spatiAlign", color_palette=custom_palette, height=4)
    spatiAlign_adata.obs['true_labels'] = spatiAlign_adata.obs.apply(utils.rename_GT_labels, axis=1)
    #utils.visualize_latent(spatiAlign_adata, use_rep_key='correct', result_dir=latent_dir, method="spatiAlign", color_key=['timepoint', 'celltype', 'true_labels', 'spatiAlign_pred_labels'])

    SPIRAL_adata = analyze_SPIRAL_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    utils.visualize_spatial(SPIRAL_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="SPIRAL")
    #utils.visualize_latent(SPIRAL_adata, use_rep_key='spiral', result_dir=latent_dir, method="SPIRAL")

    STADIA_adata = analyze_STADIA_domains(result_dir, result_prefix, label_key='label',
                        rename_dict={'label': "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
    utils.visualize_spatial(STADIA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="STADIA")
    #utils.visualize_latent(STADIA_adata, use_rep_key='X', result_dir=latent_dir, method="STADIA")
                            
    # --- combine results from all domains
    method_results_list = []
    for method in methods:
        if method == 'GATVAE-SVGP-VAE':
            method_result_dir = result_dir+os.sep+"explore_"+method+os.sep+spaDOT_result_prefix
        else:
            method_result_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
        method_results = pd.read_csv(method_result_dir+os.sep+method
        +"_domain_metrics.csv", index_col=0)
        method_results = method_results.reset_index()
        method_results_long = pd.melt(method_results, id_vars=['index'], var_name='timepoint', value_name='value')
        method_results_long['method'] = method
        method_results_list.append(method_results_long)
    final_method_results = pd.concat(method_results_list)
    final_method_results.to_csv(figure_dir+os.sep+"final_simulation_domain_metrics.csv")
