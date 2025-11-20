# conda activate ~/envs/spatialATAC
# on mulan
import os
import argparse
import sys
import pandas as pd

from simulation_prop_method_utils import *
import simulation_prop_utils

result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis"

parser = argparse.ArgumentParser(description="Summarize results from simulation proportion scenarios")
parser.add_argument('--DEG', type=int, help='Simulation number of DEGs')
parser.add_argument('--mixedNPCs', action='store_true', help='Whether the NPCs is a mixture')
args = parser.parse_args()

methods = ["SpaDOT", "spatialPCA", "GraphST", "SEDR", "spaGCN", "spaVAE", "STAligner", "PRECAST", "spatiAlign", "SPIRAL", "STADIA", "STAMP"]
annotation_key = "annotation"
seeds = [1993, 2024, 1105]
combinations = ["prop_0.9", "prop_0.8", "prop_0.7", "prop_0.6", "prop_0.5"]
simulation_combinations = [f"{seed}/{comb}" for seed in seeds for comb in combinations]
for simulation_comb in simulation_combinations:
    print(simulation_comb)
    seed, prop = simulation_comb.split("/")
    # create figure directory
    if args.mixedNPCs:
        figure_dir = result_dir+os.sep+f"summarize_benchmark/figures/simulation_prop_figures_mixedNPCs_DEG{args.DEG}/{seed}/{prop}"
        result_prefix = f"simulation_prop_results_mixedNPCs_DEG{args.DEG}"+os.sep+simulation_comb
    else:
        figure_dir = result_dir+os.sep+f"summarize_benchmark/figures/simulation_prop_figures_DEG{args.DEG}/{seed}/{prop}"
        result_prefix = f"simulation_prop_results_DEG{args.DEG}"+os.sep+simulation_comb
    os.makedirs(figure_dir, exist_ok=True)
    spatial_dir = figure_dir+os.sep+"spatial_layout"
    os.makedirs(spatial_dir, exist_ok=True)
    latent_dir = figure_dir+os.sep+"latent_layout"
    os.makedirs(latent_dir, exist_ok=True)

    spaDOT_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="SpaDOT",
                                label_key=annotation_key,
                                rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    # --- split spaDOT results into two parts
    spaDOT_SVGP_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="SpaDOT",
                                        label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"},
                                        latent_index=list(range(10)))
    simulation_prop_utils.visualize_spatial(spaDOT_SVGP_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaDOT_SVGP", height=6)
    #utils.visualize_latent(spaDOT_SVGP_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_SVGPVAE")
    spaDOT_SVGP_adata.obs.to_csv(result_dir+os.sep+"explore_SpaDOT"+os.sep+result_prefix+os.sep+"SVGPVAE_predictions.csv")

    spaDOT_GAT_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="SpaDOT",
                                        label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"},
                                        latent_index=list(range(11, 20)))
    simulation_prop_utils.visualize_spatial(spaDOT_GAT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaDOT_GAT", height=6)
    #utils.visualize_latent(spaDOT_GAT_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_GATVAE")
    spaDOT_GAT_adata.obs.to_csv(result_dir+os.sep+"explore_SpaDOT"+os.sep+result_prefix+os.sep+"GATVAE_predictions.csv")

    # obtain results
    spaDOT_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="SpaDOT",
                                    label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    spaDOT_adata.obs['pred_labels'] = spaDOT_adata.obs['pred_labels'].astype(str)
    spaDOT_adata.obs['spaDOT_pred_labels'] = spaDOT_adata.obs['timepoint'].astype(str) + '_' + spaDOT_adata.obs['pred_labels'].astype(str)
    spaDOT_adata.obs['true_labels'] = spaDOT_adata.obs.apply(simulation_prop_utils.rename_GT_labels, axis=1)
    simulation_prop_utils.visualize_spatial(spaDOT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", 
        pred_label_key="true_labels", result_dir=spatial_dir, method="GT", height=4, 
        color_palette=simulation_prop_utils.GT_custom_palette)
    #utils.visualize_latent(spaDOT_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT",
    #color_key=['timepoint', 'celltype', 'true_labels', 'spaDOT_pred_labels'])

    spatialPCA_adata = analyze_spatialPCA_domains(result_dir, result_prefix,
                           label_key='label', 
                           rename_dict={'label': "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
    simulation_prop_utils.visualize_spatial(spatialPCA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spatialPCA")
    #utils.visualize_latent(spatialPCA_adata, use_rep_key='X', result_dir=latent_dir, method="spatialPCA")

    GraphST_adata = analyze_GraphST_domains(result_dir, result_prefix,
                                        label_key=annotation_key, 
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    simulation_prop_utils.visualize_spatial(GraphST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="GraphST")
    #utils.visualize_latent(GraphST_adata, use_rep_key='emb', result_dir=latent_dir, method="GraphST")

    SEDR_adata = analyze_SEDR_domains(result_dir, result_prefix, label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    simulation_prop_utils.visualize_spatial(SEDR_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="SEDR")
    #utils.visualize_latent(SEDR_adata, use_rep_key='SEDR.Harmony', result_dir=latent_dir, method="SEDR")

    spaGCN_adata = analyze_spaGCN_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    simulation_prop_utils.visualize_spatial(spaGCN_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaGCN")
    #utils.visualize_latent(spaGCN_adata, use_rep_key='spaGCN_embed', result_dir=latent_dir, method="spaGCN")

    # --- integration methods
    spaVAE_adata = analyze_spaVAE_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    simulation_prop_utils.visualize_spatial(spaVAE_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spaVAE")
    #utils.visualize_latent(spaVAE_adata, use_rep_key='X', result_dir=latent_dir, method="spaVAE")

    STAligner_adata = analyze_STAligner_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    simulation_prop_utils.visualize_spatial(STAligner_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="STAligner")
    #utils.visualize_latent(STAligner_adata, use_rep_key='STAligner', result_dir=latent_dir, method="STAligner")

    PRECAST_adata = analyze_PRECAST_domains(result_dir, result_prefix, label_key='label',
                            rename_dict={'label': "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
    simulation_prop_utils.visualize_spatial(PRECAST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="PRECAST")
    #utils.visualize_latent(PRECAST_adata, use_rep_key='X', result_dir=latent_dir, method="PRECAST")

    spatiAlign_adata = analyze_spatiAlign_domains(result_dir, result_prefix)
    simulation_prop_utils.visualize_spatial(spatiAlign_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="spatiAlign")
    # utils.visualize_latent(spatiAlign_adata, use_rep_key='correct', result_dir=latent_dir, method="spatiAlign", color_key=['timepoint', 'celltype', 'true_labels', 'spatiAlign_pred_labels'])

    SPIRAL_adata = analyze_SPIRAL_domains(result_dir, result_prefix, label_key=annotation_key,
                                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    simulation_prop_utils.visualize_spatial(SPIRAL_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="SPIRAL")
    # utils.visualize_latent(SPIRAL_adata, use_rep_key='spiral', result_dir=latent_dir, method="SPIRAL")

    STADIA_adata = analyze_STADIA_domains(result_dir, result_prefix, label_key='label',
                        rename_dict={'label': "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
    simulation_prop_utils.visualize_spatial(STADIA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="STADIA")
    #utils.visualize_latent(STADIA_adata, use_rep_key='X', result_dir=latent_dir, method="STADIA")

    STAMP_adata = analyze_STAMP_domains(result_dir, result_prefix, label_key=annotation_key,
                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
    simulation_prop_utils.visualize_spatial(STAMP_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=spatial_dir, method="STAMP")
    #utils.visualize_latent(STAMP_adata, use_rep_key='X', result_dir=latent_dir, method="STAMP")
                            
    # --- combine results from all domains
    method_results_list = []
    for method in methods:
        method_result_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
        method_results = pd.read_csv(method_result_dir+os.sep+method+"_domain_metrics.csv", index_col=0)
        method_results = method_results.reset_index()
        method_results_long = pd.melt(method_results, id_vars=['index'], var_name='timepoint', value_name='value')
        method_results_long['method'] = method
        method_results_list.append(method_results_long)
    final_method_results = pd.concat(method_results_list)
    final_method_results.to_csv(figure_dir+os.sep+"final_simulation_domain_metrics.csv")
