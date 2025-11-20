# conda activate /media/sda/wenjinma/envs/spatialATAC
import os
import pandas as pd

import sys
sys.path.append('/net/mulan/home/wenjinma/projects/SpaDOT_analysis/summarize_benchmark/scripts')

from method_utils import *
import utils

base_dir = "/net/mulan/home/wenjinma"
result_dir = base_dir+os.sep+"projects/SpaDOT_analysis"
figure_dir = base_dir+os.sep+"projects/SpaDOT_analysis/summarize_benchmark/figures/ChickenHeart_figures/spatial_layout"
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)
latent_dir = base_dir+os.sep+"projects/SpaDOT_analysis/summarize_benchmark/figures/ChickenHeart_figures/latent_layout"
if not os.path.exists(latent_dir):
    os.makedirs(latent_dir)

methods = ["SpaDOT", "spatialPCA", "GraphST", "SEDR", "spaGCN", "spaVAE", "STAligner", "PRECAST", "spatiAlign", "SPIRAL", "STADIA", 'STAMP']
annotation_key = "annotation"

# --- analyze data complexity
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_ChickenHeart"
result_prefix = "ChickenHeart_results"

# --- split spaDOT results into two parts
# spaDOT_SVGP_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="GATVAE-SVGP-VAE",
#                                     label_key=annotation_key,
#                                     rename_dict={annotation_key: "true_labels"},
#                                     latent_index=list(range(10)))
# utils.visualize_spatial(spaDOT_SVGP_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaDOT_SVGP", height=6)
# utils.visualize_latent(spaDOT_SVGP_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_SVGPVAE")
# spaDOT_SVGP_adata.obs.to_csv(result_dir+os.sep+"explore_GATVAE-SVGP-VAE"+os.sep+result_prefix+os.sep+"SVGPVAE_predictions.csv")

# spaDOT_GAT_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="GATVAE-SVGP-VAE",
#                                     label_key=annotation_key,
#                                     rename_dict={annotation_key: "true_labels"},
#                                     latent_index=list(range(11, 20)))
# utils.visualize_spatial(spaDOT_GAT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaDOT_GAT", height=6)
# utils.visualize_latent(spaDOT_GAT_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_GATVAE")
# spaDOT_GAT_adata.obs.to_csv(result_dir+os.sep+"explore_GATVAE-SVGP-VAE"+os.sep+result_prefix+os.sep+"GATVAE_predictions.csv")

# --- start analyzing results
spaDOT_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="SpaDOT",
                                    label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels"})
# utils.visualize_spatial(spaDOT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="GATVAE-SVGP-VAE")
# # visualize ground truth
# utils.visualize_spatial(spaDOT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="true_labels", result_dir=figure_dir, method="GT", height=6)
# utils.visualize_latent(spaDOT_adata, use_rep_key='X', result_dir=latent_dir, method="GATVAE-SVGP-VAE")



# --- clustering with multi-slice methods
spatialPCA_adata = analyze_spatialPCA_domains(result_dir, result_prefix,
                           label_key=annotation_key, 
                           rename_dict={annotation_key: "true_labels"})
# utils.visualize_spatial(spatialPCA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spatialPCA")
# utils.visualize_latent(spatialPCA_adata, use_rep_key='X', result_dir=latent_dir, method="spatialPCA")


GraphST_adata = analyze_GraphST_domains(result_dir, result_prefix,
                                      label_key=annotation_key, 
                                      rename_dict={annotation_key: "true_labels"})
# utils.visualize_spatial(GraphST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="GraphST")
# utils.visualize_latent(GraphST_adata, use_rep_key='emb', result_dir=latent_dir, method="GraphST")

SEDR_adata = analyze_SEDR_domains(result_dir, result_prefix, label_key=annotation_key,
                                rename_dict={annotation_key: "true_labels"})
# utils.visualize_spatial(SEDR_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="SEDR")
# utils.visualize_latent(SEDR_adata, use_rep_key='SEDR.Harmony', result_dir=latent_dir, method="SEDR")

spaGCN_adata = analyze_spaGCN_domains(result_dir, result_prefix, label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels"})
# utils.visualize_spatial(spaGCN_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaGCN")
# utils.visualize_latent(spaGCN_adata, use_rep_key='spaGCN_embed', result_dir=latent_dir, method="spaGCN")

# --- integration methods
spaVAE_adata = analyze_spaVAE_domains(result_dir, result_prefix, label_key=annotation_key)
# utils.visualize_spatial(spaVAE_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaVAE")
# utils.visualize_latent(spaVAE_adata, use_rep_key='X', result_dir=latent_dir, method="spaVAE")

STAligner_adata = analyze_STAligner_domains(result_dir, result_prefix, label_key=annotation_key)
# utils.visualize_spatial(STAligner_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="STAligner")
# utils.visualize_latent(STAligner_adata, use_rep_key='STAligner', result_dir=latent_dir, method="STAligner")

PRECAST_adata = analyze_PRECAST_domains(result_dir, result_prefix, label_key=annotation_key, rename_dict={annotation_key: "true_labels"})
# utils.visualize_spatial(PRECAST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="PRECAST")
# utils.visualize_latent(PRECAST_adata, use_rep_key='X', result_dir=latent_dir, method="PRECAST")

spatiAlign_adata = analyze_spatiAlign_domains(result_dir, result_prefix)
# utils.visualize_spatial(spatiAlign_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spatiAlign")
# utils.visualize_latent(spatiAlign_adata, use_rep_key='correct', result_dir=latent_dir, method="spatiAlign")

SPIRAL_adata = analyze_SPIRAL_domains(result_dir, result_prefix, label_key=annotation_key)
# utils.visualize_spatial(SPIRAL_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="SPIRAL")
# utils.visualize_latent(SPIRAL_adata, use_rep_key='spiral', result_dir=latent_dir, method="SPIRAL")

STADIA_adata = analyze_STADIA_domains(result_dir, result_prefix, label_key=annotation_key, rename_dict={annotation_key: "true_labels"})
# utils.visualize_spatial(STADIA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="STADIA")
# utils.visualize_latent(STADIA_adata, use_rep_key='X', result_dir=latent_dir, method="STADIA")

STAMP_adata = analyze_STAMP_domains(result_dir, result_prefix, label_key=annotation_key,
                       rename_dict={annotation_key: "true_labels"})


# --- combine results from all domains
# method_results_list = []
# for method in methods:
#     method_result_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
#     method_results = pd.read_csv(method_result_dir+os.sep+method
#     +"_kmeans_domain_metrics.csv", index_col=0)
#     method_results = method_results.reset_index()
#     method_results_long = pd.melt(method_results, id_vars=['index'], var_name='timepoint', value_name='value')
#     method_results_long['method'] = method
#     method_results_list.append(method_results_long)
# final_method_results = pd.concat(method_results_list)
# final_method_results.to_csv(os.path.dirname(figure_dir)+os.sep+"final_ChickenHeart_domain_metrics.csv")