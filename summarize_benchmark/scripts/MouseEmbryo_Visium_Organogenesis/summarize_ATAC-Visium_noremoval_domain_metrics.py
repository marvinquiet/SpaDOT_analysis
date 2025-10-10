# conda activate /media/sda/wenjinma/envs/spatialATAC
import os
import pandas as pd

from method_utils import *
import utils

base_dir = "/net/mulan/home/wenjinma"
result_dir = base_dir+os.sep+"projects/SpaDOT_analysis"
figure_dir = base_dir+os.sep+"projects/SpaDOT_analysis/summarize_benchmark/figures/ATAC-Visium_noremoval_figures/spatial_layout"
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

latent_dir = base_dir+os.sep+"projects/SpaDOT_analysis/summarize_benchmark/figures/ATAC-Visium_noremoval_figures/latent_layout"
if not os.path.exists(latent_dir):
    os.makedirs(latent_dir)

methods = ["SpaDOT", "spatialPCA", "GraphST", "SEDR", "spaGCN", "spaVAE", "STAligner", "PRECAST", "spatiAlign", "SPIRAL", "STADIA", "STAMP"]
annotation_key = "annotation"

# --- analyze data complexity
data_dir = "/net/zootopia/disk1/wenjinma/data/spatial_scRNAseq_GSE214989"
# analyze_ATACVisium_complexity(data_dir)
result_prefix = "ATAC-Visium_noremoval_results"
# spaDOT_res = analyze_spaDOT_domains(result_dir, result_prefix, method="SpaDOT",
#                                     label_key=annotation_key,
#                                     rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
# utils.visualize_ATACVisium_spatial(spaDOT_res, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaDOT")
# visualize ground truth
# spaDOT_res['seurat_clusters_harmony'] = spaDOT_res['seurat_clusters_harmony'].astype('category')
# utils.visualize_ATACVisium_spatial(spaDOT_res, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="seurat_clusters_harmony", result_dir=figure_dir, method="GT_refined")
# utils.visualize_ATACVisium_spatial(spaDOT_res, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="true_labels", result_dir=figure_dir, method="GT_combined")

# --- split spaDOT results into two parts
# spaDOT_SVGP_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="GATVAE-SVGP-VAE",
#                                     label_key=annotation_key,
#                                     rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"},
#                                     latent_index=list(range(10)))
# utils.visualize_ATACVisium_spatial(spaDOT_SVGP_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaDOT_SVGP", height=6)
# utils.visualize_latent(spaDOT_SVGP_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_SVGPVAE")
# spaDOT_SVGP_adata.obs.to_csv(result_dir+os.sep+"explore_GATVAE-SVGP-VAE"+os.sep+result_prefix+os.sep+"SVGPVAE_predictions.csv")

# spaDOT_GAT_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="GATVAE-SVGP-VAE",
#                                     label_key=annotation_key,
#                                     rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"},
#                                     latent_index=list(range(11, 20)))
# utils.visualize_ATACVisium_spatial(spaDOT_GAT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaDOT_GAT", height=6)
# # utils.visualize_latent(spaDOT_GAT_adata, use_rep_key='X', result_dir=latent_dir, method="spaDOT_GATVAE")
# spaDOT_GAT_adata.obs.to_csv(result_dir+os.sep+"explore_GATVAE-SVGP-VAE"+os.sep+result_prefix+os.sep+"GATVAE_predictions.csv")

# --- start analyzing results
spaDOT_adata = analyze_spaDOT_domains(result_dir, result_prefix, method="SpaDOT",
                                    label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
# utils.visualize_ATACVisium_spatial(spaDOT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method=method)
# utils.visualize_latent(spaDOT_adata, use_rep_key='X', result_dir=latent_dir, method=method)
# visualize ground truth
# utils.visualize_ATACVisium_spatial(spaDOT_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="true_labels", result_dir=figure_dir, method="GT", height=6)


# # --- calculate spatial continuity metrics for GT
# timepoints = list(spaDOT_res["timepoint"].unique())
# timepoints.sort()
# res_dict = dict()
# for tp in timepoints:
#     tp_res = spaDOT_res[spaDOT_res["timepoint"]==tp]
#     # --- domain continuity
#     loc = np.array(tp_res[["pixel_x", "pixel_y"]])
#     pred = np.array(tp_res["true_labels"])
#     tp_ASW = utils._compute_ASW(loc, pred)
#     tp_PAS = utils._compute_PAS(loc, pred)
#     tp_CHAOS = utils._compute_CHAOS(loc, pred)
#     #print(f"{method} at timepoint {tp}: ASW={tp_ASW}, PAS={tp_PAS}, CHAOS={tp_CHAOS}")
#     # save results
#     res_dict[tp] = {}
#     res_dict[tp]["ASW"], res_dict[tp]["PAS"], res_dict[tp]["CHAOS"] = tp_ASW, tp_PAS, tp_CHAOS
# res_df = pd.DataFrame(res_dict)
# res_df['avg'] = res_df.mean(axis=1)
# print(res_df)
# res_df.to_csv(result_dir+os.sep+'GT_domain_metrics.csv')
# # --- these metrics based on GT label is not performed well.
 
# --- clustering with multi-slice methods
spatialPCA_adata = analyze_spatialPCA_domains(result_dir, result_prefix,
                           label_key=annotation_key, 
                           rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
# utils.visualize_ATACVisium_spatial(spatialPCA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spatialPCA")
# utils.visualize_latent(spatialPCA_adata, use_rep_key='X', result_dir=latent_dir, method="spatialPCA")

GraphST_adata = analyze_GraphST_domains(result_dir, result_prefix,
                                      label_key=annotation_key, 
                                      rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
# utils.visualize_ATACVisium_spatial(GraphST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="GraphST")
# utils.visualize_latent(GraphST_adata, use_rep_key='emb', result_dir=latent_dir, method="GraphST")

SEDR_adata = analyze_SEDR_domains(result_dir, result_prefix, label_key=annotation_key,
                                rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
utils.visualize_ATACVisium_spatial(SEDR_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="SEDR")
# utils.visualize_latent(SEDR_adata, use_rep_key='SEDR.Harmony', result_dir=latent_dir, method="SEDR")

spaGCN_adata = analyze_spaGCN_domains(result_dir, result_prefix, label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
utils.visualize_ATACVisium_spatial(spaGCN_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaGCN")
# utils.visualize_latent(spaGCN_adata, use_rep_key='spaGCN_embed', result_dir=latent_dir, method="spaGCN")

# --- integration methods
spaVAE_adata = analyze_spaVAE_domains(result_dir, result_prefix, label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
utils.visualize_ATACVisium_spatial(spaVAE_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spaVAE")
# utils.visualize_latent(spaVAE_adata, use_rep_key='X', result_dir=latent_dir, method="spaVAE")

STAligner_adata = analyze_STAligner_domains(result_dir, result_prefix, label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
utils.visualize_ATACVisium_spatial(STAligner_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="STAligner")
# utils.visualize_latent(STAligner_adata, use_rep_key='STAligner', result_dir=latent_dir, method="STAligner")

PRECAST_adata = analyze_PRECAST_domains(result_dir, result_prefix, label_key=annotation_key,
                        rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
utils.visualize_ATACVisium_spatial(PRECAST_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="PRECAST")
# utils.visualize_latent(PRECAST_adata, use_rep_key='X', result_dir=latent_dir, method="PRECAST")

spatiAlign_adata = analyze_spatiAlign_domains(result_dir, result_prefix)
utils.visualize_ATACVisium_spatial(spatiAlign_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="spatiAlign")
# utils.visualize_latent(spatiAlign_adata, use_rep_key='correct', result_dir=latent_dir, method="spatiAlign")

SPIRAL_adata = analyze_SPIRAL_domains(result_dir, result_prefix, label_key=annotation_key,
                                    rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})
utils.visualize_ATACVisium_spatial(SPIRAL_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="SPIRAL")
# utils.visualize_latent(SPIRAL_adata, use_rep_key='spiral', result_dir=latent_dir, method="SPIRAL")

STADIA_adata = analyze_STADIA_domains(result_dir, result_prefix, label_key=annotation_key,
                       rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"})
utils.visualize_ATACVisium_spatial(STADIA_adata.obs, coord_x_key="pixel_x", coord_y_key="pixel_y", pred_label_key="pred_labels", result_dir=figure_dir, method="STADIA")
# utils.visualize_latent(STADIA_adata, use_rep_key='X', result_dir=latent_dir, method="STADIA")

STAMP_adata = analyze_STAMP_domains(result_dir, result_prefix, label_key=annotation_key,
                       rename_dict={annotation_key: "true_labels", "kmeans": "pred_labels"})



# --- combine results from all domains
method_results_list = []
for method in methods:
    method_result_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
    method_results = pd.read_csv(method_result_dir+os.sep+method
    +"_domain_metrics.csv", index_col=0)
    method_results = method_results.reset_index()
    method_results_long = pd.melt(method_results, id_vars=['index'], var_name='timepoint', value_name='value')
    method_results_long['method'] = method
    method_results_list.append(method_results_long)
final_method_results = pd.concat(method_results_list)
final_method_results.to_csv(os.path.dirname(figure_dir)+os.sep+"final_ATAC-Visium_noremoval_domain_metrics.csv")