# conda activate /media/sda/wenjinma/envs/spatialATAC
import os
import pandas as pd

from method_utils import *
import utils

base_dir = "/media/sda/wenjinma"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore"
figure_dir = result_dir+os.sep+"/summarize_benchmark/MOSTA_brain_figures/spatial_layout"
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

methods = ["GATVAE-SVGP-VAE", "spatialPCA", "GraphST", "SEDR", "spaGCN", "spaVAE", "STAligner", "PRECAST", "spatiAlign", "SPIRAL", "STADIA"]
annotation_key = "annotation"
result_prefix = "MOSTA_brain_results"

# --- analyze results
for method in methods:
    method_result_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
    if method == 'GATVAE-SVGP-VAE': method = 'spaDOT'
    method_pred_df = pd.read_csv(method_result_dir+os.sep+method+"_pred.csv", index_col=0)
    res_dict = dict()
    for tp in set(method_pred_df["timepoint"]):
        tp_res = method_pred_df[method_pred_df["timepoint"] == tp]
        loc = np.array(tp_res[["pixel_x", "pixel_y"]])
        pred = np.array(tp_res["pred_labels"])
        tp_ASW = utils._compute_ASW(loc, pred)
        tp_PAS = utils._compute_PAS(loc, pred)
        tp_CHAOS = utils._compute_CHAOS(loc, pred)
        res_dict[tp] = {}
        res_dict[tp]["ASW"], res_dict[tp]["PAS"], res_dict[tp]["CHAOS"] = tp_ASW, tp_PAS, tp_CHAOS
    res_df = pd.DataFrame.from_dict(res_dict)
    res_df['avg'] = res_df.mean(axis=1)
    res_df.to_csv(method_result_dir+os.sep+method+'_domain_metrics.csv')

method_results_list = []
for method in methods:
    method_result_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
    if method == 'GATVAE-SVGP-VAE': 
        method_results = pd.read_csv(method_result_dir+os.sep+'spaDOT_domain_metrics.csv', index_col=0)
    else:
        method_results = pd.read_csv(method_result_dir+os.sep+method+'_domain_metrics.csv', index_col=0)
    method_results = method_results.reset_index()
    method_results_long = pd.melt(method_results, id_vars=['index'], var_name='timepoint', value_name='value')
    method_results_long['method'] = method
    method_results_list.append(method_results_long)
final_method_results = pd.concat(method_results_list)
final_method_results.to_csv(os.path.dirname(figure_dir)+os.sep+"final_MOSTA_brain_domain_metrics.csv")