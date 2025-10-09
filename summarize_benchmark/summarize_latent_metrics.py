# conda activate /media/sda/wenjinma/envs/spatialATAC
import os, argparse
from method_utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', type=str, help='Dataset name')
args = parser.parse_args()
dataset = args.dataset

result_dir = "/media/sda/wenjinma/projects/spatialATAC-explore"
figure_dir = result_dir+os.sep+"summarize_benchmark"+os.sep+dataset+"_figures/latent_layout"
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)

def analyze_latent(result_dir, result_prefix):
    # --- analyze latent metrics
    analyze_spaDOT_latent(result_dir, result_prefix, method="GATVAE-SVGP-VAE")

    # --- clustering with multi-slice methods
    analyze_spatialPCA_latent(result_dir, result_prefix)
    analyze_GraphST_latent(result_dir, result_prefix)
    analyze_SEDR_latent(result_dir, result_prefix)
    analyze_spaGCN_latent(result_dir, result_prefix)
    # --- integration methods
    analyze_spaVAE_latent(result_dir, result_prefix)
    analyze_STAligner_latent(result_dir, result_prefix)
    analyze_PRECAST_latent(result_dir, result_prefix)
    analyze_spatiAlign_latent(result_dir, result_prefix)
    analyze_SPIRAL_latent(result_dir, result_prefix)
    analyze_STADIA_latent(result_dir, result_prefix)


if dataset == 'ATAC-Visium' or dataset == 'ATAC-Visium_noremoval' or dataset == 'Axolotl' or dataset == 'ChickenHeart':
    result_prefix = dataset + "_results"
    analyze_latent(result_dir, result_prefix)
elif dataset == 'simulation':
    simulation_combinations = ["666/10_30", "666/20_30", "666/30_30", "666/10_50", "666/20_50", "666/30_50"]
    for comb in simulation_combinations:
        result_prefix = "simulation_results" + os.sep + comb
        analyze_latent(result_dir, result_prefix)


