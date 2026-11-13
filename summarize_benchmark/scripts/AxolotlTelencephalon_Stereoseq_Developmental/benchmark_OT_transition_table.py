"""Score archived Axolotl OT benchmarks or recompute paper-setting spot OT.

Running without arguments retains the original archived comparison.
Use the compute subcommand with SpaDOT 1.0.7 to generate transition tables
from final 20-dimensional embeddings. For example, from the repository root:

    python summarize_benchmark/scripts/AxolotlTelencephalon_Stereoseq_Developmental/benchmark_OT_transition_table.py compute --latents INPUT --output OUTPUT

INPUT contains AxolotlTelencephalon/latent.h5ad. Use --dataset all to also
read ChickenHeart/latent.h5ad and MouseEmbryo/latent.h5ad. Each file uses
obs.timepoint (t0, t1, ...) and obs.annotation in the saved spot order.
For the threshold sweep, pass OUTPUT/pairs to score_thresholds.py.
For the manuscript panels, pass OUTPUT/inputs to render_transitions.R.
"""
import argparse
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from transition_reference_edges import edges


def accuracy(df, G):
    correct = 0
    wrong = 0
    for cell_i in df.index:
        for cell_j in df.columns:
            if (cell_i, cell_j) in list(G.edges):
                correct += df[cell_j][cell_i]
            else:
                wrong += df[cell_j][cell_i]
    return correct / (correct + wrong)

def evaluate_using_curated_transitions(df, G, cutoff=0.2):
    child_ct = list(df.columns)
    parent_ct = list(df.index)
    known_edges=list(G.edges)
    known_transition=0
    unknown_transition=0
    M=df.values
    for i,p in zip(range(len(parent_ct)), parent_ct):
        for j,c in zip(range(len(child_ct)),child_ct):
            edge_weight=M[i,j]
            if edge_weight>cutoff:
                if (p, c) in known_edges or c==p:
                    known_transition+=edge_weight
                else:
                    unknown_transition+=edge_weight
    all_transitions=unknown_transition+known_transition
    accuracy=known_transition/all_transitions
    import pandas as pd

    dfr=pd.DataFrame({'Accuracy':[accuracy], 'Total_weight':all_transitions})
    return(dfr)


def compare_archived_tables():
    import anndata
    import networkx as nx
    import numpy as np
    import pandas as pd

    result_dir = 'OT_transition_comparisons/'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    base_dir = '/net/mulan/home/wenjinma/projects/spatialATAC-explore'
    SpaDOT_output_dir = '/net/zootopia/disk1/wenjinma/minicluster/minicluster_backup_06102025/wenjinma/projects/spatialATAC-explore/explore_GATVAE-SVGP-VAE/Axolotl_results/'
    moscot_output_dir = f'{base_dir}/explore_moscot/Axolotl_results/'
    DeST_OT_output_dir = f'{base_dir}/explore_DeST-OT/Axolotl_results/'
    TOAST_output_dir = f'{base_dir}/explore_TOAST/Axolotl_results'

    latent_adata = anndata.read_h5ad(SpaDOT_output_dir+os.sep+'final_latent.h5ad')
    latent_adata.obs['day'] = latent_adata.obs['timepoint'].astype('category').cat.codes

    # --- load SpaDOT/moscot/DeST-OT transition tables and evaluate
    days = latent_adata.obs['day'].unique()
    days.sort()
    comparison_dfs = []
    for tp_i in range(len(days)-1):
        G = nx.DiGraph(edges("AxolotlTelencephalon", tp_i, ()))
        prev_day = days[tp_i]
        next_day = days[tp_i+1]
        # SpaDOT
        SpaDOT_transition_df = pd.read_csv(result_dir+os.sep+'SpaDOT_transition_prob_'+str(prev_day)+'_'+str(next_day)+'_min_norm.csv', index_col=0, header=0)
        # moscot -> need to specify GPU
        moscot_transition_row_norm_df = pd.read_csv(moscot_output_dir+os.sep+f'{prev_day}_{next_day}_row_normalized.csv', index_col=0, header=0)
        moscot_transition_col_norm_df = pd.read_csv(moscot_output_dir+os.sep+f'{prev_day}_{next_day}_column_normalized.csv', index_col=0, header=0)
        moscot_transition_min_df = np.minimum(moscot_transition_row_norm_df.values, moscot_transition_col_norm_df.values)
        moscot_transition_df = pd.DataFrame(moscot_transition_min_df, index=moscot_transition_row_norm_df.index, columns=moscot_transition_col_norm_df.columns)
        # DeST-OT
        DeST_OT_transition_row_norm_df = pd.read_csv(DeST_OT_output_dir+os.sep+f't{prev_day}_t{next_day}_row_normalized.csv', index_col=0, header=0)
        DeST_OT_transition_col_norm_df = pd.read_csv(DeST_OT_output_dir+os.sep+f't{prev_day}_t{next_day}_column_normalized.csv', index_col=0, header=0)
        DeST_OT_transition_min_df = np.minimum(DeST_OT_transition_row_norm_df.values, DeST_OT_transition_col_norm_df.values)
        DeST_OT_transition_df = pd.DataFrame(DeST_OT_transition_min_df, index=DeST_OT_transition_row_norm_df.index, columns=DeST_OT_transition_col_norm_df.columns)
        # TOAST
        TOAST_transition_row_norm_df = pd.read_csv(TOAST_output_dir+os.sep+f'{prev_day}_{next_day}_row_normalized.csv', index_col=0, header=0)
        TOAST_transition_col_norm_df = pd.read_csv(TOAST_output_dir+os.sep+f'{prev_day}_{next_day}_column_normalized.csv', index_col=0, header=0)
        TOAST_transition_min_df = np.minimum(TOAST_transition_row_norm_df.values, TOAST_transition_col_norm_df.values)
        TOAST_transition_df = pd.DataFrame(TOAST_transition_min_df, index=TOAST_transition_row_norm_df.index, columns=TOAST_transition_col_norm_df.columns)
        TOAST_transition_df.index.name = None
        # --- evalute accuracy
        print(f'Accuracy for transition from {prev_day} to {next_day}:')
        print('SpaDOT:', accuracy(SpaDOT_transition_df, G))
        print('Moscot:', accuracy(moscot_transition_df, G))
        print('DeST-OT:', accuracy(DeST_OT_transition_df, G))
        # --- evaluate based on cutoff
        SpaDOT_accuracy_df = evaluate_using_curated_transitions(SpaDOT_transition_df, G, cutoff=0.05)
        SpaDOT_accuracy_df['method'] = 'SpaDOT'
        moscot_accuracy_df = evaluate_using_curated_transitions(moscot_transition_df, G, cutoff=0.05)
        moscot_accuracy_df['method'] = 'Moscot'
        DeST_OT_accuracy_df = evaluate_using_curated_transitions(DeST_OT_transition_df, G, cutoff=0.05)
        DeST_OT_accuracy_df['method'] = 'DeST-OT'
        TOAST_accuracy_df = evaluate_using_curated_transitions(TOAST_transition_df, G, cutoff=0.05)
        TOAST_accuracy_df['method'] = 'TOAST'
        # --- results
        SpaDOT_accuracy_df['naive_accuracy'] = accuracy(SpaDOT_transition_df, G)
        moscot_accuracy_df['naive_accuracy'] = accuracy(moscot_transition_df, G)
        DeST_OT_accuracy_df['naive_accuracy'] = accuracy(DeST_OT_transition_df, G)
        TOAST_accuracy_df['naive_accuracy'] = accuracy(TOAST_transition_df, G)
        comparison_df = pd.concat([SpaDOT_accuracy_df, moscot_accuracy_df, DeST_OT_accuracy_df, TOAST_accuracy_df], axis=0)
        comparison_df['Day'] = f'{prev_day}->{next_day}'
        comparison_dfs.append(comparison_df)
    comparison_df = pd.concat(comparison_dfs, axis=0)
    comparison_df['Data'] = "Axolotl Telencephalon"
    comparison_df.to_csv(result_dir+os.sep+'OT_method_comparison.csv')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command")
    compute = commands.add_parser("compute", help="Compute paper-setting spot and domain OT tables.")
    compute.add_argument("--latents", type=Path, required=True, help="Root containing DATASET/latent.h5ad.")
    compute.add_argument("--output", type=Path, required=True, help="New output directory.")
    compute.add_argument("--dataset", choices=("AxolotlTelencephalon", "MouseEmbryo", "ChickenHeart", "all"),
                         default="AxolotlTelencephalon")
    compute.add_argument("--domain-keys", type=Path, help="Root containing DATASET/domain_key.csv for plotting.")
    compute.add_argument("--save-spot-plans", action="store_true")
    args = parser.parse_args()
    if args.command == "compute":
        from OT_loss.spot_transition_tables import compute_tables

        compute_tables(args)
    else:
        compare_archived_tables()


if __name__ == "__main__":
    main()
