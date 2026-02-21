# conda activate ~/envs/spatialATAC
import os
import numpy as np
import pandas as pd
import wot
import anndata
import networkx as nx

result_dir = 'OT_transition_comparisons/'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
base_dir = '/net/mulan/home/wenjinma/projects/spatialATAC-explore'
SpaDOT_output_dir = '/net/zootopia/disk1/wenjinma/minicluster/minicluster_backup_06102025/wenjinma/projects/spatialATAC-explore/explore_GATVAE-SVGP-VAE/ChickenHeart_results/'
moscot_output_dir = f'{base_dir}/explore_moscot/ChickenHeart_results/'
DeST_OT_output_dir = f'{base_dir}/explore_DeST-OT/ChickenHeart_results/'
TOAST_output_dir = f'{base_dir}/explore_TOAST/ChickenHeart_results'

# SpaDOT normalization based on latent and annotation
latent_adata = anndata.read_h5ad(SpaDOT_output_dir+os.sep+'final_latent.h5ad')
latent_adata.obs['day'] = latent_adata.obs['timepoint'].astype('category').cat.codes
latent_adata.obs['cell_growth_rate'] = 1
# ot_model = wot.ot.OTModel(latent_adata, epsilon = 0.05, epsilon0 = 1, lambda1 = 0.1,lambda2 = 5, growth_iters=3) # imbalanced OT with imbalanced row and comparatiely balanced target
# ot_model.compute_all_transport_maps(tmap_out=result_dir+os.sep+'OT') # compute transport maps
# tmap_model = wot.tmap.TransportMapModel.from_directory(result_dir+os.sep+'OT')
# latent_cell_sets = latent_adata.obs.groupby('annotation').apply(lambda x: x.index.tolist()).to_dict()
# days = latent_adata.obs['day'].unique()
# days.sort()
# for tp_i in range(len(days)-1):  
#     prev_day = days[tp_i]
#     next_day = days[tp_i+1]
#     prev_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=prev_day)
#     next_day_populations = tmap_model.population_from_cell_sets(latent_cell_sets, at_time=next_day)
#     transition_table = tmap_model.transition_table(prev_day_populations, next_day_populations) # aggregated OT matrix
#     transition_df = pd.DataFrame(transition_table.X, index=transition_table.obs_names, columns=transition_table.var_names)
#     # --- normalize by column sum
#     transition_prob = transition_table.X
#     transition_prob = transition_prob/transition_prob.sum(axis=0, keepdims=True)
#     transition_prob_col_norm_df = pd.DataFrame(transition_prob, index=transition_table.obs_names, columns=transition_table.var_names)
#     # --- normalize by row sum
#     transition_prob = transition_table.X
#     transition_prob = transition_prob/transition_prob.sum(axis=1, keepdims=True)
#     transition_prob_row_norm_df = pd.DataFrame(transition_prob, index=transition_table.obs_names, columns=transition_table.var_names)
#     # --- obtain minimum of both normalized transition probabilities
#     transition_prob_min_df = np.minimum(transition_prob_col_norm_df.values, transition_prob_row_norm_df.values)
#     SpaDOT_transition_prob_min_df = pd.DataFrame(transition_prob_min_df, index=transition_table.obs_names, columns=transition_table.var_names)
#     SpaDOT_transition_prob_min_df.to_csv(result_dir+os.sep+'SpaDOT_transition_prob_'+str(prev_day)+'_'+str(next_day)+'_min_norm.csv')

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
    dfr=pd.DataFrame({'Accuracy':[accuracy], 'Total_weight':all_transitions})
    return(dfr)


# --- in this scenario, there is only one transition according to the annotation
regions = latent_adata.obs['annotation'].unique().tolist()
t0_t1_G = nx.DiGraph()
t0_t1_G.add_nodes_from(['Epicardium- like', 'Outflow tract', 'Ventricle', 'Valves', 'Atria', 'Epicardium', 'Endothelium', 
                    'Trabecular LV and \nendocardium', 'Right ventricle', 'Compact LV and \ninter-ventricular septum'])
t0_t1_G.add_edges_from( # suggested by chatGPT based on name match
    [('Epicardium- like', 'Epicardium'),
    ('Ventricle', 'Trabecular LV and \nendocardium'),
    ('Ventricle', 'Right ventricle'),
    ('Ventricle', 'Compact LV and \ninter-ventricular septum'),
    ('Atria', 'Atria'),
    #('Outflow tract', 'Valves'),
    ('Valves', 'Valves')
])
print(t0_t1_G.edges)

t1_t2_G = nx.DiGraph()
t1_t2_G.add_nodes_from(['Valves', 'Atria', 'Epicardium', 'Endothelium', 
                    'Trabecular LV and \nendocardium', 'Right ventricle', 'Compact LV and \ninter-ventricular septum',
                    'Outflow tract'])
t1_t2_G.add_edges_from( # suggested by chatGPT
    [('Valves', 'Valves'),
    ('Atria', 'Atria'),
    ('Epicardium', 'Epicardium'),
    ('Trabecular LV and \nendocardium', 'Trabecular LV and \nendocardium'),
    ('Right ventricle', 'Right ventricle'),
    ('Compact LV and \ninter-ventricular septum', 'Compact LV and \ninter-ventricular septum')
])
print(t1_t2_G.edges)

t2_t3_G = nx.DiGraph()
t2_t3_G.add_nodes_from(['Epicardium', 'Endocardium', 'Right ventricle', 'Compact LV and \ninter-ventricular septum',
                    'Atria', 'Valves'])
t2_t3_G.add_edges_from( # suggested by chatGPT
    [('Epicardium', 'Epicardium'),
    ('Right ventricle', 'Right ventricle'),
    ('Trabecular LV and \nendocardium', 'Endocardium'), 
    ('Compact LV and \ninter-ventricular septum', 'Compact LV and \ninter-ventricular septum'),
    ('Atria', 'Atria'),
    ('Valves', 'Valves')
])
print(t2_t3_G.edges)


# --- load SpaDOT/moscot/DeST-OT transition tables and evaluate
days = latent_adata.obs['day'].unique()
days.sort()
GT_Gs = [t0_t1_G, t1_t2_G, t2_t3_G]
comparison_dfs = []
for tp_i in range(len(days)-1):
    G = GT_Gs[tp_i]
    prev_day = days[tp_i]
    next_day = days[tp_i+1]
    # SpaDOT
    SpaDOT_transition_df = pd.read_csv(result_dir+os.sep+'SpaDOT_transition_prob_'+str(prev_day)+'_'+str(next_day)+'_min_norm.csv', index_col=0, header=0)
    # # moscot -> cannot run, out of memory.. what happened
    moscot_transition_row_norm_df = pd.read_csv(moscot_output_dir+os.sep+f'{prev_day}_{next_day}_row_normalized.csv', index_col=0, header=0)
    moscot_transition_col_norm_df = pd.read_csv(moscot_output_dir+os.sep+f'{prev_day}_{next_day}_column_normalized.csv', index_col=0, header=0)
    moscot_transition_min_df = np.minimum(moscot_transition_row_norm_df.values, moscot_transition_col_norm_df.values)
    moscot_transition_df = pd.DataFrame(moscot_transition_min_df, index=moscot_transition_row_norm_df.index, columns=moscot_transition_col_norm_df.columns)
    # DeST-OT
    DeST_OT_transition_row_norm_df = pd.read_csv(DeST_OT_output_dir+os.sep+f'{prev_day}_{next_day}_row_normalized.csv', index_col=0, header=0)
    DeST_OT_transition_col_norm_df = pd.read_csv(DeST_OT_output_dir+os.sep+f'{prev_day}_{next_day}_column_normalized.csv', index_col=0, header=0)
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
comparison_df['Data'] = "Chicken Heart"
comparison_df.to_csv(result_dir+os.sep+'OT_method_comparison.csv')


## --- seaborn barplot for each Transition day
#import seaborn as sns
#import matplotlib.pyplot as plt
#plt.figure(figsize=(6,6))
#sns.barplot(data=comparison_df, x='Day', y='naive_accuracy', hue='method', hue_order=["SpaDOT", "DeST-OT", "Moscot"], palette="Set2")
#plt.ylim(0,1)
#plt.ylabel('Accuracy based on transitions (moscot metrics)', fontsize=12)
#plt.xlabel('Transitions', fontsize=12)
#plt.title('Comparison of OT methods on Chicken Heart data', fontsize=14)
#plt.legend(title='Method', fontsize=12)
#plt.tight_layout()
#plt.savefig(result_dir+os.sep+'OT_method_comparison_naive_accuracy.png', bbox_inches='tight')
#
#plt.figure(figsize=(6,6))
#sns.barplot(data=comparison_df, x='Day', y='Accuracy', hue='method', hue_order=["SpaDOT", "DeST-OT", "Moscot"], palette="Set2")                                                                                                                                                                                
#plt.ylim(0,1)
#plt.ylabel('Accuracy based on transitions (moscot metrics)', fontsize=12)
#plt.xlabel('Transitions', fontsize=12)
#plt.title('Comparison of OT methods on Chicken Heart data', fontsize=14)
#plt.legend(title='Method', fontsize=12)
#plt.tight_layout()
#plt.savefig(result_dir+os.sep+'OT_method_comparison_curated_accuracy.png', bbox_inches='tight')
