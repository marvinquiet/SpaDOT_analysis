import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import *
import sklearn
from sklearn.preprocessing import *
from sklearn.metrics import *
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from tqdm import trange

from OT_loss.ot_solvers import compute_transport_map

# --- analyze domain metrics including clustering metrics and domain continuity metrics
def analyze_domain_metrics(res, result_dir, method):
    timepoints = list(res["timepoint"].unique())
    timepoints.sort()
    res_dict = dict()
    for tp in timepoints:
        tp_res = res[res["timepoint"]==tp]
        # --- domain labels
        tp_ARI = adjusted_rand_score(tp_res["true_labels"], tp_res["pred_labels"])
        tp_NMI = normalized_mutual_info_score(tp_res["true_labels"], tp_res["pred_labels"])
        tp_AMI = adjusted_mutual_info_score(tp_res["true_labels"], tp_res["pred_labels"])
        tp_H = homogeneity_score(tp_res["true_labels"], tp_res["pred_labels"])
        #print(f"{method} at timepoint {tp}: ARI={tp_ARI}, NMI={tp_NMI}, AMI={tp_AMI}, H={tp_H}")
        # --- domain continuity
        loc = np.array(tp_res[["pixel_x", "pixel_y"]])
        pred = np.array(tp_res["pred_labels"])
        tp_ASW = _compute_ASW(loc, pred)
        tp_true_ASW = _compute_ASW(loc, np.array(tp_res["true_labels"]))
        tp_PAS = _compute_PAS(loc, pred)
        tp_true_PAS = _compute_PAS(loc, np.array(tp_res["true_labels"]))
        tp_CHAOS = _compute_CHAOS(loc, pred)
        tp_true_CHAOS = _compute_CHAOS(loc, np.array(tp_res["true_labels"]))
        # tp_SCS = _compute_SCS(loc, pred)
        #print(f"{method} at timepoint {tp}: ASW={tp_ASW}, PAS={tp_PAS}, CHAOS={tp_CHAOS}")
        # save results
        res_dict[tp] = {}
        res_dict[tp]["ARI"], res_dict[tp]["NMI"], res_dict[tp]["AMI"], res_dict[tp]["H"] = tp_ARI, tp_NMI, tp_AMI, tp_H
        # res_dict[tp]["ASW"], res_dict[tp]["PAS"], res_dict[tp]["CHAOS"], res_dict[tp]['SCS'] = tp_ASW, tp_PAS, tp_CHAOS, tp_SCS
        res_dict[tp]["ASW"], res_dict[tp]["PAS"], res_dict[tp]["CHAOS"] = tp_ASW, tp_PAS, tp_CHAOS
        res_dict[tp]["true_ASW"], res_dict[tp]["true_PAS"], res_dict[tp]["true_CHAOS"] = tp_true_ASW, tp_true_PAS, tp_true_CHAOS
    res_df = pd.DataFrame(res_dict)
    res_df['avg'] = res_df.mean(axis=1)
    print(res_df)
    res_df.to_csv(result_dir+os.sep+method+'_domain_metrics.csv')
    return res_df

# --- analyze latent metrics
def analyze_latent_metrics(latent, batch, result_dir, method):
    latent.to_csv(result_dir+os.sep+method+'_latent.csv')
    batch.to_csv(result_dir+os.sep+method+'_batch.csv')
    # run kBET and tLISI
    os.system(f"Rscript latent_metrics_utils.R --result_dir {result_dir} --method {method}")
    # run tASW
    tASW = -silhouette_score(latent, batch)
    with open(result_dir+os.sep+method+'_tASW.txt', 'w') as f:
        f.write(str(tASW)+'\n')

# --- visualize latent
def visualize_latent(adata, use_rep_key, result_dir, method, color_key=['timepoint', 'true_labels', 'pred_labels']):
    # visualize tsne as well
    sc.tl.tsne(adata, use_rep=use_rep_key)
    sc.pl.tsne(adata, color=color_key)
    plt.savefig(result_dir+os.sep+method+'_tsne.png')
    plt.close()
    
    sc.pp.neighbors(adata, use_rep=use_rep_key)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=color_key)
    plt.savefig(result_dir+os.sep+method+'_umap.png')
    plt.close()
    # --- for each timepoint
    timepoints = list(adata.obs['timepoint'].unique())
    timepoints.sort()
    for tp in timepoints:
        tp_adata = adata[adata.obs['timepoint']==tp]
        sc.pl.umap(tp_adata, color=color_key)
        plt.savefig(result_dir+os.sep+method+'_'+str(tp)+'_umap.png')
        plt.close()

# --- calculate and visualize OT transport matrix
def compute_OT_heatmap(kmeans_dict, result_dir, method=""):
    '''Compute OT heatmap
    '''
    ot_config = {
        "growth_iters": 3,
        "epsilon": 0.05, # larger epsilon can result in a blurred transport plan
        # "lambda1": 1, # row imbalance
        "lambda1": 0.1, # row imbalance, larger to make results more balanced along rows
        "lambda2": 5, # column imbalance, larger to make results more balanced along columns
        "epsilon0": 1,
        "tau": 1000,
        "scaling_iter": 3000,
        "inner_iter_max": 50,
        "tolerance": 1e-8,
        "max_iter": 1e7,
        "batch_size": 5,
        "extra_iter": 1000,
        "numItermax": 1000000,
        "use_Py": False,
        "use_C": True,
        "profiling": False,
        "method": "waddington" # Default
    }  # specified in tutorial
    timepoints = list(kmeans_dict.keys())
    timepoints.sort()
    for tp_i in range(len(timepoints)-1):
        cur_tp = timepoints[tp_i]
        next_tp = timepoints[tp_i+1]
        # regularization on C
        C = sklearn.metrics.pairwise.pairwise_distances(kmeans_dict[cur_tp].cluster_centers_, kmeans_dict[next_tp].cluster_centers_, metric="sqeuclidean", n_jobs=-1)
        C = C**2
        C = C / np.median(C)
        gamma = compute_transport_map(kmeans_dict[cur_tp].cluster_centers_, kmeans_dict[next_tp].cluster_centers_, ot_config, C=C, G=None)
        ## write OT matrix
        gamma_df = pd.DataFrame(gamma)
        gamma_df.to_csv(result_dir+os.sep+f'{method}_gamma_tp{cur_tp}_tp{next_tp}.csv')
        plt.figure(figsize=(8,7))
        sns.heatmap(gamma, linewidth=0.5, annot=True, cmap="Reds")
        plt.savefig(result_dir+os.sep+f'{method}_gamma_tp{cur_tp}_tp{next_tp}_final_heatmap.png')
        plt.close()


# --- visualize spatial
### ATACVisium
def edit_ATACVisium_coord(result_df, coord_x_key, coord_y_key, tp_key='timepoint'):
    timepoints = list(result_df[tp_key].unique())
    timepoints.sort()
    tp_results = []
    result_df[coord_x_key+'_bak'] = result_df[coord_x_key]
    result_df[coord_y_key+'_bak'] = result_df[coord_y_key]
    for tp in timepoints:
        tp_res = result_df[result_df[tp_key]==tp]
        if tp == 't1':
            tp_res[coord_x_key] = tp_res[coord_y_key+'_bak']
            tp_res[coord_y_key] = -tp_res[coord_x_key+'_bak']
        elif tp == 't0' or tp == 't2':
            tp_res[coord_x_key] = -tp_res[coord_y_key+'_bak']
            tp_res[coord_y_key] = -tp_res[coord_x_key+'_bak']
        tp_results.append(tp_res)
    result_df = pd.concat(tp_results)
    return result_df

def visualize_ATACVisium_spatial(result_df, coord_x_key, coord_y_key, pred_label_key, result_dir, method, tp_key='timepoint', height=4):
    result_df = edit_ATACVisium_coord(result_df, coord_x_key, coord_y_key, tp_key=tp_key)
    visualize_spatial(result_df, coord_x_key, coord_y_key, pred_label_key, result_dir, method, tp_key=tp_key, height=height)

#### STRS
def edit_STRS_coord(result_df, coord_x_key, coord_y_key, tp_key='timepoint'):
    timepoints = list(result_df[tp_key].unique())
    timepoints.sort()
    tp_results = []
    result_df[coord_x_key+'_bak'] = result_df[coord_x_key]
    result_df[coord_y_key+'_bak'] = result_df[coord_y_key]
    for tp in timepoints:
        tp_res = result_df[result_df[tp_key]==tp]
        if tp == 't1':
            tp_res[coord_x_key] = -tp_res[coord_y_key+'_bak']
            tp_res[coord_y_key] = tp_res[coord_x_key+'_bak']
        elif tp == 't0' or tp == 't2':
            tp_res[coord_x_key] = tp_res[coord_y_key+'_bak']
            tp_res[coord_y_key] = -tp_res[coord_x_key+'_bak']
        tp_results.append(tp_res)
    result_df = pd.concat(tp_results)
    return result_df

def visualize_STRS_spatial(result_df, coord_x_key, coord_y_key, pred_label_key, result_dir, method, tp_key='timepoint', height=4):
    result_df = edit_STRS_coord(result_df, coord_x_key, coord_y_key, tp_key=tp_key)
    visualize_spatial(result_df, coord_x_key, coord_y_key, pred_label_key, result_dir, method, tp_key=tp_key, height=height)


def edit_ChichenHeart_coord(result_df, coord_x_key, coord_y_key, tp_key='timepoint'):
    timepoints = list(result_df[tp_key].unique())
    timepoints.sort()
    tp_results = []
    result_df[coord_x_key+'_bak'] = result_df[coord_x_key]
    result_df[coord_y_key+'_bak'] = result_df[coord_y_key]
    for tp in timepoints:
        tp_res = result_df[result_df[tp_key]==tp]
        tp_res[coord_x_key] = -tp_res[coord_y_key+'_bak']
        tp_res[coord_y_key] = -tp_res[coord_x_key+'_bak']
        tp_results.append(tp_res)
    result_df = pd.concat(tp_results)
    return result_df

def visualize_ChickenHeart_spatial(result_df, coord_x_key, coord_y_key, pred_label_key, result_dir, method, tp_key='timepoint', height=4):
    result_df = edit_ChichenHeart_coord(result_df, coord_x_key, coord_y_key, tp_key=tp_key)
    visualize_spatial(result_df, coord_x_key, coord_y_key, pred_label_key, result_dir, method, tp_key=tp_key, height=height)


#### general visualization
def visualize_spatial(result_df, coord_x_key, coord_y_key, pred_label_key, result_dir, method, tp_key='timepoint',
                              height=4, color_palette=None):
    timepoints = list(result_df[tp_key].unique())
    timepoints.sort()
    g = sns.FacetGrid(result_df, col=tp_key, sharex=False, sharey=False, col_order=timepoints,
                      height=height, aspect=0.8)  # Facet the grid by "timepoint"
    if color_palette is None:
        g.map_dataframe(sns.scatterplot, x=coord_x_key, y=coord_y_key, hue=pred_label_key, s=3)
    else:
        g.map_dataframe(sns.scatterplot, x=coord_x_key, y=coord_y_key, hue=pred_label_key, s=3, palette=color_palette)
    g.set(xticklabels=[], yticklabels=[], xlabel=None, ylabel=None)  
    g.tick_params(bottom=False, left=False)
    # g.add_legend()
    for ax in g.axes.flat:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc="upper center",
              bbox_to_anchor=(0.5, -0.15), ncol=3, frameon=False)
    plt.tight_layout()
    plt.savefig(result_dir+os.sep+method+'_spatial.png')
    plt.close()

# rename true labels to indicate domain labels
def rename_GT_labels(row):
    if row['timepoint'] == 't1':
        if row['true_labels'] == 1:
            return "2"
        if row['true_labels'] == 2:
            return "1-2"
        if row['true_labels'] == 3:
            return "1-1"
        if row['true_labels'] == 4:
            return "3"
    if row['timepoint'] == 't2':
        if row['true_labels'] == 1:
            return "2"
        if row['true_labels'] == 2:
            return "3-1"
        if row['true_labels'] == 3:
            return "1-2"
        if row['true_labels'] == 4:
            return "1-1"
        if row['true_labels'] == 5:
            return "3-2"
    return str(row['true_labels'])

GT_custom_palette = {
    '1': '#EA3423',
    '2': '#76FC4C',
    '3': '#FEFC53',
    '1-1': '#891AF5',
    '1-2': '#EA33F4',
    '3-1': '#716A1E',
    '3-2': '#F5C242'
}

# --- quantify data complexity
def _get_spatial_entropy(C, C_sum):
    H = 0
    for i in range(len(C)):
        for j in range(len(C)):
            z = C[i, j]
            if z != 0:
                H += -(z / C_sum) * np.log(z / C_sum)
    return H

def spatial_entropy(k_neighbors, labels, degree=4):
    """
    Calculates spatial entropy of graph
    """
    S = np.broadcast_to(labels[:, None], (len(labels), degree))
    N = labels[k_neighbors]
    cluster_names = np.unique(labels)
    cluster_nums = len(cluster_names)
    C = np.zeros((cluster_nums, cluster_nums))
    for i in range(cluster_nums):
        for j in range(cluster_nums):
            # C[i, j] = np.sum(np.logical_and(N == i, S == j))
            C[i, j] = np.sum(np.logical_and(S == cluster_names[i], N == cluster_names[j]))
    C_sum = C.sum()
    return _get_spatial_entropy(C, C_sum)

def spatial_coherence_score(adata, annotation_key, degree=4, rep_time=1000, seed=0):
    spatial_coords = adata.obsm['spatial']
    origin_labels = adata.obs[annotation_key].values
    # Use kneighbors_graph to get the adjacency matrix
    neigh = NearestNeighbors(n_neighbors=degree, metric='euclidean').fit(spatial_coords)
    k_neighbors = neigh.kneighbors(n_neighbors=degree, return_distance=False)
    true_entropy = spatial_entropy(k_neighbors, origin_labels, degree=degree)
    entropies = []
    rng = np.random.default_rng(seed)
    shuffled_labels = origin_labels.copy()
    for _ in trange(rep_time):
        rng.shuffle(shuffled_labels)
        entropies.append(spatial_entropy(k_neighbors, shuffled_labels, degree=degree))
    return (true_entropy - np.mean(entropies)) / np.std(entropies), true_entropy, entropies

# --- domain continuity metrics
# https://github.com/maiziezhoulab/BenchmarkST/blob/main/metrics/coherence_metrics.py
def _compute_CHAOS(location, clusterlabel):
    """
    Calculate the CHAOS score for a given set of spatial coordinates and predicted labels.
    param: X - spatial coordinates
    param: pred_labels - predicted labels
    return: CHAOS score
    """
    # Standardize the spatial coordinates
    location = StandardScaler().fit_transform(location)
    # Get the unique cluster labels
    cluster_labels = np.unique(clusterlabel)
    # Initialize the distance value and count
    dist_val = 0.
    count = 0
    # Iterate through each cluster
    for k in cluster_labels:
        # Get the spatial coordinates for the current cluster
        cluster_coords = location[clusterlabel == k, :]
        # Check if there are at least 2 spatial coordinates in the cluster
        if len(cluster_coords) <= 2:
            continue
        else:
            count += len(cluster_coords)
        # Calculate the distance to the nearest neighbor for each spatial coordinate in the cluster
        nbrs = NearestNeighbors(n_neighbors=1).fit(cluster_coords)
        distances, _ = nbrs.kneighbors()
        # Sum the distances
        dist_val = dist_val + np.sum(distances)
    # Calculate the CHAOS score
    return dist_val / count

def _compute_PAS(location, clusterlabel, k=10):
    """
    Calculate the PAS score for a given set of spatial coordinates and predicted labels.
    param: X - spatial coordinates
    param: pred_labels - predicted labels
    param: k - number of nearest neighbors to consider
    return: PAS score
    """
    # Use NearestNeighbors to find the nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=k).fit(location)
    indices = nbrs.kneighbors(return_distance=False)
    # Calculate the PAS score
    return ((clusterlabel.reshape(-1, 1) != clusterlabel[indices]).sum(1) > k / 2).mean()

def _compute_ASW(location, clusterlabel):
    return silhouette_score(X=location,labels=clusterlabel)

def _compute_SCS(location, clusterlabel, degree=4, rep_time=1000, seed=0):
    neigh = NearestNeighbors(n_neighbors=degree, metric='euclidean').fit(location)
    k_neighbors = neigh.kneighbors(n_neighbors=degree, return_distance=False)
    true_entropy = spatial_entropy(k_neighbors, clusterlabel, degree=degree)
    entropies = []
    rng = np.random.default_rng(seed)
    shuffled_labels = clusterlabel.copy()
    for _ in trange(rep_time):
        rng.shuffle(shuffled_labels)
        entropies.append(spatial_entropy(k_neighbors, shuffled_labels, degree=degree))
    return (true_entropy - np.mean(entropies)) / np.std(entropies)