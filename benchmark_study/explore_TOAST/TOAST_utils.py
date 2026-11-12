
import numpy as np
import pandas as pd
from spatial_OT.OT import *
from spatial_OT.utils import *
from scipy.spatial import distance
import numpy as np
import scanpy as sc
import anndata as ad
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import jensenshannon, cdist
import ot



def compute_spatial_fgw_alignment(slice1, slice2, alpha, epsilon, n_comps, k, annotation_col):
    """Compute FGW and FGW-SN alignment and return accuracy and JS divergence scores."""
    joint_adata = sc.concat([slice1, slice2])
    sc.pp.normalize_total(joint_adata, inplace=True)
    sc.pp.log1p(joint_adata)
    sc.pp.pca(joint_adata, n_comps=n_comps)
    joint_datamatrix = joint_adata.obsm['X_pca']

    X = joint_datamatrix[:slice1.shape[0], :]
    Y = joint_datamatrix[slice1.shape[0]:, :]
    
    # Compute spatial graphs
    coords1 = pd.DataFrame(slice1.obsm["spatial"], columns=["x", "y"])
    X_df = pd.DataFrame(X, columns=[f"PC{i+1}" for i in range(X.shape[1])])
    X_df["x"], X_df["y"], X_df["cell_type"] = coords1["x"].values, coords1["y"].values, slice1.obs[annotation_col].values

    G1 = build_knn_graph_from2d(X_df, k=k)
    X_df["spatial_entropy"] = X_df.index.map(compute_spatial_entropy(G1))
    slice1_avg_expr = compute_average_neighbor_expression(G1, pd.DataFrame(X))
    
    coords2 = pd.DataFrame(slice2.obsm["spatial"], columns=["x", "y"])
    Y_df = pd.DataFrame(Y, columns=[f"PC{i+1}" for i in range(Y.shape[1])])
    Y_df["x"], Y_df["y"], Y_df["cell_type"] = coords2["x"].values, coords2["y"].values, slice2.obs[annotation_col].values

    G2 = build_knn_graph_from2d(Y_df, k=k)
    Y_df["spatial_entropy"] = Y_df.index.map(compute_spatial_entropy(G2))
    slice2_avg_expr = compute_average_neighbor_expression(G2, pd.DataFrame(Y))
    
    # Compute cost matrices
    M = distance.cdist(X, Y).astype(float)
    C1 = distance.cdist(slice1.obsm["spatial"], slice1.obsm["spatial"]).astype(float)
    C2 = distance.cdist(slice2.obsm["spatial"], slice2.obsm["spatial"]).astype(float)
    C3 = np.abs(X_df["spatial_entropy"].values[:, np.newaxis] - Y_df["spatial_entropy"].values[np.newaxis, :])
    C4 = distance.cdist(slice1_avg_expr.values, slice2_avg_expr.values).astype(float)

    # Normalize matrices
    for mat in [M, C1, C2, C3, C4]:
        mat /= mat.max() if mat.max() > 0 else 1  # Avoid division by zero
    
    # Compute transport maps
    p, q = ot.unif(X.shape[0]), ot.unif(Y.shape[0])
    G0 = np.outer(p, q)

    FGW_SN = compute_transport(G0, epsilon, alpha, C1, C2, p, q, M, C3, C4)
    return FGW_SN