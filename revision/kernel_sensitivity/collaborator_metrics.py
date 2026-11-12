"""Metric functions. graphs(), accuracy(), evaluate_using_curated_transitions() and min_norm() are verbatim from the
manuscript benchmark (SpaDOT_analysis/summarize_benchmark/scripts/benchmark_OT_transition_table.py); PAS and CHAOS are from
summarize_benchmark/scripts/utils.py."""
import numpy as np, pandas as pd, networkx as nx
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
norm = lambda s: str(s).replace("\n", " ")

def graphs(ds, regions):
    """Curated region correspondences of the manuscript (node names newline-normalised)."""
    if ds == "MouseEmbryo":
        g01 = nx.DiGraph(); g12 = nx.DiGraph()
        for r in regions:
            g01.add_node(r); g01.add_edge(r, r); g12.add_node(r)
            if r != "Liver": g12.add_edge(r, r)
        return [g01, g12]
    if ds == "ChickenHeart":
        g01 = nx.DiGraph(); g01.add_nodes_from(map(norm, ['Epicardium- like', 'Outflow tract', 'Ventricle', 'Valves', 'Atria', 'Epicardium', 'Endothelium', 'Trabecular LV and \nendocardium', 'Right ventricle', 'Compact LV and \ninter-ventricular septum']))
        g01.add_edges_from([(norm(a), norm(b)) for a, b in [('Epicardium- like', 'Epicardium'), ('Ventricle', 'Trabecular LV and \nendocardium'), ('Ventricle', 'Right ventricle'), ('Ventricle', 'Compact LV and \ninter-ventricular septum'), ('Atria', 'Atria'), ('Valves', 'Valves')]])
        g12 = nx.DiGraph(); g12.add_nodes_from(map(norm, ['Valves', 'Atria', 'Epicardium', 'Endothelium', 'Trabecular LV and \nendocardium', 'Right ventricle', 'Compact LV and \ninter-ventricular septum', 'Outflow tract']))
        g12.add_edges_from([(norm(a), norm(b)) for a, b in [('Valves', 'Valves'), ('Atria', 'Atria'), ('Epicardium', 'Epicardium'), ('Trabecular LV and \nendocardium', 'Trabecular LV and \nendocardium'), ('Right ventricle', 'Right ventricle'), ('Compact LV and \ninter-ventricular septum', 'Compact LV and \ninter-ventricular septum')]])
        g23 = nx.DiGraph(); g23.add_nodes_from(map(norm, ['Epicardium', 'Endocardium', 'Right ventricle', 'Compact LV and \ninter-ventricular septum', 'Atria', 'Valves']))
        g23.add_edges_from([(norm(a), norm(b)) for a, b in [('Epicardium', 'Epicardium'), ('Right ventricle', 'Right ventricle'), ('Trabecular LV and \nendocardium', 'Endocardium'), ('Compact LV and \ninter-ventricular septum', 'Compact LV and \ninter-ventricular septum'), ('Atria', 'Atria'), ('Valves', 'Valves')]])
        return [g01, g12, g23]
    if ds == "Axolotl":
        common = ['VZ', 'CP', 'VLMC', 'Septum', 'Striatum']
        g01 = nx.DiGraph(); g12 = nx.DiGraph()
        for r in common: g01.add_edge(r, r); g12.add_edge(r, r)
        g01.add_edge('MP&DP&LP', 'MP&DP'); g01.add_edge('MP&DP&LP', 'LP')
        g12.add_edge('MP&DP', 'MP'); g12.add_edge('MP&DP', 'DP'); g12.add_edge('LP', 'LP')
        return [g01, g12]
    raise ValueError(f"no curated graph for {ds}")

def graphs_from_edges(edge_lists):
    out = []
    for edges in edge_lists:
        g = nx.DiGraph(); g.add_edges_from([(str(a), str(b)) for a, b in edges]); out.append(g)
    return out

def transition_graphs(ds, spec, regions):
    if "curated" in spec: return graphs(spec["curated"], regions)
    if "edges" in spec: return graphs_from_edges(spec["edges"])
    raise ValueError(f"{ds}: transitions must have 'curated' or 'edges'")

def accuracy(df, G):
    """naive accuracy = mass on known edges / total mass"""
    correct = wrong = 0; edges = list(G.edges)
    for ci in df.index:
        for cj in df.columns:
            if (ci, cj) in edges: correct += df[cj][ci]
            else: wrong += df[cj][ci]
    return correct / (correct + wrong)

def evaluate_using_curated_transitions(df, G, cutoff=0.05):
    """curated accuracy: entries > cutoff only; same-name counts as known"""
    known = unknown = 0; M = df.values; edges = list(G.edges)
    for i, p in enumerate(df.index):
        for j, c in enumerate(df.columns):
            w = M[i, j]
            if w > cutoff:
                if (p, c) in edges or c == p: known += w
                else: unknown += w
    return known / (known + unknown), known + unknown

def min_norm(T):
    """T: raw coupling mass (rows = source, cols = target) -> min(row-normalised, column-normalised)"""
    X = T.values.astype(float)
    return pd.DataFrame(np.minimum(X / X.sum(0, keepdims=True), X / X.sum(1, keepdims=True)), index=T.index, columns=T.columns)

def compute_CHAOS(location, clusterlabel):
    location = StandardScaler().fit_transform(location); dist_val = 0.; count = 0
    for k in np.unique(clusterlabel):
        cc = location[clusterlabel == k, :]
        if len(cc) <= 2: continue
        count += len(cc); dist_val += np.sum(NearestNeighbors(n_neighbors=1).fit(cc).kneighbors()[0])
    return dist_val / count

def compute_PAS(location, clusterlabel, k=10):
    idx = NearestNeighbors(n_neighbors=k).fit(location).kneighbors(return_distance=False)
    return ((clusterlabel.reshape(-1, 1) != clusterlabel[idx]).sum(1) > k / 2).mean()

def flat_corr(A, B):
    """Pearson r between two transition tables over the union of their rows/columns (missing entries = 0)"""
    idx = sorted(set(A.index) | set(B.index)); col = sorted(set(A.columns) | set(B.columns))
    a = A.reindex(index=idx, columns=col).fillna(0).values.ravel(); b = B.reindex(index=idx, columns=col).fillna(0).values.ravel()
    return np.corrcoef(a, b)[0, 1]
