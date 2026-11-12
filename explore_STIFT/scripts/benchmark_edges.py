import numpy as np

COMMON = ['VZ', 'CP', 'VLMC', 'Septum', 'Striatum']
CHICKEN_EDGES = [[('Epicardium- like', 'Epicardium'), ('Ventricle', 'Trabecular LV and \nendocardium'), ('Ventricle', 'Right ventricle'), ('Ventricle', 'Compact LV and \ninter-ventricular septum'), ('Atria', 'Atria'), ('Valves', 'Valves')], [(x, x) for x in ['Valves', 'Atria', 'Epicardium', 'Trabecular LV and \nendocardium', 'Right ventricle', 'Compact LV and \ninter-ventricular septum']], [('Epicardium', 'Epicardium'), ('Right ventricle', 'Right ventricle'), ('Trabecular LV and \nendocardium', 'Endocardium'), ('Compact LV and \ninter-ventricular septum', 'Compact LV and \ninter-ventricular septum'), ('Atria', 'Atria'), ('Valves', 'Valves')]]

def edges(dataset, i, labels):
    if dataset == 'MouseEmbryo':
        return {(x, x) for x in labels if i == 0 or x != 'Liver'}
    if dataset == 'ChickenHeart':
        return set(CHICKEN_EDGES[i])
    return {(x, x) for x in COMMON} | (set([('MP&DP&LP', 'MP&DP'), ('MP&DP&LP', 'LP')]) if i == 0 else set([('MP&DP', 'MP'), ('MP&DP', 'DP'), ('LP', 'LP')]))

def sym(x):
    assert (x.sum(axis=0) > 0).all() and (x.sum(axis=1) > 0).all()
    return np.minimum(x / x.sum(axis=1, keepdims=True), x / x.sum(axis=0, keepdims=True))

def score(m, rows, cols, graph):
    mask = np.array([[(a, b) in graph for b in cols] for a in rows])
    return float(m[mask].sum() / m.sum())

def onehot(values, names):
    return (np.asarray(values)[:, None] == np.asarray(names)[None, :]).astype(float)
