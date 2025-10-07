import os, glob
import anndata
import pandas as pd

import scanpy as sc
import matplotlib.pyplot as plt

# --- load my scripts
from utils import *

# random_seed = 1234
# np.random.seed(random_seed)
# import rpy2.robjects as robjects
# robjects.r.library("mclust")

# import rpy2.robjects.numpy2ri
# rpy2.robjects.numpy2ri.activate()
# r_random_seed = robjects.r['set.seed']
# r_random_seed(random_seed)
# rmclust = robjects.r['Mclust']

def analyze_simulation_complexity(data_dir, data_prefix):
    # load data
    import random
    adata = anndata.read_h5ad(data_dir+os.sep+data_prefix+os.sep+'adata.h5ad')
    random.seed(1234)
    sampled_idx = random.sample(range(adata.shape[0]), 3000*3)
    adata = adata[sampled_idx].copy()
    adata.obs.rename(columns={'x_coor': 'pixel_x', 'y_coor': 'pixel_y', 'label': 'annotation'}, inplace=True)
    adata.obsm['spatial'] = adata.obs[['pixel_x', 'pixel_y']].values
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()
    res_dict = dict()
    for tp in timepoints:
        tp_adata = adata[adata.obs['timepoint'] == tp]
        tp_SCS, tp_true_entropy, tp_entropies = spatial_coherence_score(tp_adata, annotation_key='annotation')
        print(f"Timepoint {tp}: SCS={tp_SCS}, true_entropy={tp_true_entropy}")
        res_dict[tp] = {}
        res_dict[tp]["SCS"] = tp_SCS
        res_dict[tp]["true_entropy"] = tp_true_entropy
        # res_dict[tp]["tp_entropies"] = tp_entropies
    res_df = pd.DataFrame(res_dict)
    res_df.to_csv(data_dir+os.sep+data_prefix+os.sep+'simulation_complexity_SCS.csv')

def analyze_ATACVisium_complexity(data_dir):
    # --- load data
    section_dict = {1: 't2', 2: 't1', 3: 't0'}
    data_dict = {1: 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz', # E15.5
                2: 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz', # E13.5
                3: 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz' # E12.5
                }
    adata = anndata.read_h5ad(data_dir+os.sep+'exprs.h5ad')
    adata.obs['x'] = None
    adata.obs['y'] = None
    adata.obs['pixel_x'] = None
    adata.obs['pixel_y'] = None
    loc_metadata_list = []
    for key, value in data_dict.items():
        metadata = pd.read_csv(data_dir+os.sep+value, sep=',', index_col=0, header=None)
        metadata.columns = ['tissue', 'y', 'x', 'pixel_y', 'pixel_x']
        metadata['barcodes'] = metadata.index
        metadata.index = metadata['barcodes']+'_'+str(key)
        key_barcodes = adata[adata.obs['section'] == key].obs_names.tolist()
        common_barcodes = list(set(metadata.index.tolist()).intersection(set(key_barcodes)))
        section_metadata = metadata.loc[common_barcodes, :]
        loc_metadata_list.append(section_metadata)
    loc_metadata = pd.concat(loc_metadata_list, axis=0)
    common_barcodes = list(set(loc_metadata.index.tolist()).intersection(set(adata.obs_names.tolist())))
    adata = adata[common_barcodes].copy()
    adata.obs.loc[common_barcodes, 'pixel_x'] = loc_metadata.loc[common_barcodes, 'pixel_x'].values
    adata.obs.loc[common_barcodes, 'pixel_y'] = loc_metadata.loc[common_barcodes, 'pixel_y'].values
    adata.obs['annotation'] = adata.obs['seurat_clusters_harmony']
    adata.obs['annotation'].replace({3: 'CNS', 2: 'CNS', 4: 'CNS',
                                    1: 'Men/PNS', 9: 'Men/PNS', 11: 'Men/PNS',
                                    0: 'Mesenchyme', 5: 'Mesenchyme', 6: 'Mesenchyme', 8: 'Mesenchyme',
                                    7: 'Liver', 10: 'Heart'}, inplace=True)
    adata.obs['annotation'] = adata.obs['annotation'].astype('category')
    adata.obs['timepoint'] = [section_dict[_] for _ in adata.obs['section']]
    adata.obsm['spatial'] = adata.obs[['pixel_x', 'pixel_y']].values
    # remove those cells that are not in the image
    with open(data_dir+os.sep+'t2_removed_cells.txt', 'r') as f:
        t2_removed_cells = f.readlines()
    t2_removed_cells = [_.strip() for _ in t2_removed_cells]
    adata = adata[~adata.obs_names.isin(t2_removed_cells)]
    # analyze complexity
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()
    res_dict = dict()
    for tp in timepoints:
        tp_adata = adata[adata.obs['timepoint'] == tp]
        tp_SCS, tp_true_entropy, tp_entropies = spatial_coherence_score(tp_adata, annotation_key='annotation')
        print(f"Timepoint {tp}: SCS={tp_SCS}, true_entropy={tp_true_entropy}")
        res_dict[tp] = {}
        res_dict[tp]["SCS"] = tp_SCS
        res_dict[tp]["true_entropy"] = tp_true_entropy
        # res_dict[tp]["tp_entropies"] = tp_entropies
    res_df = pd.DataFrame(res_dict)
    res_df.to_csv(data_dir+os.sep+'complexity_SCS.csv')

def analyze_STRS_complexity(data_dir):
    # --- load data
    import scipy
    exprs = scipy.io.mmread(data_dir+os.sep+'STRS.mtx.gz').toarray()
    genes = np.loadtxt(data_dir+os.sep+'STRS_genes.txt', dtype='str')
    metadata = pd.read_csv(data_dir+os.sep+'STRS_metadata.csv', header=0, index_col=0)
    adata = anndata.AnnData(X=exprs.T, obs=metadata, var=pd.DataFrame(index=genes, data={'gene_id':genes}))
    adata.obs.rename(columns={'injury.zones': 'annotation'}, inplace=True)
    loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
    adata.obsm['spatial'] = loc
    # --- end loading data
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()
    res_dict = dict()
    for tp in timepoints:
        tp_adata = adata[adata.obs['timepoint'] == tp]
        tp_SCS, tp_true_entropy, tp_entropies = spatial_coherence_score(tp_adata, annotation_key='annotation')
        print(f"Timepoint {tp}: SCS={tp_SCS}, true_entropy={tp_true_entropy}")
        res_dict[tp] = {}
        res_dict[tp]["SCS"] = tp_SCS
        res_dict[tp]["true_entropy"] = tp_true_entropy
        # res_dict[tp]["tp_entropies"] = tp_entropies
    res_df = pd.DataFrame(res_dict)
    res_df.to_csv(data_dir+os.sep+'complexity_SCS.csv')

# Generate clustering metrics for different methods: SpatialPCA, GraphST, SEDR, spaGCN, spaVAE, STAligner, PRECAST, spatiAlign, SPIRAL, STADIAN, spaDOT
def analyze_spaDOT_domains(result_dir, result_prefix,
                           label_key = "annotation", 
                           rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"},
                           method='spaDOT', latent_index=None):
    # --- spaDOT --- 
    spaDOT_results_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
    spaDOT_adata = anndata.read_h5ad(spaDOT_results_dir+os.sep+'final_latent.h5ad')
    if latent_index is not None:
        spaDOT_adata = spaDOT_adata[:, latent_index].copy()
    spaDOT_results_list = []
    kmeans_dict = dict()
    for tp in set(spaDOT_adata.obs['timepoint']):
        tp_adata = spaDOT_adata[spaDOT_adata.obs['timepoint']==tp]
        tp_emb = tp_adata.X
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        # run mclust_R, does not have a better performance
        # tp_mclust = rmclust(tp_emb, len(set(tp_adata.obs[label_key])), 'EEE') # EEE is the model name
        # tp_adata.obs['mclust'] = np.array(tp_mclust[-2]).astype('int').astype(str)
        spaDOT_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    spaDOT_results = pd.concat(spaDOT_results_list)
    spaDOT_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(spaDOT_results, spaDOT_results_dir, method)
    try:
        compute_OT_heatmap(kmeans_dict, spaDOT_results_dir, method=method)
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    spaDOT_adata.obs = spaDOT_results.loc[spaDOT_adata.obs_names.tolist(), :]
    spaDOT_adata.write_h5ad(spaDOT_results_dir+os.sep+'prediction_adata.h5ad')
    return spaDOT_adata


def analyze_spatialPCA_domains(result_dir, result_prefix, 
                               label_key="label",
                               rename_dict={"label": "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"}):
    # --- SpatialPCA --- 
    spatialPCA_results_dir = result_dir+os.sep+"explore_spatialPCA"+os.sep+result_prefix
    spatialPCA_results_df = pd.read_csv(spatialPCA_results_dir+os.sep+'spatialPCA_results.csv', index_col=0)
    spatialPCA_results_list = []
    kmeans_dict = dict()
    for tp in set(spatialPCA_results_df['timepoint']):
        tp_spatialPCA_results = spatialPCA_results_df[spatialPCA_results_df['timepoint']==tp]
        tp_emb = tp_spatialPCA_results.iloc[:, :20].to_numpy()
        tp_kmeans = KMeans(n_clusters=len(set(tp_spatialPCA_results[label_key])), random_state=0).fit(tp_emb)
        tp_spatialPCA_results['kmeans'] = tp_kmeans.labels_.astype(str)
        spatialPCA_results_list.append(tp_spatialPCA_results)
        kmeans_dict[tp] = tp_kmeans
    spatialPCA_results = pd.concat(spatialPCA_results_list)
    spatialPCA_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(spatialPCA_results, spatialPCA_results_dir, "spatialPCA")
    try:
        compute_OT_heatmap(kmeans_dict, spatialPCA_results_dir, method='spatialPCA')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    spatialPCA_adata = anndata.AnnData(spatialPCA_results.iloc[:, :20].to_numpy(), obs=spatialPCA_results)
    spatialPCA_adata.write_h5ad(spatialPCA_results_dir+os.sep+'prediction_adata.h5ad')
    return spatialPCA_adata

def analyze_GraphST_domains(result_dir, result_prefix,
                            label_key="annotation",
                            rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"}):
    # --- GraphST ---
    GraphST_results_dir = result_dir+os.sep+"explore_GraphST"+os.sep+result_prefix
    GraphST_adata = anndata.read_h5ad(GraphST_results_dir+os.sep+'GraphST_adata.h5ad')
    GraphST_results_list = []
    kmeans_dict = dict()
    for tp in set(GraphST_adata.obs['timepoint']):
        tp_adata = GraphST_adata[GraphST_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['emb']
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        GraphST_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    GraphST_results = pd.concat(GraphST_results_list)
    GraphST_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(GraphST_results, GraphST_results_dir, "GraphST")
    try:
        compute_OT_heatmap(kmeans_dict, GraphST_results_dir, method='GraphST')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    GraphST_adata.obs = GraphST_results.loc[GraphST_adata.obs_names.tolist(), :]
    GraphST_adata.write_h5ad(GraphST_results_dir+os.sep+'prediction_adata.h5ad')
    return GraphST_adata

def analyze_SEDR_domains(result_dir, result_prefix,
                         label_key="annotation",
                         rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"}):
    # --- SEDR ---
    SEDR_results_dir = result_dir+os.sep+"explore_SEDR"+os.sep+result_prefix
    SEDR_adata = anndata.read_h5ad(SEDR_results_dir+os.sep+'SEDR_adata.h5ad')
    SEDR_results_list = []
    kmeans_dict = dict()
    for tp in set(SEDR_adata.obs['timepoint']):
        tp_adata = SEDR_adata[SEDR_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['SEDR.Harmony']
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        SEDR_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    SEDR_results = pd.concat(SEDR_results_list)
    SEDR_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(SEDR_results, SEDR_results_dir, "SEDR")
    # capture Runtime Error
    try:
        compute_OT_heatmap(kmeans_dict, SEDR_results_dir, method='SEDR')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    SEDR_adata.obs = SEDR_results.loc[SEDR_adata.obs_names.tolist(), :]
    SEDR_adata.write_h5ad(SEDR_results_dir+os.sep+'prediction_adata.h5ad')
    return SEDR_adata

def analyze_spaGCN_domains(result_dir, result_prefix,
                           label_key="annotation",
                           rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"}):
    # --- spaGCN ---
    spaGCN_results_dir = result_dir+os.sep+"explore_spaGCN"+os.sep+result_prefix
    spaGCN_adata = anndata.read_h5ad(spaGCN_results_dir+os.sep+'spaGCN.h5ad')
    spaGCN_results_list = []
    kmeans_dict = dict()
    for tp in set(spaGCN_adata.obs['timepoint']):
        tp_adata = spaGCN_adata[spaGCN_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['spaGCN_embed']
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        spaGCN_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    spaGCN_results = pd.concat(spaGCN_results_list)
    spaGCN_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(spaGCN_results, spaGCN_results_dir, "spaGCN")
    try:
        compute_OT_heatmap(kmeans_dict, spaGCN_results_dir, method='spaGCN')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    spaGCN_adata.obs = spaGCN_results.loc[spaGCN_adata.obs_names.tolist(), :]
    spaGCN_adata.write_h5ad(spaGCN_results_dir+os.sep+'prediction_adata.h5ad')
    return spaGCN_adata

def analyze_spaVAE_domains(result_dir, result_prefix,
                           label_key="annotation",
                           rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"}):
    # --- spaVAE ---
    spaVAE_results_dir = result_dir+os.sep+"explore_spaVAE"+os.sep+result_prefix
    metadata_files = [f for f in os.listdir(spaVAE_results_dir) if f.endswith("_metadata.csv")]
    assert len(metadata_files) == 1
    metadata = pd.read_csv(spaVAE_results_dir+os.sep+metadata_files[0], index_col=0)
    latent = np.loadtxt(spaVAE_results_dir+os.sep+'final_latent.txt', delimiter=',')
    latent_adata = anndata.AnnData(X=latent, obs=metadata)
    latent_adata.obs[label_key] = latent_adata.obs[label_key].astype('category')
    spaVAE_results_list = []
    kmeans_dict = dict()
    for tp in set(latent_adata.obs['timepoint']):
        tp_adata = latent_adata[latent_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.X
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        spaVAE_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    spaVAE_results = pd.concat(spaVAE_results_list)
    spaVAE_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(spaVAE_results, spaVAE_results_dir, "spaVAE")
    try:
        compute_OT_heatmap(kmeans_dict, spaVAE_results_dir, method='spaVAE')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    latent_adata.obs = spaVAE_results.loc[latent_adata.obs_names.tolist(), :]
    latent_adata.write_h5ad(spaVAE_results_dir+os.sep+'prediction_adata.h5ad')
    return latent_adata

def analyze_STAligner_domains(result_dir, result_prefix,
                              label_key="annotation",
                              rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"}):
    # --- STAligner ---
    STAligner_results_dir = result_dir+os.sep+"explore_STAligner"+os.sep+result_prefix
    STAligner_adata = anndata.read_h5ad(STAligner_results_dir+os.sep+'STAligner_adata.h5ad')
    STAligner_results_list = []
    kmeans_dict = dict()
    for tp in set(STAligner_adata.obs['timepoint']):
        tp_adata = STAligner_adata[STAligner_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['STAligner']
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        STAligner_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    STAligner_results = pd.concat(STAligner_results_list)
    STAligner_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(STAligner_results, STAligner_results_dir, "STAligner")
    try:
        compute_OT_heatmap(kmeans_dict, STAligner_results_dir, method='STAligner')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    STAligner_adata.obs = STAligner_results.loc[STAligner_adata.obs_names.tolist(), :]
    STAligner_adata.write_h5ad(STAligner_results_dir+os.sep+'prediction_adata.h5ad')
    return STAligner_adata

def analyze_PRECAST_domains(result_dir, result_prefix,
                            label_key="annotation",
                            rename_dict={"label": "true_labels", "kmeans": "pred_labels", "row": "pixel_x", "col": "pixel_y"}
                            ):
    # --- PRECAST ---
    PRECAST_results_dir = result_dir+os.sep+"explore_PRECAST"+os.sep+result_prefix
    PRECAST_results_df = pd.read_csv(PRECAST_results_dir+os.sep+'PRECAST_results.csv', index_col=0)
    PRECAST_results_list = []
    kmeans_dict = dict()
    for tp in set(PRECAST_results_df['timepoint']):
        tp_PRECAST_results = PRECAST_results_df[PRECAST_results_df['timepoint']==tp]
        tp_emb = tp_PRECAST_results.iloc[:, :15].to_numpy()
        tp_kmeans = KMeans(n_clusters=len(set(tp_PRECAST_results[label_key])), random_state=0).fit(tp_emb)
        tp_PRECAST_results['kmeans'] = tp_kmeans.labels_.astype(str)
        PRECAST_results_list.append(tp_PRECAST_results)
        kmeans_dict[tp] = tp_kmeans
    PRECAST_results = pd.concat(PRECAST_results_list)
    PRECAST_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(PRECAST_results, PRECAST_results_dir, "PRECAST")
    try:
        compute_OT_heatmap(kmeans_dict, PRECAST_results_dir, method="PRECAST")
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    PRECAST_adata = anndata.AnnData(PRECAST_results.iloc[:, :15].to_numpy(), obs=PRECAST_results)
    PRECAST_adata.write_h5ad(PRECAST_results_dir+os.sep+'prediction_adata.h5ad')
    return PRECAST_adata

def analyze_spatiAlign_domains(result_dir, result_prefix,
                               label_key="annotation",
                               rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"}):
    # --- spatiAlign ---
    spatiAlign_results_dir = result_dir+os.sep+"explore_spatiAlign"+os.sep+result_prefix
    # list all res/correct_data*.h5ad
    correct_adata_files = glob.glob(spatiAlign_results_dir+os.sep+"res/correct_data*.h5ad")
    correct_adata_list = list()
    for correc_adata_file in correct_adata_files:
        correct_adata = sc.read_h5ad(correc_adata_file)
        correct_adata.var['genes'] = correct_adata.var_names.tolist()
        correct_adata_list.append(correct_adata)
    merge_data = anndata.AnnData.concatenate(*correct_adata_list)
    spatiAlign_results_list = []
    kmeans_dict = dict()
    for tp in set(merge_data.obs['timepoint']):
        tp_adata = merge_data[merge_data.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['correct']
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        spatiAlign_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    spatiAlign_results = pd.concat(spatiAlign_results_list)
    spatiAlign_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(spatiAlign_results, spatiAlign_results_dir, "spatiAlign")
    try:
        compute_OT_heatmap(kmeans_dict, spatiAlign_results_dir, method="spatiAlign")
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    merge_data.obs = spatiAlign_results.loc[merge_data.obs_names.tolist(), :]
    merge_data.write_h5ad(spatiAlign_results_dir+os.sep+'prediction_adata.h5ad')
    return merge_data

def analyze_SPIRAL_domains(result_dir, result_prefix,
                           label_key="annotation",
                           rename_dict={"annotation": "true_labels", "kmeans": "pred_labels"}):
    # --- SPIRAL ---
    SPIRAL_results_dir = result_dir+os.sep+"explore_SPIRAL"+os.sep+result_prefix
    SPIRAL_adata = anndata.read_h5ad(SPIRAL_results_dir+os.sep+'SPIRAL_corrected.h5ad')
    SPIRAL_results_list = []
    kmeans_dict = dict()
    for tp in set(SPIRAL_adata.obs['timepoint']):
        tp_adata = SPIRAL_adata[SPIRAL_adata.obs['timepoint'] == tp]
        tp_emb = tp_adata.obsm['spiral']
        tp_kmeans = KMeans(n_clusters=len(set(tp_adata.obs[label_key])), random_state=0).fit(tp_emb)
        tp_adata.obs['kmeans'] = tp_kmeans.labels_.astype(str)
        SPIRAL_results_list.append(tp_adata.obs)
        kmeans_dict[tp] = tp_kmeans
    SPIRAL_results = pd.concat(SPIRAL_results_list)
    SPIRAL_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(SPIRAL_results, SPIRAL_results_dir, "SPIRAL")
    try:
        compute_OT_heatmap(kmeans_dict, SPIRAL_results_dir, method='SPIRAL')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    SPIRAL_adata.obs = SPIRAL_results.loc[SPIRAL_adata.obs_names.tolist(), :]
    SPIRAL_adata.write_h5ad(SPIRAL_results_dir+os.sep+'prediction_adata.h5ad')
    return SPIRAL_adata

def analyze_STADIA_domains(result_dir, result_prefix,
                           label_key="label",
                           rename_dict={"label": "true_labels", "cluster": "pred_labels", "row": "pixel_x", "col": "pixel_y"}):
    # --- STADIA ---
    STADIA_results_dir = result_dir+os.sep+"explore_STADIA"+os.sep+result_prefix
    STADIA_results_df = pd.read_csv(STADIA_results_dir+os.sep+'STADIA_results.csv', index_col=0)
    STADIA_results_list = []
    kmeans_dict = dict()
    for tp in set(STADIA_results_df['timepoint']):
        tp_STADIA_results = STADIA_results_df[STADIA_results_df['timepoint']==tp]
        tp_emb = tp_STADIA_results.iloc[:, :35].to_numpy()
        tp_kmeans = KMeans(n_clusters=len(set(tp_STADIA_results[label_key])), random_state=0).fit(tp_emb)
        tp_STADIA_results['kmeans'] = tp_kmeans.labels_.astype(str)
        STADIA_results_list.append(tp_STADIA_results)
        kmeans_dict[tp] = tp_kmeans
    STADIA_results = pd.concat(STADIA_results_list)
    STADIA_results.rename(columns=rename_dict, inplace=True)
    analyze_domain_metrics(STADIA_results, STADIA_results_dir, "STADIA")
    try:
        compute_OT_heatmap(kmeans_dict, STADIA_results_dir, method='STADIA')
    except RuntimeError as e:
        print(f"Runtime Error: {e}")
    STADIA_adata = anndata.AnnData(STADIA_results.iloc[:, :35].to_numpy(), obs=STADIA_results)
    STADIA_adata.write_h5ad(STADIA_results_dir+os.sep+'prediction_adata.h5ad')
    return STADIA_adata

# --- analyze latent space of each method
def analyze_spatialPCA_latent(result_dir, result_prefix):
    # --- SpatialPCA latent --- 
    spatialPCA_results_dir = result_dir+os.sep+"explore_spatialPCA"+os.sep+result_prefix
    spatialPCA_results = pd.read_csv(spatialPCA_results_dir+os.sep+'spatialPCA_results.csv', index_col=0)
    spatialPCA_latent = spatialPCA_results.iloc[:, :20]
    # analyze latent
    analyze_latent_metrics(spatialPCA_latent, spatialPCA_results[['timepoint']], spatialPCA_results_dir, "spatialPCA")

def analyze_spaDOT_latent(result_dir, result_prefix, method='spaDOT'):
    # --- spaDOT latent ---
    spaDOT_results_dir = result_dir+os.sep+"explore_"+method+os.sep+result_prefix
    spaDOT_adata = anndata.read_h5ad(spaDOT_results_dir+os.sep+'final_latent.h5ad')
    spaDOT_latent = pd.DataFrame(spaDOT_adata.X)
    spaDOT_latent.index = spaDOT_adata.obs.index
    # analyze latent
    analyze_latent_metrics(spaDOT_latent, spaDOT_adata.obs[['timepoint']], spaDOT_results_dir, "spaDOT")

def analyze_GraphST_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- GraphST latent ---
    GraphST_results_dir = result_dir+os.sep+"explore_GraphST"+os.sep+result_prefix
    GraphST_adata = anndata.read_h5ad(GraphST_results_dir+os.sep+'GraphST_adata.h5ad')
    GraphST_latent = pd.DataFrame(GraphST_adata.obsm['emb'])
    GraphST_latent.index = GraphST_adata.obs.index
    # analyze latent
    analyze_latent_metrics(GraphST_latent, GraphST_adata.obs[['timepoint']], GraphST_results_dir, "GraphST")


def analyze_SEDR_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- SEDR latent ---
    SEDR_results_dir = result_dir+os.sep+"explore_SEDR"+os.sep+result_prefix
    SEDR_adata = anndata.read_h5ad(SEDR_results_dir+os.sep+'SEDR_adata.h5ad')
    SEDR_latent = pd.DataFrame(SEDR_adata.obsm['SEDR.Harmony'])
    SEDR_latent.index = SEDR_adata.obs.index
    # analyze latent
    analyze_latent_metrics(SEDR_latent, SEDR_adata.obs[['timepoint']], SEDR_results_dir, "SEDR")

def analyze_spaGCN_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- spaGCN latent ---
    spaGCN_results_dir = result_dir+os.sep+"explore_spaGCN"+os.sep+result_prefix
    spaGCN_adata = anndata.read_h5ad(spaGCN_results_dir+os.sep+'spaGCN.h5ad')
    spaGCN_latent = pd.DataFrame(spaGCN_adata.obsm['spaGCN_embed'])
    spaGCN_latent.index = spaGCN_adata.obs.index
    # analyze latent
    analyze_latent_metrics(spaGCN_latent, spaGCN_adata.obs[['timepoint']], spaGCN_results_dir, "spaGCN")

def analyze_spaVAE_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- spaVAE latent ---
    spaVAE_results_dir = result_dir+os.sep+"explore_spaVAE"+os.sep+result_prefix
    metadata_files = [f for f in os.listdir(spaVAE_results_dir) if f.endswith("_metadata.csv")]
    assert len(metadata_files) == 1
    metadata = pd.read_csv(spaVAE_results_dir+os.sep+metadata_files[0], index_col=0)
    latent = np.loadtxt(spaVAE_results_dir+os.sep+'final_latent.txt', delimiter=',')
    spaVAE_latent = pd.DataFrame(latent)
    spaVAE_latent.index = metadata.index
    # analyze latent
    analyze_latent_metrics(spaVAE_latent, metadata[['timepoint']], spaVAE_results_dir, "spaVAE")

def analyze_STAligner_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- STAligner latent ---
    STAligner_results_dir = result_dir+os.sep+"explore_STAligner"+os.sep+result_prefix
    STAligner_adata = anndata.read_h5ad(STAligner_results_dir+os.sep+'STAligner_adata.h5ad')
    STAligner_latent = pd.DataFrame(STAligner_adata.obsm['STAligner'])
    STAligner_latent.index = STAligner_adata.obs.index
    # analyze latent
    analyze_latent_metrics(STAligner_latent, STAligner_adata.obs[['timepoint']], STAligner_results_dir, "STAligner")

def analyze_PRECAST_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- PRECAST latent ---
    PRECAST_results_dir = result_dir+os.sep+"explore_PRECAST"+os.sep+result_prefix
    PRECAST_results = pd.read_csv(PRECAST_results_dir+os.sep+'PRECAST_results.csv', index_col=0)
    PRECAST_latent = PRECAST_results.iloc[:, :15]
    # analyze latent
    analyze_latent_metrics(PRECAST_latent, PRECAST_results[['timepoint']], PRECAST_results_dir, "PRECAST")

def analyze_spatiAlign_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- spatiAlign latent ---
    spatiAlign_results_dir = result_dir+os.sep+"explore_spatiAlign"+os.sep+result_prefix
    # list all res/correct_data*.h5ad
    correct_adata_files = glob.glob(spatiAlign_results_dir+os.sep+"res/correct_data*.h5ad")
    correct_adata_list = list()
    for correc_adata_file in correct_adata_files:
        correct_adata = sc.read_h5ad(correc_adata_file)
        correct_adata.var['genes'] = correct_adata.var_names.tolist()
        correct_adata_list.append(correct_adata)
    merge_data = anndata.AnnData.concatenate(*correct_adata_list)
    spatiAlign_latent = pd.DataFrame(merge_data.obsm['correct'])
    spatiAlign_latent.index = merge_data.obs.index
    # analyze latent
    analyze_latent_metrics(spatiAlign_latent, merge_data.obs[['timepoint']], spatiAlign_results_dir, "spatiAlign")

def analyze_SPIRAL_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- SPIRAL latent ---
    SPIRAL_results_dir = result_dir+os.sep+"explore_SPIRAL"+os.sep+result_prefix
    SPIRAL_adata = anndata.read_h5ad(SPIRAL_results_dir+os.sep+'SPIRAL_corrected.h5ad')
    SPIRAL_latent = pd.DataFrame(SPIRAL_adata.obsm['spiral'])
    SPIRAL_latent.index = SPIRAL_adata.obs.index
    # analyze latent
    analyze_latent_metrics(SPIRAL_latent, SPIRAL_adata.obs[['timepoint']], SPIRAL_results_dir, "SPIRAL")

def analyze_STADIA_latent(result_dir, result_prefix, color_key=['timepoint', 'annotation']):
    # --- STADIA latent ---
    STADIA_results_dir = result_dir+os.sep+"explore_STADIA"+os.sep+result_prefix
    STADIA_results = pd.read_csv(STADIA_results_dir+os.sep+'STADIA_results.csv', index_col=0)
    STADIA_latent = STADIA_results.iloc[:, :35]
    # analyze latent
    analyze_latent_metrics(STADIA_latent, STADIA_results[['timepoint']], STADIA_results_dir, "STADIA")
