# conda activate /media/sda/wenjinma/envs/Siyu_CCC/
import os
import commot as ct
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# load data
import anndata

base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_scRNAseq_GSE214989"

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
adata.obs.loc[common_barcodes, 'x'] = loc_metadata.loc[common_barcodes, 'x'].values
adata.obs.loc[common_barcodes, 'y'] = loc_metadata.loc[common_barcodes, 'y'].values
adata.obs.loc[common_barcodes, 'pixel_x'] = loc_metadata.loc[common_barcodes, 'pixel_x'].values
adata.obs.loc[common_barcodes, 'pixel_y'] = loc_metadata.loc[common_barcodes, 'pixel_y'].values
adata.obs['annotation'] = adata.obs['seurat_clusters_harmony']
adata.obs['annotation'].replace({3: 'CNS', 2: 'CNS', 4: 'CNS',
                                 1: 'Men/PNS', 9: 'Men/PNS', 11: 'Men/PNS',
                                 0: 'Mesenchyme', 5: 'Mesenchyme', 6: 'Mesenchyme', 8: 'Mesenchyme',
                                 7: 'Liver', 10: 'Heart'}, inplace=True)
adata.obs['annotation'] = adata.obs['annotation'].astype('category')
adata.obs['timepoint'] = [section_dict[_] for _ in adata.obs['section']]
adata.obs['timepoint_numeric'] = adata.obs['timepoint'].astype('category').cat.codes
loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
adata.obsm['spatial'] = loc
adata.obs['data'] = adata.obs['timepoint']

result_dir = "/media/sda/wenjinma/projects/spatialATAC-explore/explore_GATVAE-SVGP-VAE/ATAC-Visium_noremoval_results"
latent_adata = anndata.read_h5ad(result_dir+os.sep+'prediction_adata.h5ad')
common_cells = list(set(latent_adata.obs_names).intersection(set(adata.obs_names)))
adata = adata[common_cells].copy()
latent_adata = latent_adata[common_cells].copy()
adata.obs['pred_labels'] = latent_adata.obs['pred_labels'].values

# preprocess data for each time point
df_cellchat=ct.pp.ligand_receptor_database(database='CellChat', species='mouse')
for tp in adata.obs['timepoint'].unique():
    # tp_adata = adata[adata.obs['timepoint'] == tp].copy()
    # sc.pp.normalize_total(tp_adata, inplace=True)
    # sc.pp.log1p(tp_adata)
    # ct.tl.spatial_communication(tp_adata,
    # database_name='cellchat', df_ligrec=df_cellchat, dis_thr=1000, heteromeric=True, pathway_sum=True)
    # tp_adata.obs['x'] = tp_adata.obs['x'].astype(int)
    # tp_adata.obs['y'] = tp_adata.obs['y'].astype(int)
    # tp_adata.obs['pixel_x'] = tp_adata.obs['pixel_x'].astype(int)
    # tp_adata.obs['pixel_y'] = tp_adata.obs['pixel_y'].astype(int)
    # tp_adata.write(result_dir+os.sep+'commot_adata_'+tp+'.h5ad')
    tp_adata = anndata.read_h5ad(result_dir+os.sep+'commot_adata_'+tp+'.h5ad')
    # total pathway
    ct.tl.communication_direction(tp_adata, database_name='cellchat', k=5)
    ct.pl.plot_cell_communication(tp_adata, database_name='cellchat', plot_method='grid', background_legend=True,
    scale=0.000001, ndsize=8, grid_density=0.4, summary='sender', clustering='pred_labels', 
    normalize_v = True, normalize_v_quantile=0.995)
    plt.savefig(result_dir+os.sep+'total_communication_sender_'+tp+'.png', dpi=300)

    ct.pl.plot_cell_communication(tp_adata, database_name='cellchat', plot_method='grid', background_legend=True,
    scale=0.000001, ndsize=8, grid_density=0.4, summary='receiver', clustering='pred_labels', 
    normalize_v = True, normalize_v_quantile=0.995)
    plt.savefig(result_dir+os.sep+'total_communication_receiver_'+tp+'.png', dpi=300)

    # WNT can be a very important signaling pathway
    pathway_names = ['WNT', 'BMP', 'TGFb', 'FGF', 'PDGF', 'VEGF', 'ncWNT']
    for pathway in pathway_names:
        ct.tl.communication_direction(tp_adata, database_name='cellchat', pathway_name=pathway, k=5)
        ct.pl.plot_cell_communication(tp_adata, database_name='cellchat', pathway_name=pathway, 
                                    plot_method='grid', background_legend=True,
            scale=0.0000005, ndsize=8, grid_density=0.4, summary='sender', clustering='pred_labels', normalize_v = True, normalize_v_quantile=0.995)
        plt.savefig(result_dir+os.sep+pathway+'_communication_sender_'+tp+'.png', dpi=300)
        ct.pl.plot_cell_communication(tp_adata, database_name='cellchat', pathway_name=pathway, 
                                    plot_method='grid', background_legend=True,
            scale=0.0000005, ndsize=8, grid_density=0.4, summary='receiver', clustering='pred_labels', normalize_v = True, normalize_v_quantile=0.995)
        plt.savefig(result_dir+os.sep+pathway+'_communication_receiver_'+tp+'.png', dpi=300)



    ct.tl.communication_direction(tp_adata, database_name='cellchat', pathway_name='TGFb', k=5)
    ct.pl.plot_cell_communication(tp_adata, database_name='cellchat', pathway_name='TGFb', 
                                  plot_method='grid', background_legend=True,
        scale=0.0000005, ndsize=8, grid_density=0.4, summary='sender', clustering='pred_labels', normalize_v = True, normalize_v_quantile=0.995)
    plt.savefig(result_dir+os.sep+'TGFb_communication_sender_'+tp+'.png', dpi=300)
    ct.pl.plot_cell_communication(tp_adata, database_name='cellchat', pathway_name='TGFb', 
                                  plot_method='grid', background_legend=True,
        scale=0.0000005, ndsize=8, grid_density=0.4, summary='receiver', clustering='pred_labels', normalize_v = True, normalize_v_quantile=0.995)
    plt.savefig(result_dir+os.sep+'TGFb_communication_receiver_'+tp+'.png', dpi=300)

    ct.tl.cluster_communication(tp_adata, database_name='cellchat', clustering='pred_labels', n_permutations=1000)
    ct.tl.cluster_position(tp_adata, clustering='pred_labels')
    ct.pl.plot_cluster_communication_network(tp_adata, uns_names=['commot_cluster-pred_labels-cellchat-total-total'], clustering='pred_labels',
    nx_node_pos='cluster', nx_pos_idx=np.array([0, 1]), nx_bg_ndsize=0.25, 
    filename=result_dir+os.sep+'total_communication_cluster_'+tp+'.pdf', nx_node_cmap='Light24')


    ct.tl.cluster_communication(tp_adata, database_name='cellchat', pathway_name='TGFb', clustering='pred_labels', n_permutations=1000)
    ct.pl.plot_cluster_communication_network(tp_adata, uns_names=['commot_cluster-pred_labels-cellchat-total-total'], clustering='pred_labels',
    nx_node_pos='cluster', nx_pos_idx=np.array([0, 1]), nx_bg_ndsize=0.25, 
    filename=result_dir+os.sep+'TGFb_communication_cluster_'+tp+'.pdf', nx_node_cmap='Light24')

