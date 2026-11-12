import os
import h5py
import anndata
import numpy as np
import pandas as pd
base_dir = "/media/sda/wenjinma"
data_dir = base_dir+os.sep+"data/spatial_scRNAseq_GSE214989"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore/explore_spaVAE/ATAC-Visium_noremoval_results"
os.makedirs(result_dir, exist_ok=True)

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
adata.obs.to_csv(result_dir+os.sep+'combined_ATAC-Visium_metadata.csv')
print("Data size:", adata.shape)

loc = np.array(adata.obs[['pixel_x', 'pixel_y']]).astype('float64')
batch_info = np.array(adata.obs['timepoint_numeric']).astype('int')
batch_mat = np.zeros((batch_info.size, batch_info.max()+1))
batch_mat[np.arange(batch_info.size), batch_info] = 1
 
with h5py.File(result_dir+os.sep+'combined_ATAC_visium.h5', 'w') as file:
    file.create_dataset('X', data=adata.X.toarray())
    file.create_dataset('pos', data=loc)
    file.create_dataset('batch', data=batch_mat)
