import os
import anndata
prefix_dir = "/net/mulan/home/wenjinma"
data_dir = prefix_dir+os.sep+'data/spatial_Artista'
adata = anndata.read_h5ad(data_dir+os.sep+'Development.h5ad')

# turn adata.X to adata counts
adata.X = adata.layers['counts']
adata.uns = {}
del adata.layers
adata.var_names = adata.var['Gene'].values
adata.obs['status'] = [_.split('_')[0] for _ in  adata.obs['Batch']]
adata.obs['timepoint'] = adata.obs['status'].replace({'Stage44': 't0',
        'Stage54': 't1', 'Stage57': 't2', 'Injury': 't3', 'Batch1': 't4', 'Meta': 't5'})
adata.obs_names = adata.obs_names+'_'+adata.obs['status']
adata.write_h5ad(data_dir+os.sep+'Development_reformatted.h5ad')

