# conda activate /media/sda/wenjinma/envs/GraphST
# --- scripts for generate edges: https://github.com/guott15/SPIRAL/blob/main/Demo/GenerateEdges.ipynb
# --- scripts for preprocessing: https://github.com/guott15/SPIRAL/blob/main/Demo/preprocess.ipynb
import os, argparse
import scanpy as sc
import pandas as pd
import numpy as np
import sklearn.neighbors
import anndata

import torch
from spiral.main import SPIRAL_integration
from spiral.layers import *
from spiral.utils import *
os.environ['CUDA_VISIBLE_DEVICES'] = '0'


def Cal_Spatial_Net(adata, rad_cutoff=None, k_cutoff=None, model='Radius', verbose=True):
    """\
    Construct the spatial neighbor networks.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    rad_cutoff
        radius cutoff when model='Radius'
    k_cutoff
        The number of nearest neighbors when model='KNN'
    model
        The network construction model. When model=='Radius', the spot is connected to spots whose distance is less than rad_cutoff. When model=='KNN', the spot is connected to its first k_cutoff nearest neighbors.
    
    Returns
    -------
    The spatial networks are saved in adata.uns['Spatial_Net']
    """
    assert(model in ['Radius', 'KNN'])
    if verbose:
        print('------Calculating spatial graph...')
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
#     coor.columns = ['imagerow', 'imagecol']
    if model == 'Radius':
        nbrs = sklearn.neighbors.NearestNeighbors(radius=rad_cutoff).fit(coor)
        distances, indices = nbrs.radius_neighbors(coor, return_distance=True)
        KNN_list = []
        for it in range(indices.shape[0]):
            KNN_list.append(pd.DataFrame(zip([it]*indices[it].shape[0], indices[it], distances[it])))
    if model == 'KNN':
        nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=k_cutoff+1).fit(coor)
        distances, indices = nbrs.kneighbors(coor)
        KNN_list = []
        for it in range(indices.shape[0]):
            KNN_list.append(pd.DataFrame(zip([it]*indices.shape[1],indices[it,:], distances[it,:])))
    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']
    Spatial_Net = KNN_df.copy()
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance']>0,]
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index), ))
    Spatial_Net['Cell1'] = Spatial_Net['Cell1'].map(id_cell_trans)
    Spatial_Net['Cell2'] = Spatial_Net['Cell2'].map(id_cell_trans)
    if verbose:
        print('The graph contains %d edges, %d cells.' %(Spatial_Net.shape[0], adata.n_obs))
        print('%.4f neighbors per cell on average.' %(Spatial_Net.shape[0]/adata.n_obs))
    adata.uns['Spatial_Net'] = Spatial_Net

def Preprocess(adata, result_dir, knn=10, rad=10):
    Batch_list = []
    timepoints = list(set(adata.obs['timepoint']))
    timepoints.sort()
    # preprocess data per timepoint
    for tp in timepoints:
        tp_adata = adata[adata.obs['timepoint']==tp].copy()
        sc.pp.normalize_total(tp_adata, target_sum=1e4)
        sc.pp.log1p(tp_adata)
        sc.pp.highly_variable_genes(tp_adata, flavor="seurat_v3", n_top_genes=5000)
        tp_adata = tp_adata[:, tp_adata.var['highly_variable']]
        Batch_list.append(tp_adata)
    # concatenate data together
    adata_concat = anndata.concat(Batch_list, label="slice_name", keys=timepoints)
    vf=np.array(adata_concat.var.index)
    for i in np.arange(len(timepoints)):
        adata1=adata_concat[adata_concat.obs['slice_name']==timepoints[i],:]
        cells=adata1.obs_names
        mat=pd.DataFrame(adata1.X,index=cells,columns=vf)
        coord=pd.DataFrame(adata1.obsm['spatial'],index=cells,columns=['x','y'])
        meta=adata1.obs[['annotation','slice_name']]
        meta.columns=['celltype','batch']
        meta.index=cells
        mat.to_csv(result_dir+os.sep+timepoints[i]+"_mat.csv")
        meta.to_csv(result_dir+os.sep+timepoints[i]+"_meta.csv")
        coord.to_csv(result_dir+os.sep+timepoints[i]+"_coord.csv")
    # calculate edeges
    for tp in timepoints:
        feat=pd.read_csv(result_dir+os.sep+tp+"_mat.csv",header=0,index_col=0,sep=',')
        coord=pd.read_csv(result_dir+os.sep+tp+"_coord.csv",header=0,index_col=0,sep=',')
        coord.columns=['x','y']
        tp_adata = sc.AnnData(feat)
        tp_adata.var_names_make_unique()
        tp_adata.X=sp.csr_matrix(tp_adata.X)
        tp_adata.obsm["spatial"] = coord.loc[:,['x','y']].to_numpy()
        Cal_Spatial_Net(tp_adata, rad_cutoff=rad, k_cutoff=knn, model='KNN', verbose=True)
        if 'highly_variable' in tp_adata.var.columns:
            adata_Vars =  tp_adata[:, tp_adata.var['highly_variable']]
        else:
            adata_Vars = tp_adata
        features = pd.DataFrame(adata_Vars.X.toarray()[:, ], index=adata_Vars.obs.index, columns=adata_Vars.var.index)
        cells = np.array(features.index)
        # cells_id_tran = dict(zip(cells, range(cells.shape[0])))
        if 'Spatial_Net' not in tp_adata.uns.keys():
            raise ValueError("Spatial_Net is not existed! Run Cal_Spatial_Net first!")
        Spatial_Net = tp_adata.uns['Spatial_Net']
        G_df = Spatial_Net.copy()
        np.savetxt(result_dir+os.sep+tp+"_edge_KNN_"+str(knn)+".csv",G_df.values[:,:2],fmt='%s')

def train_SPIRAL(result_dir, timepoints, knn=10):
    feat_file, edge_file, meta_file, coord_file = [], [], [], []
    for i in range(len(timepoints)):
        feat_file.append(result_dir+os.sep+str(timepoints[i])+"_mat.csv")
        edge_file.append(result_dir+os.sep+str(timepoints[i])+"_edge_KNN_"+str(knn)+".csv")
        meta_file.append(result_dir+os.sep+str(timepoints[i])+"_meta.csv")
        coord_file.append(result_dir+os.sep+str(timepoints[i])+"_coord.csv")
    N=pd.read_csv(feat_file[0],header=0,index_col=0).shape[1]
    if (len(timepoints)==2):
        M=1
    else:
        M=len(timepoints)

    N_WALKS=knn
    WALK_LEN=1
    N_WALK_LEN=knn
    NUM_NEG=knn
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=0, help='The seed of initialization.')
    parser.add_argument('--AEdims', type=list, default=[N,[512],32], help='Dim of encoder.')
    parser.add_argument('--AEdimsR', type=list, default=[32,[512],N], help='Dim of decoder.')
    parser.add_argument('--GSdims', type=list, default=[512,32], help='Dim of GraphSAGE.')
    parser.add_argument('--zdim', type=int, default=32, help='Dim of embedding.')
    parser.add_argument('--znoise_dim', type=int, default=4, help='Dim of noise embedding.')
    parser.add_argument('--CLdims', type=list, default=[4,[],M], help='Dim of classifier.')
    parser.add_argument('--DIdims', type=list, default=[28,[32,16],M], help='Dim of discriminator.')
    parser.add_argument('--beta', type=float, default=1.0, help='weight of GraphSAGE.')
    parser.add_argument('--agg_class', type=str, default=MeanAggregator, help='Function of aggregator.')
    parser.add_argument('--num_samples', type=str, default=knn, help='number of neighbors to sample.')
    parser.add_argument('--N_WALKS', type=int, default=N_WALKS, help='number of walks of random work for postive pairs.')
    parser.add_argument('--WALK_LEN', type=int, default=WALK_LEN, help='walk length of random work for postive pairs.')
    parser.add_argument('--N_WALK_LEN', type=int, default=N_WALK_LEN, help='number of walks of random work for negative pairs.')
    parser.add_argument('--NUM_NEG', type=int, default=NUM_NEG, help='number of negative pairs.')
    parser.add_argument('--epochs', type=int, default=100, help='Number of epochs to train.')
    parser.add_argument('--batch_size', type=int, default=1024, help='Size of batches to train.') ####512 for withon donor;1024 for across donor###
    parser.add_argument('--lr', type=float, default=1e-3, help='Initial learning rate.')
    parser.add_argument('--weight_decay', type=float, default=5e-4, help='Weight decay.')
    parser.add_argument('--alpha1', type=float, default=N, help='Weight of decoder loss.')
    parser.add_argument('--alpha2', type=float, default=1, help='Weight of GraphSAGE loss.')
    parser.add_argument('--alpha3', type=float, default=1, help='Weight of classifier loss.')
    parser.add_argument('--alpha4', type=float, default=1, help='Weight of discriminator loss.')
    parser.add_argument('--lamda', type=float, default=1, help='Weight of GRL.')
    parser.add_argument('--Q', type=float, default=10, help='Weight negative loss for sage losss.')
    params, _ =parser.parse_known_args()

    # train SPIRAL
    SPII=SPIRAL_integration(params,feat_file,edge_file,meta_file)
    SPII.train()
    SPII.model.eval()
    all_idx=np.arange(SPII.feat.shape[0])
    all_layer,all_mapping=layer_map(all_idx.tolist(),SPII.adj,len(SPII.params.GSdims))
    all_rows=SPII.adj.tolil().rows[all_layer[0]]
    all_feature=torch.Tensor(SPII.feat.iloc[all_layer[0],:].values).float().cuda()
    all_embed,ae_out,clas_out,disc_out=SPII.model(all_feature,all_layer,all_mapping,all_rows,SPII.params.lamda,SPII.de_act,SPII.cl_act)
    [ae_embed,gs_embed,embed]=all_embed
    embed=embed.cpu().detach()
    names=['GTT_'+str(i) for i in range(embed.shape[1])]
    embed1=pd.DataFrame(np.array(embed),index=SPII.feat.index,columns=names)
    if not os.path.exists(result_dir+os.sep+"gtt_output/"):
        os.makedirs(result_dir+os.sep+"gtt_output/")
    embed_file=result_dir+os.sep+"gtt_output/SPIRAL_embed_"+str(SPII.params.batch_size)+".csv"
    embed1.to_csv(embed_file)
    # obtain embeddings
    ann=anndata.AnnData(SPII.feat)
    ann.obsm['spiral']=embed1.iloc[:,SPII.params.znoise_dim:].values
    # ann.obs = adata.obs.loc[ann.obs_names, :]
    sc.pp.neighbors(ann,use_rep='spiral')
    ann.write_h5ad(result_dir+os.sep+"SPIRAL.h5ad")