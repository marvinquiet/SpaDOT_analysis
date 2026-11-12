## conda activate /net/mulan/home/wenjinma/envs/spatialATAC
Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")

.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)
library(anndata)
library(PRECAST)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)


base_dir = '/net/mulan/home/wenjinma'
project_dir = file.path(base_dir, "projects/spatialATAC-explore")
work_dir = file.path(project_dir, 'explore_PRECAST')
setwd(work_dir)
set.seed(2024)
data_dir = file.path(base_dir, 'data/spatial_Artista')
result_dir = file.path(work_dir, 'Axolotl_results')
dir.create(result_dir, showWarnings = F)


#adata = read_h5ad(file.path(data_dir, 'Development_curated.h5ad'))
adata = read_h5ad(file.path(data_dir, 'Development_addCP.h5ad'))
Seurat_obj = CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = adata$obs)
Seurat_obj$pixel_x = adata$obsm[['spatial']][, 1]
Seurat_obj$pixel_y = -adata$obsm[['spatial']][, 2]
Seurat_obj$row = Seurat_obj$pixel_x
Seurat_obj$col = Seurat_obj$pixel_y
timepoints = c('t0', 't1', 't2')
Seurat_obj = subset(Seurat_obj, subset=timepoint %in% timepoints)
Seurat_obj$annotation = Seurat_obj$region_annotation

seuList = list()
seuMetaList = list()
for (tp in unique(Seurat_obj$timepoint)) {
    tp_obj = subset(Seurat_obj, subset=timepoint==tp)
    seuList[[tp]] = tp_obj
}

# run PRECAST
preobj = CreatePRECASTObject(seuList = seuList, selectGenesMethod = "HVGs", gene.number = 2000) 
PRECASTObj = AddAdjList(preobj, platform = "Visium")
PRECASTObj = AddParSetting(PRECASTObj, Sigma_equal = TRUE, coreNum = 1, maxIter = 30, verbose = TRUE)
PRECASTObj = PRECAST(PRECASTObj, K = 10) # set 5 for now because in total there will be 5 clusters
## backup the fitting results in resList
resList = PRECASTObj@resList
PRECASTObj = SelectModel(PRECASTObj)
seuInt = IntegrateSpaData(PRECASTObj, species='Unknown')
seuInt = AddTSNE(seuInt, n_comp = 2) 
seuInt@meta.data = Seurat_obj@meta.data[rownames(seuInt@meta.data), ]

# perform KMeans clustering on the embeddings
sample_res_all = NULL
ARI_res = c()
for (tp in unique(seuInt$timepoint)) {
    tp_seu = subset(seuInt, subset=timepoint == tp)
    tp_seu_embed = Embeddings(tp_seu, reduction='PRECAST')
    tp_kmeans = kmeans(tp_seu_embed, centers=length(unique(tp_seu$annotation)))
    tp_seu$cluster = tp_kmeans$cluster
    tp_ari = ARI(tp_seu$annotation, tp_seu$cluster)
    cat(tp, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # plot embeddings
    cols_cluster = chooseColors(palettes_name = 'Nature 10', n_colors = length(unique(tp_seu$annotation)), plot_colors = F)
    p1 = dimPlot(tp_seu, item='annotation', point_size = 1,  font_family='serif')
    p2 = dimPlot(tp_seu, item='cluster', font_family='serif', cols=cols_cluster) # Times New Roman
    drawFigs(list(p1, p2), common.legend=FALSE, align='hv', layout.dim=c(2,1)) 
    ggsave(file.path(result_dir, paste0(tp, '_tSNE.png')))
    # save sample res
    sample_res = cbind(tp_seu_embed, tp_seu@meta.data)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'PRECAST_results.csv'))
names(ARI_res) = unique(seuInt$timepoint)
write.csv(data.frame(ARI_res), file.path(result_dir, 'PRECAST_ARI.csv'))

