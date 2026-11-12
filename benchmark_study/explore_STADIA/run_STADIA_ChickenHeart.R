## conda activate ~/envs/spatialATAC
Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)
library(anndata)

library(stadia)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)
library(optparse)
library(Matrix)

## on Mulan
base_dir = '/net/mulan/home/wenjinma'
project_dir = file.path(base_dir, "projects/spatialATAC-explore")
work_dir = file.path(project_dir, 'explore_STADIA')
setwd(work_dir)
set.seed(2024)
data_dir = file.path(base_dir, 'data/spatial_ChickenHeart')
result_dir = file.path(work_dir, 'ChickenHeart_results')
dir.create(result_dir, showWarnings = F)

adata = read_h5ad(file.path(data_dir, 'processed_ChickenHeart.h5ad'))
Seurat_obj = CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = adata$obs)
Seurat_obj$row = Seurat_obj$pixel_x
Seurat_obj$col = Seurat_obj$pixel_y
Seurat_obj$annotation = Seurat_obj$region

seuList = list()
seuMetaList = list()
for (tp in unique(Seurat_obj$timepoint)) {
    sub_obj = subset(Seurat_obj, subset=timepoint==tp)
    sub_assay_v3 = CreateAssayObject(counts=sub_obj[['RNA']]$'counts')
    sub_obj[['RNA']] = sub_assay_v3
    seuList[[tp]] = sub_obj
    seuMetaList[[tp]] = sub_obj@meta.data
}
seuMeta = do.call(rbind, seuMetaList)

## set parameters for STADIA model
etas = 0.15
K = 10
hyper = HyperParameters(seuList, eta = etas)
stadia = stadia(seuList, hyper, n_cluster = K, min.features = 0, min.spots = 20, # do not remove cells based on feature expression
                platform = "others", em.maxiter = 30)
save(stadia, file=file.path(result_dir, 'stadia.RData'))

# perform KMeans clustering on the embeddings
sample_res_all = NULL
ARI_res = c()
for (tp in unique(seuMeta$timepoint)) {
    # rename stadia factors
    stadia_factors = stadia$factors
    colnames(stadia_factors) = paste0(tp, '.', colnames(stadia_factors))
    tp_meta = seuMeta[seuMeta$timepoint == tp, ]
    common_cells = intersect(rownames(tp_meta), colnames(stadia_factors))
    tp_embed = t(stadia_factors)[common_cells, ]
    tp_meta = tp_meta[common_cells, ]
    tp_kmeans = kmeans(tp_embed, centers=length(unique(tp_meta$annotation)))
    tp_meta$cluster = tp_kmeans$cluster
    tp_ari = ARI(tp_meta$annotation, tp_meta$cluster)
    cat(tp, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # visualizing results
    g = ggplot(tp_meta, aes(row, col)) + 
        geom_point(size = .3, shape = 21, color = "#00000000", aes(fill = factor(cluster))) +
        labs(subtitle = paste0("ARI:", round(tp_ari, 2)))
    ggsave(file.path(result_dir, paste0(tp, '_results.png')))
    # save embeddings
    sample_res = cbind(tp_embed, tp_meta)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'STADIA_results.csv'))
names(ARI_res) = unique(seuMeta$timepoint)
write.csv(data.frame(ARI_res), file.path(result_dir, 'STADIA_ARI.csv'))

