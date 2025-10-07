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
data_dir = file.path(base_dir, 'data/spatial_MOSTA')
result_dir = file.path(work_dir, 'MOSTA_brain_results')
dir.create(result_dir, showWarnings = F)

adata = read_h5ad(file.path(data_dir, 'adata_brain_clean.h5ad'))
n_cells = 20000
subset_idx = sample(nrow(adata), n_cells)
adata = adata[subset_idx,]
adata$obs$annotation = 'Brain'
Seurat_obj = CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = adata$obs)
Seurat_obj$pixel_x = adata$obsm$spatial[,1]
Seurat_obj$pixel_y = adata$obsm$spatial[,2]
Seurat_obj$row = Seurat_obj$pixel_x
Seurat_obj$col = Seurat_obj$pixel_y
Seurat_obj$annotation = 'Brain'

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
saveRDS(stadia, file=file.path(result_dir, 'stadia.RDS'))
