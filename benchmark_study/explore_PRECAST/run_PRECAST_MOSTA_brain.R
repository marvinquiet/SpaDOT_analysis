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
saveRDS(seuInt, file.path(result_dir, 'PRECAST_res.RDS'))

