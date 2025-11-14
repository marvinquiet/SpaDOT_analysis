Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(dplyr)
library(Rtsne)
library(umap)

# adata = read_h5ad('/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop/1993/prop_0.9/adata_downsampled.h5ad')
# adata = read_h5ad('/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_DEG5000/1993/prop_0.9/adata_downsampled.h5ad')
adata = read_h5ad('/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_mixedNPCs_DEG500/1993/prop_0.9/adata_downsampled.h5ad')
obj = CreateSeuratObject(t(adata$X), meta.data=adata$obs, assay='RNA')
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = ScaleData(obj, features = rownames(obj))
Idents(obj) = obj$celltype
markers = FindAllMarkers(obj, only.pos = TRUE)
cluster_markers = markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# I have checked the markers, which are mostly the same as I simulated them (only Neurons / GABA has two markers in gene id > 2000)
