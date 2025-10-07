## conda activate /media/sda/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(SpatialPCA)
library(SeuratDisk)
library(Seurat)
#library(Signac)
library(aricode)
library(ggplot2)
library(rhdf5)

prefix_dir = "/net/mulan/home/wenjinma"
#prefix_dir = "/media/sda/wenjinma"

data_dir = file.path(prefix_dir, "data/spatial_MOSTA")
result_dir = file.path(prefix_dir, "projects/spatialATAC-explore/results/SpatialPCA_MOSTA_E12.5/")
dir.create(result_dir, showWarnings = F)

# load data
Seurat_obj = LoadH5Seurat(file.path(data_dir, 'Mouse_embryo_all_stage.h5seurat'))
annotation_counts = table(Seurat_obj$annotation)
selected_tissues = names(annotation_counts)[annotation_counts > 10000]
selected_idx = Seurat_obj$annotation %in% selected_tissues & Seurat_obj$timepoint %in% c('E12.5')
filtered_Seurat_obj = Seurat_obj[, selected_idx]
set.seed(1234)
sampled_idx = sample(ncol(filtered_Seurat_obj), size=5000)
sampled_obj = filtered_Seurat_obj[, sampled_idx]
sampled_obj[["percent.mt"]] = PercentageFeatureSet(sampled_obj, pattern = "^MT-")

## create SpatialPCA object
spatialPCA_obj = CreateSpatialPCAObject(counts=sampled_obj@assays$count@counts, location=Embeddings(sampled_obj, reduction='spatial'),
                                        project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,
                                        gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
bandwidth_common=0.1
#spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidth.set.by.user=bandwidth_common)
spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
spatialPCA_obj = SpatialPCA_EstimateLoading(spatialPCA_obj,fast=TRUE,SpatialPCnum=20)
spatialPCA_obj = SpatialPCA_SpatialPCs(spatialPCA_obj, fast=TRUE)
colnames(spatialPCA_obj@SpatialPCs) = rownames(spatialPCA_obj@location)

cluster_annotations = sampled_obj@meta.data$annotation
spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
                                              latent_dat=spatialPCA_obj@SpatialPCs,
                                              knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2])))
plot_cluster(location=Embeddings(sampled_obj, reduction='spatial'),
             clusterlabel=spatialPCA_clusterlabel,
             pointsize=1.5, text_size=20,
             title_in='E12.5', color_in=rainbow(length(unique(cluster_annotations))),
             legend='bottom')
ggsave(file.path(result_dir, 'E12.5_SpatialPCA.png'))
cat('ARI:', ARI(spatialPCA_clusterlabel, cluster_annotations), '\n')

sampled_obj = NormalizeData(sampled_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sampled_obj = FindVariableFeatures(sampled_obj, selection.method = "vst", nfeatures = 2000)
sampled_obj = ScaleData(sampled_obj, features = rownames(sampled_obj))
# add spatialPCA pcs to sampled_obj
sampled_obj[['spatialpca']] = CreateDimReducObject(embeddings = t(spatialPCA_obj@SpatialPCs), key = 'spatialPCA_', assay = 'RNA')
sampled_obj = RunUMAP(sampled_obj, reduction = "spatialpca", reduction.key='spatialumap_', dims = 1:20)
DimPlot(sampled_obj, reduction = "umap", group.by='annotation')
ggsave(file.path(result_dir, 'spatialPCA_annotation_umap.png'))

