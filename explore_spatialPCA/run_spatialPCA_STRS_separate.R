## conda activate /net/mulan/home/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(SpatialPCA)
library(Seurat)
library(aricode)
library(ggplot2)
library(rhdf5)


prefix_dir = "/net/mulan/home/wenjinma"
#prefix_dir = "/media/sda/wenjinma"

data_dir = file.path(prefix_dir, 'data/spatial_STRS/')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_STRS_separate')
dir.create(result_dir, showWarnings = F)

Seurat_obj = readRDS(file.path(data_dir, 'STRS.RDS'))
Seurat_obj = JoinLayers(Seurat_obj)
ggplot(Seurat_obj@meta.data, aes(factor(timepoint), log10(nFeature_RNA+1))) + geom_violin()
ggsave(file.path(result_dir, 'nFeature_RNA_per_timepoint_log10.png'))
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)
# add metadata obtained from spaDOT
#metadata = read.csv(file.path(data_dir, 'Development_metadata.csv'), header=T, row.names=1)
# There seems to be no annotation
DimPlot(Seurat_obj, reduction = "umap", group.by='timepoint')
ggsave(file.path(result_dir, 'timepoint_umap.png'))
DimPlot(Seurat_obj, reduction = "umap", group.by='injury.zones')
ggsave(file.path(result_dir, 'injured_region_umap.png'))

# generate a list of count/locations
pc_df = data.frame(matrix(nrow=0, ncol=24))
colnames(pc_df) = c(paste0('spatialPC-', 1:20), 'timepoint', 'annotation', 'pixel_x', 'pixel_y')
for (tp in unique(Seurat_obj$timepoint)) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    tp_counts = tp_Seurat_obj@assays$RNA@layers$counts
    rownames(tp_counts) = rownames(tp_Seurat_obj)
    colnames(tp_counts) = colnames(tp_Seurat_obj)
    tp_locations = cbind(tp_Seurat_obj$pixel_x, tp_Seurat_obj$pixel_y)
    spatialPCA_obj = CreateSpatialPCAObject(counts=tp_counts, location=tp_locations,
                                            project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,
                                            gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
    spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
    cat('Timepoint:', tp, 'Bandwidth:', spatialPCA_obj@bandwidth, '\n')
    spatialPCA_obj = SpatialPCA_EstimateLoading(spatialPCA_obj,fast=TRUE,SpatialPCnum=20)
    spatialPCA_obj = SpatialPCA_SpatialPCs(spatialPCA_obj, fast=TRUE)
    colnames(spatialPCA_obj@SpatialPCs) = rownames(spatialPCA_obj@location)
    cluster_annotations = tp_Seurat_obj$injury.zones
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
                                                  latent_dat=spatialPCA_obj@SpatialPCs,
                                                  knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2])))
    cat(tp, 'ARI:', ARI(spatialPCA_clusterlabel, cluster_annotations), '\n')
    g = plot_cluster(location=tp_locations,
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=paste("SpatialPCA", tp),
                 color_in=rainbow(length(unique(cluster_annotations))),
                 legend='bottom')
    png(file.path(result_dir, paste0(tp, '_walktrap_SpatialPCA.png')))
    print(g)
    dev.off()
    tp_pc_df = as.data.frame(t(spatialPCA_obj@SpatialPCs))
    rownames(tp_pc_df) = gsub(paste0(tp, '_'), '', rownames(tp_pc_df))
    colnames(tp_pc_df) = paste0('spatialPC-', 1:20)
    tp_pc_df$timepoint = tp
    tp_pc_df$pixel_x = tp_locations[, 1]
    tp_pc_df$pixel_y = tp_locations[, 2]
    pc_df = rbind(pc_df, tp_pc_df)
}
write.csv(pc_df, file.path(result_dir, 'spatialPCA_pcs.csv'), quote=F)

