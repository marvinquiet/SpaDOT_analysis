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

data_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/spatialPCA_simulation_brain')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_simulation_brain_separate')
dir.create(result_dir, showWarnings = F)

# load data
t0_obj = readRDS(file.path(data_dir, 't0_brain_seed_2024.RDS'))
t0_obj@meta.data$timepoint = 't0'
t1_obj = readRDS(file.path(data_dir, 't1_brain_seed_2024.RDS'))
t1_obj@meta.data$timepoint = 't1'
t2_obj = readRDS(file.path(data_dir, 't2_brain_seed_2024.RDS'))
t2_obj@meta.data$timepoint = 't2'

Seurat_obj = merge(t0_obj, c(t1_obj, t2_obj))
Seurat_obj$pixel_x = Seurat_obj$x_coor
Seurat_obj$pixel_y = Seurat_obj$y_coor
Seurat_obj$annotation = Seurat_obj$label
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)
DimPlot(Seurat_obj, reduction = "umap", group.by='annotation')
ggsave(file.path(result_dir, 'annotation_umap.png'))
DimPlot(Seurat_obj, reduction = "umap", group.by='timepoint')
ggsave(file.path(result_dir, 'timepoint_umap.png'))

pc_df = data.frame(matrix(nrow=0, ncol=24))
colnames(pc_df) = c(paste0('spatialPC-', 1:20), 'annotation', 'timepoint', 'pixel_x', 'pixel_y')
for (tp in unique(Seurat_obj$timepoint)) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint == tp)
    tp_counts = tp_Seurat_obj@assays$RNA@layers$counts
    colnames(tp_counts) = colnames(tp_Seurat_obj)
    rownames(tp_counts) = rownames(tp_Seurat_obj)
    tp_locations = cbind(tp_Seurat_obj$pixel_x, tp_Seurat_obj$pixel_y)
    spatialPCA_obj = CreateSpatialPCAObject(counts=tp_counts, location=tp_locations,
                                            project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,
                                            gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
    #bandwidth_common=1 # same as the integrated data
    #spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidth.set.by.user=bandwidth_common)
    #spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
    spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidth.set.by.user=1, sparseKernel_tol=1e-5) 
    cat('Timepoint:', tp, 'Bandwidth:', spatialPCA_obj@bandwidth, '\n')
    spatialPCA_obj = SpatialPCA_EstimateLoading(spatialPCA_obj,fast=TRUE,SpatialPCnum=20)
    spatialPCA_obj = SpatialPCA_SpatialPCs(spatialPCA_obj, fast=TRUE)
    colnames(spatialPCA_obj@SpatialPCs) = rownames(spatialPCA_obj@location)
    cluster_annotations = tp_Seurat_obj@meta.data$annotation
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
                                                  latent_dat=spatialPCA_obj@SpatialPCs,
                                                  knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2])))
    g = plot_cluster(location=tp_locations,
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=tp, color_in=rainbow(length(unique(cluster_annotations))),
                 legend='bottom')
    library(grid)
    png(file.path(result_dir, paste0(tp, '_fixed_bandwidth1_sparse_SpatialPCA.png')))
    print(g)
    dev.off()
    cat('Timepoint:', tp, 'ARI:', ARI(spatialPCA_clusterlabel, cluster_annotations), '\n')
    tp_pc_df = data.frame(t(spatialPCA_obj@SpatialPCs))
    colnames(tp_pc_df) = paste0('spatialPC-', 1:20)
    tp_pc_df$annotation = cluster_annotations
    tp_pc_df$timepoint = tp
    tp_pc_df$pixel_x = tp_Seurat_obj$pixel_x
    tp_pc_df$pixel_y = tp_Seurat_obj$pixel_y
    pc_df = rbind(pc_df, tp_pc_df)
}
write.csv(pc_df, file.path(result_dir, 'spatialPCA_pcs_fixed_bandwith1_sparse.csv'), quote=F)




