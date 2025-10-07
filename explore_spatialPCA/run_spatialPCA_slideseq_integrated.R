## conda activate /net/mulan/home/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(SpatialPCA)
library(Seurat)
library(anndata) # workaround for SeuratDisk as h5Seurat always has a issue of: Error: Missing required datasets 'levels' and 'values'
#library(Signac)
library(aricode)
library(ggplot2)
library(rhdf5)


prefix_dir = "/net/mulan/home/wenjinma"
#prefix_dir = "/media/sda/wenjinma"

data_dir = file.path(prefix_dir, 'data/spatial_Slideseq/')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_Slideseq_multiplesample_combinedlabels')
dir.create(result_dir, showWarnings = F)

#adata = read_h5ad(file.path(data_dir, 'slideseq.h5ad'))
#Seurat_obj = CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = adata$obs)
#saveRDS(Seurat_obj, file.path(data_dir, 'slideseq.RDS'))

Seurat_obj = readRDS(file.path(data_dir, 'slideseq.RDS'))
ggplot(Seurat_obj@meta.data, aes(factor(timepoint), log10(nFeature_RNA+1))) + geom_violin()
ggsave(file.path(result_dir, 'nFeature_RNA_per_timepoint_log10.png'))
Seurat_obj = subset(Seurat_obj, subset = nFeature_RNA > 100 & nFeature_RNA < 1000) # 128357 -> 123477 -> 9320 (> 200)
Seurat_obj$pixel_x = Seurat_obj$xcoord
Seurat_obj$pixel_y = Seurat_obj$ycoord
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

# generate a list of count/locations
count_list= list()
location_list = list()
seurat_list = list()
for (tp in unique(Seurat_obj$timepoint)) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    count_list[[tp]] = tp_Seurat_obj@assays$RNA@layers$counts
    rownames(count_list[[tp]]) = rownames(tp_Seurat_obj)
    colnames(count_list[[tp]]) = colnames(tp_Seurat_obj)
    location_list[[tp]] = cbind(tp_Seurat_obj$pixel_x, tp_Seurat_obj$pixel_y)
    seurat_list[[tp]] = tp_Seurat_obj
}

seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    print(paste0("Normalizing dataset"))
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
})
# select features that are repeatedly variable across datasets for integration
features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list,nfeatures = 10000)
# integration, find anchors
MultipleSample.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
# this command creates an 'integrated' data assay
MultipleSample.combined <- Seurat::IntegrateData(anchorset = MultipleSample.anchors)
DefaultAssay(MultipleSample.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
MultipleSample.combined <- Seurat::ScaleData(MultipleSample.combined, verbose = FALSE)
# obtain integrated and normalized data
integrated_data = MultipleSample.combined@assays$integrated@scale.data

# create each SpatialPCA object, this step is essentially for spatial gene selection
spatialpca_list = list()
for(count in names(seurat_list)){
    print(paste0("Creating SpatialPCA object for dataset ",count))
    spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
}

bandwidth_common=0.1
matched_data = list()
matched_location = list()
Kernal_mat = list()
for(count in names(count_list)){
    print(paste0("Creating SpatialPCA kernel for dataset ",count))
    match_gene = na.omit(match(rownames( spatialpca_list[[count]]@normalized_expr), rownames(integrated_data)))
    match_spot = grep(paste0("1_",count),colnames(integrated_data))
    matched_data[[count]] = integrated_data[match_gene,match_spot]
    matched_location[[count]] = spatialpca_list[[count]]@location
    spatialpca_list[[count]]@normalized_expr = t(scale(t(matched_data[[count]])))
    spatialpca_list[[count]]@location = scale(matched_location[[count]])
    spatialpca_list[[count]] = SpatialPCA_buildKernel(spatialpca_list[[count]], kerneltype="gaussian", bandwidth.set.by.user=bandwidth_common)
    Kernal_mat[[count]] = spatialpca_list[[count]]@kernelmat
}

# merge multiple samples
common_genes = unique(unlist(lapply(matched_data, function(x){ row.names(x)})))
MultipleSample_merge = spatialpca_list[[1]] # initialize using the first data, will replace with integrated data later
MultipleSample_merge@normalized_expr = t(scale(t(integrated_data[match(common_genes, rownames(integrated_data)),])))
MultipleSample_merge@kernelmat = Matrix::bdiag(Kernal_mat)
MultipleSample_merge@kernelmat = as.matrix(MultipleSample_merge@kernelmat)
# MultipleSample_merge@kernelmat = as(MultipleSample_merge@kernelmat, "sparseMatrix")
# MultipleSample_merge@sparseKernel = TRUE
MultipleSample_merge@params$expr = MultipleSample_merge@normalized_expr
MultipleSample_merge@location = do.call("rbind",lapply(spatialpca_list, function(x){ x@location}))
saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))

# perform SpatialPCA on integrated expression and integrated block diagonal kernel
MultipleSample_merge = SpatialPCA_EstimateLoading(MultipleSample_merge,fast=TRUE,SpatialPCnum=20)
MultipleSample_merge = SpatialPCA_SpatialPCs(MultipleSample_merge, fast=TRUE)
colnames(MultipleSample_merge@SpatialPCs) = rownames(MultipleSample_merge@location)

id_list = list()
SpatialPC_list = list()
Location_spatialpc_list = list()
for (sample in names(count_list)) {
    id_list[[sample]] = grep(sample, colnames(MultipleSample_merge@SpatialPCs))
    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]
    Location_spatialpc_list[[sample]] = MultipleSample_merge@location[id_list[[sample]],]
}
saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))
saveRDS(SpatialPC_list, file.path(result_dir, 'spatial_PCs.RDS'))
saveRDS(Location_spatialpc_list, file.path(result_dir, 'location_spatial_PCs.RDS'))

# write out integrated data in spatialPCA
write.csv(MultipleSample_merge@normalized_expr, file.path(result_dir, 'SpatialPCA_integrated_data.csv'), quote=F)
write.csv(MultipleSample_merge@kernelmat, file.path(result_dir, 'SpatialPCA_kernel.csv'), quote=F)

# split samples and match annotations
pc_df = data.frame(matrix(nrow=0, ncol=23))
colnames(pc_df) = c(paste0('spatialPC-', 1:20), 'timepoint', 'pixel_x', 'pixel_y')
for (sample in names(SpatialPC_list)) {
    cell_names = colnames(SpatialPC_list[[sample]])
    # editted information with metadata
    latent_dat = SpatialPC_list[[sample]][, cell_names]
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=5,
                                                  latent_dat=latent_dat,
                                                  knearest=round(sqrt(dim(latent_dat)[2])))
    # cluster_annotations = as.character(Seurat_obj@meta.data[cell_names, 'Annotation'])
    #spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
    #                                              latent_dat=SpatialPC_list[[sample]],
    #                                              knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    #spatialPCA_clusterlabel = louvain_clustering(clusternum=length(unique(cluster_annotations)),
    #                                              latent_dat=SpatialPC_list[[sample]],
    #                                              knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    sample_res = data.frame('cell_names'=cell_names, 'spatialPCA_labels'=spatialPCA_clusterlabel)
    g = plot_cluster(location=cbind(Seurat_obj[,cell_names]$pixel_x, Seurat_obj[,cell_names]$pixel_y),
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=paste("SpatialPCA", sample),
                 color_in=rainbow(5),
                 legend='bottom')
    library(grid)
    png(file.path(result_dir, paste0(sample, '_walktrap_SpatialPCA.png')))
    print(g)
    dev.off()
    tp_pc_df = as.data.frame(t(latent_dat))
    rownames(tp_pc_df) = gsub(paste0(sample, '_'), '', rownames(tp_pc_df))
    colnames(tp_pc_df) = paste0('spatialPC-', 1:20)
    tp_pc_df$timepoint = sample
    tp_pc_df$pixel_x = Seurat_obj[,cell_names]$pixel_x
    tp_pc_df$pixel_y = Seurat_obj[,cell_names]$pixel_y
    pc_df = rbind(pc_df, tp_pc_df)
}
write.csv(pc_df, file.path(result_dir, 'spatialPCA_pcs.csv'), quote=F)


