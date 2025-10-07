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
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_simulation_brain_integrated')
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
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list,nfeatures = 2000)
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
sample_dict = c("t0"="1", "t1"="2", "t2"="3")
for (sample in names(count_list)) {
    idx_vect = sapply(strsplit(colnames(MultipleSample_merge@SpatialPCs), split='_'), '[', 2)
    id_list[[sample]] = grep(unname(sample_dict[sample]), idx_vect)
    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]
    Location_spatialpc_list[[sample]] = MultipleSample_merge@location[id_list[[sample]],]
}
saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))
saveRDS(SpatialPC_list, file.path(result_dir, 'spatial_PCs.RDS'))
saveRDS(Location_spatialpc_list, file.path(result_dir, 'location_spatial_PCs.RDS'))

# write out integrated data in spatialPCA
#write.csv(MultipleSample_merge@normalized_expr, file.path(result_dir, 'SpatialPCA_integrated_data.csv'), quote=F)
#write.csv(MultipleSample_merge@kernelmat, file.path(result_dir, 'SpatialPCA_kernel.csv'), quote=F)

write(rownames(MultipleSample_merge@normalized_expr), file.path(result_dir, 'spatialPCA_selected_genes.txt'))

# split samples and match annotations
for (sample in names(SpatialPC_list)) {
    cell_names = colnames(SpatialPC_list[[sample]])
    cluster_annotations = as.character(Seurat_obj@meta.data[cell_names, 'annotation'])
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
                                                  latent_dat=SpatialPC_list[[sample]],
                                                  knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    #spatialPCA_clusterlabel = louvain_clustering(clusternum=length(unique(cluster_annotations)),
    #                                              latent_dat=SpatialPC_list[[sample]],
    #                                              knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    sample_res = data.frame('cell_names'=cell_names, 'annotations'=cluster_annotations,
                            'spatialPCA_labels'=spatialPCA_clusterlabel)
    g = plot_cluster(location=cbind(Seurat_obj@meta.data[cell_names, 'pixel_x'], Seurat_obj@meta.data[cell_names, 'pixel_y']),
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=paste("SpatialPCA", sample),
                 color_in=rainbow(length(unique(cluster_annotations))),
                 legend='bottom')
    library(grid)
    png(file.path(result_dir, paste0(sample, '_walktrap_SpatialPCA.png')))
    print(g)
    dev.off()
    cat(sample, 'ARI:', ARI(spatialPCA_clusterlabel, cluster_annotations), '\n')
}

pcs = MultipleSample_merge@SpatialPCs
merge_obj = Seurat_obj[, colnames(pcs)]
merge_obj = NormalizeData(merge_obj, normalization.method = "LogNormalize", scale.factor = 10000)
merge_obj = FindVariableFeatures(merge_obj, selection.method = "vst", nfeatures = 2000)
merge_obj = ScaleData(merge_obj, features = rownames(merge_obj))
# add spatialPCA pcs to merge_obj
colnames(pcs) = merge_cells
merge_obj[['spatialpca']] = CreateDimReducObject(embeddings = t(pcs), key = 'spatialPCA_', assay = 'RNA')
merge_obj = RunUMAP(merge_obj, reduction = "spatialpca", reduction.key='spatialumap_', dims = 1:20)
DimPlot(merge_obj, reduction = "umap", group.by='annotation')
ggsave(file.path(result_dir, 'spatialPCA_annotation_umap.png'))
DimPlot(merge_obj, reduction = "umap", group.by='timepoint')
ggsave(file.path(result_dir, 'spatialPCA_timepoint_umap.png'))
pc_df = as.data.frame(t(pcs))
pc_df$annotation = merge_obj$annotation
pc_df$timepoint = merge_obj$timepoint
pc_df$pixel_x = merge_obj$pixel_x
pc_df$pixel_y = merge_obj$pixel_y
write.csv(pc_df, file.path(result_dir, 'spatialPCA_pcs.csv'), quote=F)

