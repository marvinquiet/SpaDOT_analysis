## conda activate /media/sda/wenjinma/envs/spatialATAC
# .libPaths("/net/mulan/home/wenjinma/Rlib")
Sys.setenv(RETICULATE_PYTHON = "/media/sda/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)

library(SpatialPCA)
library(Seurat)
library(anndata) 
library(aricode)
library(ggplot2)
library(rhdf5)

source('/media/sda/wenjinma/projects/spatialATAC-explore/explore_spatialPCA_simulation/utils.R')

base_dir = '/media/sda/wenjinma'
project_dir = file.path(base_dir, "projects/spatialATAC-explore")
work_dir = file.path(project_dir, 'explore_spatialPCA')
setwd(work_dir)
set.seed(2024)
data_dir = file.path(base_dir, 'data/spatial_ChickenHeart')
result_dir = file.path(work_dir, 'ChickenHeart_results')
dir.create(result_dir, showWarnings = F)

adata = read_h5ad(file.path(data_dir, 'processed_ChickenHeart.h5ad'))
Seurat_obj = CreateSeuratObject(counts = t(as.matrix(adata$X)), meta.data = adata$obs)
Seurat_obj$annotation = Seurat_obj$region
id_to_tp = c('t0'='D4-A1', 't1'='D7-B1', 't2'='D10-C1', 't3'='D14-D1')

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
    id_list[[sample]] = grep(unname(id_to_tp[sample]), colnames(MultipleSample_merge@SpatialPCs))
    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]
    Location_spatialpc_list[[sample]] = MultipleSample_merge@location[id_list[[sample]],]
}
saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))
saveRDS(SpatialPC_list, file.path(result_dir, 'spatial_PCs.RDS'))
saveRDS(Location_spatialpc_list, file.path(result_dir, 'location_spatial_PCs.RDS'))

# --- analyze results
sample_res_all = NULL
ARI_res = c()
for (sample in names(SpatialPC_list)) {
    cell_names = colnames(SpatialPC_list[[sample]])
    tp_Seurat_obj = Seurat_obj[, cell_names]
    cluster_annotations = as.character(tp_Seurat_obj$annotation)
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
                                                  latent_dat=SpatialPC_list[[sample]],
                                                  knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    tp_Seurat_obj$walktrap = spatialPCA_clusterlabel
    #spatialPCA_clusterlabel = louvain_clustering(clusternum=length(unique(cluster_annotations)),
    #                                              latent_dat=SpatialPC_list[[sample]],
    #                                              knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    # run KMeans
    tp_kmeans = kmeans(t(SpatialPC_list[[sample]]), centers=length(unique(cluster_annotations)))
    tp_Seurat_obj$kmeans = unname(tp_kmeans$cluster)
    tp_ari = ARI(cluster_annotations, tp_kmeans$cluster)
    cat(sample, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # plot results
    g = plot_cluster(location=cbind(tp_Seurat_obj$pixel_x, tp_Seurat_obj$pixel_y),
                 clusterlabel=as.factor(tp_Seurat_obj$kmeans),
                 pointsize=1.5, text_size=20,
                 title_in=paste("SpatialPCA", sample),
                 color_in=rainbow(length(unique(cluster_annotations))),
                 legend='bottom')
    library(grid)
    png(file.path(result_dir, paste0(sample, '_KMeans_SpatialPCA.png')))
    print(g)
    dev.off()
    sample_res = cbind(t(SpatialPC_list[[sample]]), tp_Seurat_obj@meta.data)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'spatialPCA_results.csv'))
names(ARI_res) = names(SpatialPC_list)
write.csv(data.frame(ARI_res), file.path(result_dir, 'spatialPCA_ARI.csv'))