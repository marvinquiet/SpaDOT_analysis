## conda activate /media/sda/wenjinma/envs/spatialATAC
# 10/12/2024, run on comp1 cluster
# .libPaths("/net/mulan/home/wenjinma/Rlib")

library(SpatialPCA)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)
library(optparse)
library(Matrix)

source('/media/sda/wenjinma/projects/spatialATAC-explore/explore_spatialPCA_simulation/utils.R')

base_dir = '/media/sda/wenjinma'
project_dir = file.path(base_dir, "projects/spatialATAC-explore")
work_dir = file.path(project_dir, 'explore_spatialPCA')
dir.create(work_dir)
set.seed(2024)

# set optparse
option_list = list(
    make_option(c("-x", "--x_DEG"), type="integer", default=20, help="x% DEGs when having domain split"),
    make_option(c("-y", "--y_DEG"), type="integer", default=50, help="x% DEGs among domains"),
    make_option(c("-s", "--seed"), type="integer", default=2024, help="random seed")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)

seed = opt$seed
x_DEG = opt$x_DEG
y_DEG = opt$y_DEG
scenario_dir = file.path(base_dir, 'data', 'spatial_simulation_brain', seed, paste(x_DEG, y_DEG, sep='_'))
result_dir = file.path(work_dir, 'simulation_results', seed, paste(x_DEG, y_DEG, sep='_'))
dir.create(result_dir, recursive=T)

# load data
t0_mtx = readMM(file.path(scenario_dir, 't0.mtx'))
t0_metadata = read.csv(file.path(scenario_dir, 't0_metadata.csv'), header=T, row.names=1)
t0_metadata['timepoint'] = 't0'
t0_metadata = t0_metadata %>% rename(row=x_coor, col=y_coor)
t0_obj = CreateSeuratObject(t0_mtx, meta.data=t0_metadata)
colnames(t0_obj) = paste0('t0_', colnames(t0_obj))
t1_mtx = readMM(file.path(scenario_dir, 't1.mtx'))
t1_metadata = read.csv(file.path(scenario_dir, 't1_metadata.csv'), header=T, row.names=1)
t1_metadata['timepoint'] = 't1'
t1_metadata = t1_metadata %>% rename(row=x_coor, col=y_coor)
t1_obj = CreateSeuratObject(t1_mtx, meta.data=t1_metadata)
colnames(t1_obj) = paste0('t1_', colnames(t1_obj))
t2_mtx = readMM(file.path(scenario_dir, 't2.mtx'))
t2_metadata = read.csv(file.path(scenario_dir, 't2_metadata.csv'), header=T, row.names=1)
t2_metadata['timepoint'] = 't2'
t2_metadata = t2_metadata %>% rename(row=x_coor, col=y_coor)
t2_obj = CreateSeuratObject(t2_mtx, meta.data=t2_metadata)
colnames(t2_obj) = paste0('t2_', colnames(t2_obj))
# dowmsample to 9,000 spots
sample_idx = sample(c(colnames(t0_obj), colnames(t1_obj), colnames(t2_obj)), 3000*3)
Seurat_obj = merge(t0_obj, c(t1_obj, t2_obj))
Seurat_obj = Seurat_obj[, sample_idx]

# generate a list of count/locations
count_list= list()
location_list = list()
seurat_list = list()
for (tp in unique(Seurat_obj$timepoint)) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    count_list[[tp]] = tp_Seurat_obj@assays$RNA@layers$counts
    rownames(count_list[[tp]]) = rownames(tp_Seurat_obj)
    colnames(count_list[[tp]]) = colnames(tp_Seurat_obj)
    location_list[[tp]] = cbind(tp_Seurat_obj$row, tp_Seurat_obj$col)
    seurat_list[[tp]] = tp_Seurat_obj
}

start_time = proc.time()
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
    if (is.null(dim(spatialpca_list[[count]]@normalized_expr)) || nrow(spatialpca_list[[count]]@normalized_expr) < 50) {
        spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",gene.type="hvg",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
    }# not enough SVGs selected
}

bandwidth_common=0.1
matched_data = list()
matched_location = list()
Kernal_mat = list()
for(count in names(count_list)){
    print(paste0("Creating SpatialPCA kernel for dataset ",count))
    match_gene = na.omit(match(rownames( spatialpca_list[[count]]@normalized_expr), rownames(integrated_data)))
    match_spot = colnames(spatialpca_list[[count]]@normalized_expr)
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
MultipleSample_merge@params$expr = MultipleSample_merge@normalized_expr
MultipleSample_merge@location = do.call("rbind",lapply(spatialpca_list, function(x){ x@location}))
# saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))

# perform SpatialPCA on integrated expression and integrated block diagonal kernel
MultipleSample_merge = SpatialPCA_EstimateLoading(MultipleSample_merge,fast=TRUE,SpatialPCnum=20)
MultipleSample_merge = SpatialPCA_SpatialPCs(MultipleSample_merge, fast=TRUE)
colnames(MultipleSample_merge@SpatialPCs) = rownames(MultipleSample_merge@location)
end_time = proc.time()
elapsed_time = end_time - start_time
write(unname(elapsed_time["elapsed"]), file.path(result_dir, 'spatialPCA_running_time.txt'))

id_list = list()
SpatialPC_list = list()
Location_spatialpc_list = list()
for (sample in names(count_list)) {
    id_list[[sample]] = colnames(spatialpca_list[[sample]]@normalized_expr)
    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]
    Location_spatialpc_list[[sample]] = MultipleSample_merge@location[id_list[[sample]],]
}
saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))
saveRDS(SpatialPC_list, file.path(result_dir, 'spatial_PCs.RDS'))
saveRDS(Location_spatialpc_list, file.path(result_dir, 'location_spatial_PCs.RDS'))
write(rownames(MultipleSample_merge@normalized_expr), file.path(result_dir, 'spatialPCA_selected_genes.txt'))

# split samples and match labels
sample_res_all = NULL
ARI_res = c()
for (sample in names(SpatialPC_list)) {
    cell_names = colnames(SpatialPC_list[[sample]])
    tp_Seurat_obj = Seurat_obj[, cell_names]
    cluster_labels = as.character(tp_Seurat_obj$label)
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_labels)),
                                                  latent_dat=SpatialPC_list[[sample]],
                                                  knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    tp_Seurat_obj$walktrap = spatialPCA_clusterlabel
    # tp_ari = ARI(spatialPCA_clusterlabel, cluster_labels)
    # run KMeans
    tp_kmeans = kmeans(t(SpatialPC_list[[sample]]), centers=length(unique(cluster_labels)))
    tp_Seurat_obj$kmeans = tp_kmeans$cluster
    tp_ari = ARI(cluster_labels, tp_kmeans$cluster)
    cat(sample, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # plot results
    g = plot_cluster(location=cbind(Seurat_obj@meta.data[cell_names, 'row'], Seurat_obj@meta.data[cell_names, 'col']),
                 clusterlabel=as.factor(tp_kmeans$cluster),
                 pointsize=1.5, text_size=20,
                 title_in=paste("SpatialPCA", sample),
                 color_in=rainbow(length(unique(cluster_labels))),
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
