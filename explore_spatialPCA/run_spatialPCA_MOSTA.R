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
result_dir = file.path(prefix_dir, "projects/spatialATAC-explore/results/SpatialPCA_MOSTA_fourtps/")
dir.create(result_dir, showWarnings = F)
prefix = 'HVG'

# load data
Seurat_obj = LoadH5Seurat(file.path(data_dir, 'Mouse_embryo_all_stage.h5seurat'))
annotation_counts = table(Seurat_obj$annotation)
selected_tissues = names(annotation_counts)[annotation_counts > 10000]
selected_idx = Seurat_obj$annotation %in% selected_tissues & Seurat_obj$timepoint %in% c('E12.5', 'E13.5', 'E14.5', 'E15.5')
filtered_Seurat_obj = Seurat_obj[, selected_idx]
set.seed(1234)
sampled_idx = sample(ncol(filtered_Seurat_obj), size=5000*4)
sampled_obj = filtered_Seurat_obj[, sampled_idx]
sampled_obj[["percent.mt"]] = PercentageFeatureSet(sampled_obj, pattern = "^MT-")

# perform filtering and clustering
sampled_obj = NormalizeData(sampled_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sampled_obj = FindVariableFeatures(sampled_obj, selection.method = "vst", nfeatures = 2000)
sampled_obj = ScaleData(sampled_obj, features = rownames(sampled_obj))
sampled_obj = RunPCA(sampled_obj, features = VariableFeatures(object = sampled_obj))
sampled_obj = RunUMAP(sampled_obj, dims = 1:20)

DimPlot(sampled_obj, reduction = "umap", group.by='annotation')
ggsave(file.path(result_dir, 'annotation_umap.png'))
DimPlot(sampled_obj, reduction = "umap", group.by='timepoint')
ggsave(file.path(result_dir, 'timepoint_umap.png'))

# construct input for SpatialPCA_multiple_sample.R
sampled_obj_genes = rownames(sampled_obj) # 23761
t1_sampled_obj = subset(sampled_obj, subset=timepoint == 'E12.5') # 6225 cells
t2_sampled_obj = subset(sampled_obj, subset=timepoint == 'E13.5') # 9958 cells
t3_sampled_obj = subset(sampled_obj, subset=timepoint == 'E14.5') # 9958 cells
t4_sampled_obj = subset(sampled_obj, subset=timepoint == 'E15.5') # 13817 cells

t1_counts = t1_sampled_obj@assays$count@counts[sampled_obj_genes, ]
t1_locations = Embeddings(t1_sampled_obj, reduction='spatial')
t2_counts = t2_sampled_obj@assays$count@counts[sampled_obj_genes, ]
t2_locations = Embeddings(t2_sampled_obj, reduction='spatial')
t3_counts = t3_sampled_obj@assays$count@counts[sampled_obj_genes, ]
t3_locations = Embeddings(t3_sampled_obj, reduction='spatial')
t4_counts = t4_sampled_obj@assays$count@counts[sampled_obj_genes, ]
t4_locations = Embeddings(t4_sampled_obj, reduction='spatial')

count_list = list(t1=t1_counts, t2=t2_counts, t3=t3_counts, t4=t4_counts)
location_list = list(t1=t1_locations, t2=t2_locations, t3=t3_locations, t4=t4_locations)

# release memory
rm(Seurat_obj, filtered_Seurat_obj, t1_sampled_obj, t2_sampled_obj, t3_sampled_obj, t4_sampled_obj)

#! error occurred when estimating loading
if(length(count_list) != length(location_list)) {
        stop("The number of count matrix should match with the number of location matrix. ")
}# end fi

# create seurat object for each dataset
seurat_list = list()
for(count in 1:length(count_list)){
    colnames(count_list[[count]]) = paste0("Sample",count,"_",colnames(count_list[[count]]))
    rownames(location_list[[count]]) = colnames(count_list[[count]])
    seurat_list[[count]] <- CreateSeuratObject(counts = count_list[[count]], project = "multiple")
}
# normalize and identify variable features for each dataset independently
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
integrated_data = MultipleSample.combined@assays$integrated@scale.data # 10000 x 20000

# create each SpatialPCA object, this step is essentially for spatial gene selection
spatialpca_list = list()
for(count in 1:length(seurat_list)){
    print(paste0("Creating SpatialPCA object for dataset ",count))
    #spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
    spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",gene.type="hvg",gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
}

bandwidth_common=0.1
matched_data = list()
matched_location = list()
Kernal_mat = list()
for(count in 1:length(count_list)){
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

saveRDS(MultipleSample_merge, file.path(result_dir, paste0(prefix, 'spatial_merge_obj.RDS')))

# perform SpatialPCA on integrated expression and integrated block diagonal kernel
MultipleSample_merge = SpatialPCA_EstimateLoading(MultipleSample_merge,fast=TRUE,SpatialPCnum=20)
MultipleSample_merge = SpatialPCA_SpatialPCs(MultipleSample_merge, fast=TRUE)
colnames(MultipleSample_merge@SpatialPCs) = rownames(MultipleSample_merge@location)

id_list = list()
SpatialPC_list = list()
Location_spatialpc_list = list()
for(sample in 1:length(count_list)){
    id_list[[sample]] = grep(paste0("Sample",sample), colnames(MultipleSample_merge@SpatialPCs))
    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]
    Location_spatialpc_list[[sample]] = spatialpca_list[[sample]]@location
}

names(SpatialPC_list) = paste0("Sample",1:length(SpatialPC_list))
names(Location_spatialpc_list) = paste0("Sample",1:length(Location_spatialpc_list))
saveRDS(MultipleSample_merge, file.path(result_dir, paste0(prefix, 'spatial_merge_obj.RDS')))
saveRDS(SpatialPC_list, file.path(result_dir, paste0(prefix, 'spatial_PCs.RDS')))
saveRDS(Location_spatialpc_list, file.path(result_dir, paste0(prefix, 'location_spatial_PCs.RDS')))

# split samples and match annotations
for (sample in names(SpatialPC_list)) {
    cell_names = gsub(paste0(sample, '_'), '', colnames(SpatialPC_list[[sample]]))
    cluster_annotations = as.character(sampled_obj@meta.data[cell_names, 'annotation'])
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
                                                  latent_dat=SpatialPC_list[[sample]],
                                                  knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    plot_cluster(location=Embeddings(sampled_obj, reduction='spatial')[cell_names, ],
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=paste("SpatialPCA", sample),
                 color_in=rainbow(length(unique(cluster_annotations))),
                 legend='bottom')
    ggsave(file.path(result_dir, paste0(prefix, sample, '_SpatialPCA.png')))
    cat(sample, 'ARI:', ARI(spatialPCA_clusterlabel, cluster_annotations), '\n')
}

MultipleSample_merge = readRDS(file.path(result_dir, paste0(prefix, 'spatial_merge_obj.RDS')))
## merged samples latent space
pcs = MultipleSample_merge@SpatialPCs
sample1_cells = gsub('Sample1_', '', colnames(pcs)[grep('Sample1', colnames(pcs))])
sample2_cells = gsub('Sample2_', '', colnames(pcs)[grep('Sample2', colnames(pcs))])
sample3_cells = gsub('Sample3_', '', colnames(pcs)[grep('Sample3', colnames(pcs))])
sample4_cells = gsub('Sample4_', '', colnames(pcs)[grep('Sample4', colnames(pcs))])
merge_cells = c(sample1_cells, sample2_cells, sample3_cells, sample4_cells)
merge_obj = sampled_obj[, merge_cells]
merge_obj = NormalizeData(merge_obj, normalization.method = "LogNormalize", scale.factor = 10000)
merge_obj = FindVariableFeatures(merge_obj, selection.method = "vst", nfeatures = 2000)
merge_obj = ScaleData(merge_obj, features = rownames(merge_obj))
# add spatialPCA pcs to merge_obj
colnames(pcs) = merge_cells
merge_obj[['spatialpca']] = CreateDimReducObject(embeddings = t(pcs), key = 'spatialPCA_', assay = 'RNA')
merge_obj = RunUMAP(merge_obj, reduction = "spatialpca", reduction.key='spatialumap_', dims = 1:20)
DimPlot(merge_obj, reduction = "umap", group.by='annotation')
ggsave(file.path(result_dir, paste0(prefix, 'spatialPCA_annotation_umap.png')))
DimPlot(merge_obj, reduction = "umap", group.by='timepoint')
ggsave(file.path(result_dir, paste0(prefix, 'spatialPCA_timepoint_umap.png')))

