## conda activate /media/sda/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(SpatialPCA)
library(SeuratDisk)
library(Seurat)
#library(Signac)
library(aricode)
library(ggplot2)
library(rhdf5)
library(tidyverse)


prefix_dir = "/net/mulan/home/wenjinma"
#prefix_dir = "/media/sda/wenjinma"

data_dir = file.path(prefix_dir, 'data/spatial_scRNAseq_GSE214989')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_GSE214989_multiplesample_combinedlabels')
dir.create(result_dir, showWarnings = F)

# load data
h5_file = file.path(data_dir, "GSE214989_counts_embryo_visium.h5")
data = Read10X_h5(h5_file)
metadata = read.csv(file.path(data_dir, 'GSE214989_metadata_embryo_visium.csv.gz'),
                    header=T, row.names=1, sep=';')
common_cells = intersect(colnames(data), rownames(metadata))
data = data[, common_cells]
metadata = metadata[common_cells, ] # 8255 cells

# perform clustering
Seurat_obj = CreateSeuratObject(counts=data, meta.data=metadata)
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)


# from data analysis code
Seurat_obj$annotation = ifelse(Seurat_obj$seurat_clusters_harmony %in% c("3", "2","4"), 
                            "CNS",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% c("1", "9", "11"),
                            "Men/PNS",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% c("0","5","6", "8"),
                            "Mesenchyme",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% "7",
                            "Liver", "Heart"))))
DimPlot(Seurat_obj, reduction = "umap", group.by='seurat_clusters_harmony')
ggsave(file.path(result_dir, 'annotation_umap.png'))
DimPlot(Seurat_obj, reduction = "umap", group.by='annotation')
ggsave(file.path(result_dir, 'combined_annotation_umap.png'))
DimPlot(Seurat_obj, reduction = "umap", group.by='section')
ggsave(file.path(result_dir, 'section_umap.png'))


# add location information
t1_locations = read.csv(file.path(data_dir, 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz'),
                        header=F, sep=',')
t2_locations = read.csv(file.path(data_dir, 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz'),
                        header=F, sep=',')
t3_locations = read.csv(file.path(data_dir, 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz'),
                        header=F, sep=',')
colnames(t1_locations) = colnames(t2_locations) = colnames(t3_locations) = c("barcodes", "tissue", "y", "x", "pixel_y", "pixel_x")
rownames(t1_locations) = paste0(t1_locations$barcodes, '_1')
rownames(t2_locations) = paste0(t2_locations$barcodes, '_2')
rownames(t3_locations) = paste0(t3_locations$barcodes, '_3')
ggplot(t1_locations, aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't1_location_pixels.png'))

common_t1_cells = intersect(rownames(t1_locations), rownames(metadata))
ggplot(t1_locations[common_t1_cells, ], aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't1_location_pixels_common_cells.png'))  # E15.5 -> section 1 in metadata

common_t2_cells = intersect(rownames(t2_locations), rownames(metadata))
ggplot(t2_locations[common_t2_cells, ], aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't2_location_pixels_common_cells.png')) # E13.5 -> section 2 in metadata

common_t3_cells = intersect(rownames(t3_locations), rownames(metadata))
ggplot(t3_locations[common_t3_cells, ], aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't3_location_pixels_common_cells.png')) # E12.5 -> section 3 in metadata

# add locations to metadata
table(metadata$section) # 1: 3546, 2: 2312, 3: 2397

## so we change the barcodes
rownames(t1_locations) = paste0(t1_locations$barcodes, '_3')
rownames(t2_locations) = paste0(t2_locations$barcodes, '_2')
rownames(t3_locations) = paste0(t3_locations$barcodes, '_1')
ggplot(t1_locations, aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't1_location_pixels.png'))

common_t1_cells = intersect(rownames(t1_locations), rownames(metadata))
ggplot(t1_locations[common_t1_cells, ], aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't1_location_pixels_common_cells_curated.png'))  # 12.5
common_t2_cells = intersect(rownames(t2_locations), rownames(metadata))
ggplot(t2_locations[common_t2_cells, ], aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't2_location_pixels_common_cells_curated.png')) # E13.5
common_t3_cells = intersect(rownames(t3_locations), rownames(metadata))
ggplot(t3_locations[common_t3_cells, ], aes(x=pixel_x, y=pixel_y)) + geom_point()
ggsave(file.path(result_dir, 't3_location_pixels_common_cells_curated.png')) # E15.5

# add locations
metadata = metadata[c(common_t1_cells, common_t2_cells, common_t3_cells), ]
metadata$pixel_x = c(t1_locations[common_t1_cells, 'pixel_x'], t2_locations[common_t2_cells, 'pixel_x'], t3_locations[common_t3_cells, 'pixel_x'])
metadata$pixel_y = c(t1_locations[common_t1_cells, 'pixel_y'], t2_locations[common_t2_cells, 'pixel_y'], t3_locations[common_t3_cells, 'pixel_y'])
Seurat_obj@meta.data$pixel_x = metadata[colnames(Seurat_obj), 'pixel_x']
Seurat_obj@meta.data$pixel_y = metadata[colnames(Seurat_obj), 'pixel_y']

# construct input for SpatialPCA_multiple_sample.R
sampled_obj_genes = rownames(Seurat_obj)
t1_Seurat_obj = subset(Seurat_obj, subset=section == '3') # 2397 cells
t2_Seurat_obj = subset(Seurat_obj, subset=section == '2') # 2312 cells
t3_Seurat_obj = subset(Seurat_obj, subset=section == '1') # 3546 cells

t1_counts = t1_Seurat_obj@assays$RNA@layers$counts
rownames(t1_counts) = rownames(t1_Seurat_obj)
colnames(t1_counts) = colnames(t1_Seurat_obj)
t1_locations = cbind(t1_Seurat_obj$pixel_x, t1_Seurat_obj$pixel_y)
t2_counts = t2_Seurat_obj@assays$RNA@layers$counts
rownames(t2_counts) = rownames(t2_Seurat_obj)
colnames(t2_counts) = colnames(t2_Seurat_obj)
t2_locations = cbind(t2_Seurat_obj$pixel_x, t2_Seurat_obj$pixel_y)
t3_counts = t3_Seurat_obj@assays$RNA@layers$counts
rownames(t3_counts) = rownames(t3_Seurat_obj)
colnames(t3_counts) = colnames(t3_Seurat_obj)
t3_locations = cbind(t3_Seurat_obj$pixel_x, t3_Seurat_obj$pixel_y)
count_list = list(t1=t1_counts, t2=t2_counts, t3=t3_counts)
location_list = list(t1=t1_locations, t2=t2_locations, t3=t3_locations)

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
integrated_data = MultipleSample.combined@assays$integrated@scale.data

# create each SpatialPCA object, this step is essentially for spatial gene selection
spatialpca_list = list()
for(count in 1:length(seurat_list)){
    print(paste0("Creating SpatialPCA object for dataset ",count))
    spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
}

# save SVGs selected by sparkx
library(SPARK)
for (count in 1:length(seurat_list)) {
    Seu = seurat_list[[count]]
    Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    #count_test_spark = count_list[[count]][na.omit(match(rownames(Seu@assays$SCT@scale.data), rownames(count_list[[count]]))),
    #                                       na.omit(match(colnames(Seu@assays$SCT@scale.data), colnames(count_list[[count]])))]
    #location_test_spark = as.data.frame(location_list[[count]][match(colnames(Seu@assays$SCT@scale.data), rownames(location_list[[count]])),])
    #location_test_spark = as.matrix(location_test_spark)
    #sparkX = sparkx(count_test_spark, location_test_spark,numCores=5)
    #significant_gene_number = sum(sparkX$res_mtest$adjustedPval<=0.05)
    #SVGs = sparkX$res_mtest[order(sparkX$res_mtest$adjustedPval),][1:significant_gene_number, ]
    #write.csv(SVGs, file.path(result_dir, paste0(count, '_SVG_sparkx.csv')), quote=F)
    # --- perform Louvain clustering
    SVGs = read.csv(file.path(result_dir, paste0(count, '_SVG_sparkx.csv')), header=T, row.names=1)
    suppressMessages(require(FNN))
    suppressMessages(require(igraph))
    suppressMessages(require(bluster))
    SVG_mat = Seu@assays$SCT@scale.data[rownames(SVGs), ]
    SVG_pcs = prcomp(SVG_mat, center=T, scale=T)
    SVG_top_pcs = SVG_pcs$x[, 1:30]
    knn.norm = FNN::get.knn(as.matrix(SVG_top_pcs), k=30)
    knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
               k=30), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
    nw.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
    nw.norm = igraph::simplify(nw.norm)
    lc.norm = igraph::cluster_louvain(nw.norm)
    merged = bluster::mergeCommunities(nw.norm, lc.norm$membership, number=10)
    clusterlabel = as.character(as.integer(as.factor(paste0("cluster",merged))))
    SVGs$cluster = clusterlabel
    write.csv(SVGs, file.path(result_dir, paste0(count, '_SVG_sparkx_clustered_louvain.csv')), quote=F)
    # --- perform hierarchical clustering
    #gene_dist = dist(Seu@assays$SCT@scale.data[rownames(SVGs), ])
    ## dist_mat = as.matrix(gene_dist)
    #gene_hclust = hclust(gene_dist, method = "complete")
    #gene_cluster = cutree(gene_hclust, k = 10) %>% 
    #      # turn the named vector into a tibble
    #      enframe() %>% 
    #      # rename some of the columns
    #      rename(gene = name, cluster = value)
    #write.csv(as.data.frame(gene_cluster), file.path(result_dir, paste0(count, '_SVG_sparkx_clustered.csv')), quote=F)
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

saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))

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
saveRDS(MultipleSample_merge, file.path(result_dir, 'spatial_merge_obj.RDS'))
saveRDS(SpatialPC_list, file.path(result_dir, 'spatial_PCs.RDS'))
saveRDS(Location_spatialpc_list, file.path(result_dir, 'location_spatial_PCs.RDS'))

# write out integrated data in spatialPCA
write.csv(MultipleSample_merge@normalized_expr, file.path(result_dir, 'SpatialPCA_integrated_data.csv'), quote=F)
write.csv(MultipleSample_merge@kernelmat, file.path(result_dir, 'SpatialPCA_kernel.csv'), quote=F)

# split samples and match annotations
clusternum = c('Sample1'=4, 'Sample2'=5, 'Sample3'=5)
spatialPCA_res = data.frame(matrix(nrow=0, ncol=3))
colnames(spatialPCA_res) = c('cell_names', 'annotations', 'spatialPCA_labels')
for (sample in names(SpatialPC_list)) {
    cell_names = gsub(paste0(sample, '_'), '', colnames(SpatialPC_list[[sample]]))
    cluster_annotations = as.character(Seurat_obj@meta.data[cell_names, 'annotation'])
    #spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=unname(clusternum[sample]), # test with same number of clusters in spaDOT
                                                  latent_dat=SpatialPC_list[[sample]],
                                                  knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    #spatialPCA_clusterlabel = louvain_clustering(clusternum=length(unique(cluster_annotations)),
    #                                              latent_dat=SpatialPC_list[[sample]],
    #                                              knearest=round(sqrt(dim(SpatialPC_list[[sample]])[2])))
    sample_res = data.frame('cell_names'=cell_names, 'annotations'=cluster_annotations,
                            'spatialPCA_labels'=spatialPCA_clusterlabel)
    spatialPCA_res = rbind(spatialPCA_res, sample_res)
    g = plot_cluster(location=cbind(Seurat_obj@meta.data[cell_names, 'pixel_x'], Seurat_obj@meta.data[cell_names, 'pixel_y']),
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=paste("SpatialPCA", sample),
                 color_in=rainbow(length(unique(cluster_annotations))),
                 legend='bottom')
    if (sample == 'Sample1' || sample == 'Sample3') {
        g = g + scale_y_reverse()
        library(grid)
        png(file.path(result_dir, paste0(sample, '_louvain_SpatialPCA_editted.png')))
        print(g, vp=viewport(angle=-90))
        dev.off()
    } else {
        library(grid)
        png(file.path(result_dir, paste0(sample, '_louvain_SpatialPCA_editted.png')))
        print(g, vp=viewport(angle=-90))
        dev.off()
    }
    cat(sample, 'ARI:', ARI(spatialPCA_clusterlabel, cluster_annotations), '\n')
}

MultipleSample_merge = readRDS(file.path(result_dir, 'spatial_merge_obj.RDS'))
## merged samples latent space
pcs = MultipleSample_merge@SpatialPCs
sample1_cells = gsub('Sample1_', '', colnames(pcs)[grep('Sample1', colnames(pcs))])
sample2_cells = gsub('Sample2_', '', colnames(pcs)[grep('Sample2', colnames(pcs))])
sample3_cells = gsub('Sample3_', '', colnames(pcs)[grep('Sample3', colnames(pcs))])
merge_cells = c(sample1_cells, sample2_cells, sample3_cells)
merge_obj = Seurat_obj[, merge_cells]
merge_obj = NormalizeData(merge_obj, normalization.method = "LogNormalize", scale.factor = 10000)
merge_obj = FindVariableFeatures(merge_obj, selection.method = "vst", nfeatures = 2000)
merge_obj = ScaleData(merge_obj, features = rownames(merge_obj))
# add spatialPCA pcs to merge_obj
colnames(pcs) = merge_cells
merge_obj[['spatialpca']] = CreateDimReducObject(embeddings = t(pcs), key = 'spatialPCA_', assay = 'RNA')
merge_obj = RunUMAP(merge_obj, reduction = "spatialpca", reduction.key='spatialumap_', dims = 1:20)
DimPlot(merge_obj, reduction = "umap", group.by='annotation')
ggsave(file.path(result_dir, 'spatialPCA_annotation_umap.png'))
DimPlot(merge_obj, reduction = "umap", group.by='section')
ggsave(file.path(result_dir, 'spatialPCA_timepoint_umap.png'))

# save PCs
spatialPCs = readRDS(file.path(result_dir, 'spatial_PCs.RDS'))
pc_df = data.frame(matrix(nrow=0, ncol=24))
colnames(pc_df) = c(paste0('spatialPC-', 1:20), 'annotation', 'timepoint', 'pixel_x', 'pixel_y')
for (sample in names(spatialPCs)) {
    tp_pc_df = as.data.frame(t(spatialPCs[[sample]]))
    rownames(tp_pc_df) = gsub(paste0(sample, '_'), '', rownames(tp_pc_df))
    colnames(tp_pc_df) = paste0('spatialPC-', 1:20)
    tp_pc_df$annotation = Seurat_obj[, rownames(tp_pc_df)]$annotation
    if (sample == 'Sample1') timepoint = 't0'
    if (sample == 'Sample2') timepoint = 't1'
    if (sample == 'Sample3') timepoint = 't2'
    tp_pc_df$timepoint = timepoint
    tp_pc_df$pixel_x = Seurat_obj[, rownames(tp_pc_df)]$pixel_x
    tp_pc_df$pixel_y = Seurat_obj[, rownames(tp_pc_df)]$pixel_y
    pc_df = rbind(pc_df, tp_pc_df)
}
write.csv(pc_df, file.path(result_dir, 'spatialPCA_pcs.csv'), quote=F)

# explore gene clusters
