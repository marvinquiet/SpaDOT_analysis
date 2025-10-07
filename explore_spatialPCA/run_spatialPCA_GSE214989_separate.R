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

data_dir = file.path(prefix_dir, 'data/spatial_scRNAseq_GSE214989')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_GSE214989_separatesample')
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

DimPlot(Seurat_obj, reduction = "umap", group.by='seurat_clusters_harmony')
ggsave(file.path(result_dir, 'annotation_umap.png'))
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
pc_df = data.frame(matrix(nrow=0, ncol=24))
colnames(pc_df) = c(paste0('spatialPC-', 1:20), 'annotation', 'timepoint', 'pixel_x', 'pixel_y')
for (sec in unique(Seurat_obj$section)) {
    if (sec == 3) {tp = 't1'}
    if (sec == 2) {tp = 't2'}
    if (sec == 1) {tp = 't3'}
    tp_Seurat_obj = subset(Seurat_obj, subset=section == sec)
    tp_counts = tp_Seurat_obj@assays$RNA@counts[sampled_obj_genes, ]
    tp_locations = cbind(tp_Seurat_obj$pixel_x, tp_Seurat_obj$pixel_y)
    spatialPCA_obj = CreateSpatialPCAObject(counts=tp_counts, location=tp_locations,
                                            project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,
                                            gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
    bandwidth_common=0.1 # same as the integrated data
    #spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidth.set.by.user=bandwidth_common)
    spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
    cat('Timepoint:', tp, 'Bandwidth:', spatialPCA_obj@bandwidth, '\n')
    spatialPCA_obj = SpatialPCA_EstimateLoading(spatialPCA_obj,fast=FALSE,SpatialPCnum=20)
    spatialPCA_obj = SpatialPCA_SpatialPCs(spatialPCA_obj, fast=FALSE)
    colnames(spatialPCA_obj@SpatialPCs) = rownames(spatialPCA_obj@location)
    cluster_annotations = tp_Seurat_obj@meta.data$seurat_clusters_harmony
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_annotations)),
                                                  latent_dat=spatialPCA_obj@SpatialPCs,
                                                  knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2])))
    g = plot_cluster(location=tp_locations,
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=tp, color_in=rainbow(length(unique(cluster_annotations))),
                 legend='bottom')
    if (tp == 't1' || tp == 't3') {
        g = g + scale_y_reverse()
        library(grid)
        png(file.path(result_dir, paste0(tp, '_SpatialPCA.png')))
        print(g, vp=viewport(angle=-90))
        dev.off()
    } else {
        library(grid)
        png(file.path(result_dir, paste0(tp, '_SpatialPCA.png')))
        print(g, vp=viewport(angle=-90))
        dev.off()
    }
    #ggsave(file.path(result_dir, 'E12.5_SpatialPCA.png'))
    cat('Timepoint:', tp, 'ARI:', ARI(spatialPCA_clusterlabel, cluster_annotations), '\n')
    tp_pc_df = data.frame(t(spatialPCA_obj@SpatialPCs))
    colnames(tp_pc_df) = paste0('spatialPC-', 1:20)
    tp_pc_df$annotation = cluster_annotations
    tp_pc_df$timepoint = tp
    tp_pc_df$pixel_x = tp_Seurat_obj$pixel_x
    tp_pc_df$pixel_y = tp_Seurat_obj$pixel_y
    pc_df = rbind(pc_df, tp_pc_df)
}
write.csv(pc_df, file.path(result_dir, 'spatialPCA_pcs.csv'), quote=F)




