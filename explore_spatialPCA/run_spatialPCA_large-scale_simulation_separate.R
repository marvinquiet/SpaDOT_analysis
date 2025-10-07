## conda activate /media/sda/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(SpatialPCA)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)
library(optparse)
library(Matrix)

source('/net/mulan/home/wenjinma/projects/spatialATAC-explore/explore_spatialPCA_simulation/utils.R')

## on Mulan
base_dir = '/net/mulan/home/wenjinma'
project_dir = file.path(base_dir, "projects/spatialATAC-explore")
work_dir = file.path(project_dir, 'explore_spatialPCA_separate')
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
scenario_dir = file.path(project_dir, 'results', 'spatialPCA_simulation_brain', seed, paste(x_DEG, y_DEG, sep='_'))
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
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)
DimPlot(Seurat_obj, reduction = "umap", group.by='label')
ggsave(file.path(result_dir, 'label_umap.png'))
DimPlot(Seurat_obj, reduction = "umap", group.by='timepoint')
ggsave(file.path(result_dir, 'timepoint_umap.png'))

# generate a list of count/locations
sample_res_all = NULL
ARI_res = c()
for (tp in unique(Seurat_obj$timepoint)) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    tp_counts = tp_Seurat_obj@assays$RNA@layers$counts
    rownames(tp_counts) = rownames(tp_Seurat_obj)
    colnames(tp_counts) = colnames(tp_Seurat_obj)
    tp_locations = cbind(tp_Seurat_obj$row, tp_Seurat_obj$col)
    spatialPCA_obj = CreateSpatialPCAObject(counts=tp_counts, location=tp_locations,
                                            project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,
                                            gene.number=3000, customGenelist=NULL,min.loctions=0, min.features=0)
    spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
    cat('Timepoint:', tp, 'Bandwidth:', spatialPCA_obj@bandwidth, '\n')
    spatialPCA_obj = SpatialPCA_EstimateLoading(spatialPCA_obj,fast=TRUE,SpatialPCnum=20)
    spatialPCA_obj = SpatialPCA_SpatialPCs(spatialPCA_obj, fast=TRUE)
    colnames(spatialPCA_obj@SpatialPCs) = rownames(spatialPCA_obj@location)
    cluster_labels = tp_Seurat_obj$label
    spatialPCA_clusterlabel = walktrap_clustering(clusternum=length(unique(cluster_labels)),
                                                  latent_dat=spatialPCA_obj@SpatialPCs,
                                                  knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2])))
    g = plot_cluster(location=tp_locations,
                 clusterlabel=spatialPCA_clusterlabel,
                 pointsize=1.5, text_size=20,
                 title_in=tp, color_in=rainbow(length(unique(cluster_labels))),
                 legend='bottom')
    library(grid)
    png(file.path(result_dir, paste0(tp, '_walktrap_SpatialPCA.png')))
    print(g)
    dev.off()
    tp_ari = ARI(spatialPCA_clusterlabel, cluster_labels)
    cat('Timepoint:', tp, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    tp_Seurat_obj$cluster = spatialPCA_clusterlabel
    sample_res = cbind(t(spatialPCA_obj@SpatialPCs), tp_Seurat_obj@meta.data)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'spatialPCA_results.csv'))
names(ARI_res) = unique(Seurat_obj$timepoint)
write.csv(data.frame(ARI_res), file.path(result_dir, 'spatialPCA_ARI.csv'))

