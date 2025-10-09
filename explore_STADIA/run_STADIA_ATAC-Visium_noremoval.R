## conda activate ~/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(stadia)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)
library(optparse)
library(Matrix)

## on Mulan
base_dir = '/net/mulan/home/wenjinma'
project_dir = file.path(base_dir, "projects/spatialATAC-explore")
work_dir = file.path(project_dir, 'explore_STADIA')
setwd(work_dir)
set.seed(2024)
data_dir = file.path(base_dir, 'data/spatial_scRNAseq_GSE214989')
result_dir = file.path(work_dir, 'ATAC-Visium_noremoval_results')
dir.create(result_dir, showWarnings = F)

h5_file = file.path(data_dir, "GSE214989_counts_embryo_visium.h5")
data = Read10X_h5(h5_file)
metadata = read.csv(file.path(data_dir, 'GSE214989_metadata_embryo_visium.csv.gz'),
                    header=T, row.names=1, sep=';')
common_cells = intersect(colnames(data), rownames(metadata))
data = data[, common_cells]
metadata = metadata[common_cells, ] # 8255 cells
# perform clustering
Seurat_obj = CreateSeuratObject(counts=data, meta.data=metadata)
Seurat_obj$annotation = ifelse(Seurat_obj$seurat_clusters_harmony %in% c("3", "2","4"), 
                            "CNS",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% c("1", "9", "11"),
                            "Men/PNS",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% c("0","5","6", "8"),
                            "Mesenchyme",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% "7",
                            "Liver", "Heart"))))
metadata = metadata[colnames(Seurat_obj), ]

# add location information
t1_locations = read.csv(file.path(data_dir, 'GSM6619680_220420_sATAC_V10B01-031_B1_tissue_positions_list.csv.gz'),
                        header=F, sep=',')
t2_locations = read.csv(file.path(data_dir, 'GSM6619681_211007_V10S29-086_D1_tissue_positions_list.csv.gz'),
                        header=F, sep=',')
t3_locations = read.csv(file.path(data_dir, 'GSM6619682_V10B01-135_D1_tissue_positions_list.csv.gz'),
                        header=F, sep=',')
colnames(t1_locations) = colnames(t2_locations) = colnames(t3_locations) = c("barcodes", "tissue", "y", "x", "pixel_y", "pixel_x")
rownames(t1_locations) = paste0(t1_locations$barcodes, '_3')
rownames(t2_locations) = paste0(t2_locations$barcodes, '_2')
rownames(t3_locations) = paste0(t3_locations$barcodes, '_1')
common_t1_cells = intersect(rownames(t1_locations), rownames(metadata))
common_t2_cells = intersect(rownames(t2_locations), rownames(metadata))
common_t3_cells = intersect(rownames(t3_locations), rownames(metadata))
# add locations
metadata = metadata[c(common_t1_cells, common_t2_cells, common_t3_cells), ]
metadata$pixel_x = c(t1_locations[common_t1_cells, 'pixel_x'], t2_locations[common_t2_cells, 'pixel_x'], t3_locations[common_t3_cells, 'pixel_x'])
metadata$pixel_y = c(t1_locations[common_t1_cells, 'pixel_y'], t2_locations[common_t2_cells, 'pixel_y'], t3_locations[common_t3_cells, 'pixel_y'])
Seurat_obj@meta.data$row = metadata[colnames(Seurat_obj), 'pixel_x']
Seurat_obj@meta.data$col = metadata[colnames(Seurat_obj), 'pixel_y']

# construct input for STADIA
t1_Seurat_obj = subset(Seurat_obj, subset=section == '3') # 2397 cells
t1_Seurat_obj$timepoint = 't0'
t2_Seurat_obj = subset(Seurat_obj, subset=section == '2') # 2312 cells
t2_Seurat_obj$timepoint = 't1'
t3_Seurat_obj = subset(Seurat_obj, subset=section == '1') # 3353 cells
t3_Seurat_obj$timepoint = 't2'
seuVec = list(t1_Seurat_obj, t2_Seurat_obj, t3_Seurat_obj)
seuList = list()
seuMetaList = list()
for (i in 1:length(seuVec)) {
    sub_obj = seuVec[[i]]
    sub_assay_v3 = CreateAssayObject(counts=sub_obj[['RNA']]$'counts')
    sub_obj[['RNA']] = sub_assay_v3
    seuList[[i]] = sub_obj
    seuMetaList[[i]] = sub_obj@meta.data
}
seuMeta = do.call(rbind, seuMetaList)

## set parameters for STADIA model
etas = 0.15
K = 5
hyper = HyperParameters(seuList, eta = etas)
stadia = stadia(seuList, hyper, n_cluster = K, min.features = 0, min.spots = 20, # do not remove cells based on feature expression
                platform = "others", em.maxiter = 30)
save(stadia, file=file.path(result_dir, 'stadia.RData'))

# perform KMeans clustering on the embeddings
sample_res_all = NULL
ARI_res = c()
for (tp in unique(seuMeta$timepoint)) {
    tp_meta = seuMeta[seuMeta$timepoint == tp, ]
    common_cells = intersect(rownames(tp_meta), colnames(stadia$factors))
    tp_embed = t(stadia$factors)[common_cells, ]
    tp_meta = tp_meta[common_cells, ]
    tp_kmeans = kmeans(tp_embed, centers=length(unique(tp_meta$annotation)))
    tp_meta$cluster = tp_kmeans$cluster
    tp_ari = ARI(tp_meta$annotation, tp_meta$cluster)
    cat(tp, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # visualizing results
    g = ggplot(tp_meta, aes(row, col)) + 
        geom_point(size = .3, shape = 21, color = "#00000000", aes(fill = factor(cluster))) +
        labs(subtitle = paste0("ARI:", round(tp_ari, 2)))
    ggsave(file.path(result_dir, paste0(tp, '_results.png')))
    # save embeddings
    sample_res = cbind(tp_embed, tp_meta)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'STADIA_results.csv'))
names(ARI_res) = unique(seuMeta$timepoint)
write.csv(data.frame(ARI_res), file.path(result_dir, 'STADIA_ARI.csv'))

