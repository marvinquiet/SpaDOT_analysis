## conda activate /net/mulan/home/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(PRECAST)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)


base_dir = '/net/mulan/home/wenjinma'
project_dir = file.path(base_dir, "projects/spatialATAC-explore")
work_dir = file.path(project_dir, 'explore_PRECAST')
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

# construct input for SpatialPCA_multiple_sample.R
t1_Seurat_obj = subset(Seurat_obj, subset=section == '3') # 2397 cells
t2_Seurat_obj = subset(Seurat_obj, subset=section == '2') # 2312 cells
t3_Seurat_obj = subset(Seurat_obj, subset=section == '1') # 3353 cells
seurat_list = list(t1_Seurat_obj, t2_Seurat_obj, t3_Seurat_obj)

# run PRECAST
preobj = CreatePRECASTObject(seuList = seurat_list, selectGenesMethod = "HVGs", gene.number = 2000) 
PRECASTObj = AddAdjList(preobj, platform = "Visium")
PRECASTObj = AddParSetting(PRECASTObj, Sigma_equal = TRUE, coreNum = 1, maxIter = 30, verbose = TRUE)
PRECASTObj = PRECAST(PRECASTObj, K = 5) # set 5 for now because in total there will be 5 clusters
## backup the fitting results in resList
resList = PRECASTObj@resList
PRECASTObj = SelectModel(PRECASTObj)
seuInt = IntegrateSpaData(PRECASTObj, species='Unknown')
seuInt = AddTSNE(seuInt, n_comp = 2) 
seuInt@meta.data = Seurat_obj@meta.data[rownames(seuInt@meta.data), ]
seuInt@meta.data = seuInt@meta.data %>%
    mutate(timepoint = case_when(
        section == '3' ~ "t0",
        section == '2' ~ "t1",
        section == '1' ~ "t2"
))

# perform KMeans clustering on the embeddings
sample_res_all = NULL
ARI_res = c()
for (tp in unique(seuInt$timepoint)) {
    tp_seu = subset(seuInt, subset=timepoint == tp)
    tp_seu_embed = Embeddings(tp_seu, reduction='PRECAST')
    tp_kmeans = kmeans(tp_seu_embed, centers=length(unique(tp_seu$annotation)))
    tp_seu$cluster = tp_kmeans$cluster
    tp_ari = ARI(tp_seu$annotation, tp_seu$cluster)
    cat(tp, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # plot embeddings
    cols_cluster = chooseColors(palettes_name = 'Nature 10', n_colors = length(unique(tp_seu$annotation)), plot_colors = F)
    p1 = dimPlot(tp_seu, item='annotation', point_size = 1,  font_family='serif')
    p2 = dimPlot(tp_seu, item='cluster', font_family='serif', cols=cols_cluster) # Times New Roman
    drawFigs(list(p1, p2), common.legend=FALSE, align='hv', layout.dim=c(2,1)) 
    ggsave(file.path(result_dir, paste0(tp, '_tSNE.png')))
    # save sample res
    sample_res = cbind(tp_seu_embed, tp_seu@meta.data)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'PRECAST_results.csv'))
names(ARI_res) = unique(seuInt$timepoint)
write.csv(data.frame(ARI_res), file.path(result_dir, 'PRECAST_ARI.csv'))

