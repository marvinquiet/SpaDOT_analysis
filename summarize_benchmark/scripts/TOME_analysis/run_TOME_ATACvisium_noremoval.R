Sys.setenv(RETICULATE_PYTHON = "/media/sda/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(ClusterR)
library(dplyr)
library(future)
library(future.apply)
library(reshape2)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)

source('help_code.R')

base_dir = "/media/sda/wenjinma"
result_dir = file.path(base_dir, "projects/spatialATAC-explore")
downstream_dir = file.path(result_dir, "summarize_benchmark/ATAC-Visium_noremoval_figures/ATAC-Visium_noremoval_downstream", 'TOME')
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)
latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$pixel_x_bak = latent_adata$obs$pixel_x
latent_adata$obs$pixel_y_bak = latent_adata$obs$pixel_y

# --- draw ground truth
latent_adata$obs = latent_adata$obs %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "E12.5",
        timepoint == 't1' ~ "E13.5",
        timepoint == 't2' ~ "E15.5"
    ))
latent_adata$obs = latent_adata$obs %>%
    mutate(pixel_x = case_when(
        timepoint == 't0' ~ -pixel_y_bak,
        timepoint == 't2' ~ -pixel_y_bak,
        timepoint == 't1' ~ pixel_y_bak
    )) %>%
    mutate(pixel_y = -pixel_x_bak) 

latent_adata$obs$spaDOT_pred_labels = paste(latent_adata$obs$timepoint_original, latent_adata$obs$pred_labels, sep='_')
tp_color_palette = c("E12.5_1"="#F8766D", "E12.5_4"="#B79F00", "E12.5_3"="#00BA38", "E12.5_2"="#4EAEF0", "E12.5_0"="#D772EC", "E13.5_0"="#4EAEF0", "E13.5_1"="#F8766D", "E13.5_2"="#00BA38", "E13.5_3"="#D772EC", "E13.5_4"="#B79F00", "E15.5_1"="#F8766D", "E15.5_3"="#4EAEF0", "E15.5_4"="#B79F00", "E15.5_0"="#990099", "E15.5_2"="#FF99FF")

data_dir = file.path('/media/sda/wenjinma/data/spatial_scRNAseq_GSE214989')
h5_file = file.path(data_dir, "GSE214989_counts_embryo_visium.h5")
data = Read10X_h5(h5_file)
common_cells = intersect(colnames(data), rownames(latent_adata$obs))
data = data[, common_cells]
latent_adata = latent_adata[common_cells, ]
# --- create Seurat object for analysis
Seurat_obj = CreateSeuratObject(counts=data, meta.data=latent_adata$obs)
Seurat_obj[['spaDOT_latent']] = CreateAssayObject(data = t(latent_adata$X))
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)

# --- separate by timepoint
timepoints = as.vector(unique(Seurat_obj$timepoint))
timepoints = timepoints[order(timepoints)]
# --- generate data
for (tp_i in 1:(length(timepoints)-1)) {
    cur_tp = timepoints[tp_i]
    next_tp = timepoints[tp_i+1]
    cur_tp_obj = subset(Seurat_obj, subset=timepoint == cur_tp)
    next_tp_obj = subset(Seurat_obj, subset=timepoint == next_tp)

    # obj.list = list(cur_tp_obj, next_tp_obj)
    # obj.list = future_lapply(X = obj.list, FUN = function(x) {
    #     x <- NormalizeData(x, verbose = FALSE)
    #     x <- FindVariableFeatures(x, verbose = FALSE)
    # })
    # features <- SelectIntegrationFeatures(object.list = obj.list)
    # obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    #     x <- ScaleData(x, features = features, verbose = FALSE)
    #     x <- RunPCA(x, features = features, verbose = FALSE)
    # })
    # anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
    #                                   dims = 1:50)
    # obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
    
    # obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
    # obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
    # obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, n.components = 3, min.dist = 0.75)
    # saveRDS(obj.integrated, paste0(downstream_dir, "/", cur_tp, "_", next_tp, ".rds"))
    # emb = data.frame(Embeddings(object = obj.integrated, reduction = "umap"))
    # saveRDS(emb, file=paste0(downstream_dir, "/", cur_tp, "_", next_tp, "_umap3.rds"))

    emb = t(cbind(cur_tp_obj@assays$spaDOT_latent@data, next_tp_obj@assays$spaDOT_latent@data))

    # obtain annotation
    anno1 = cur_tp_obj@meta.data
    anno1$Anno = anno1$pred_labels
    anno1$day = "pre"
    anno1$stage = cur_tp
    anno2 = next_tp_obj@meta.data
    anno2$Anno = anno2$pred_labels
    anno2$day = "nex"
    anno2$stage = next_tp
    anno = rbind(anno1, anno2)
    anno = anno[rownames(emb),]

    res = createLineage_Knn(emb, anno,  k_neigh = 30) #### createLineage_Knn function was in help_code.R
    # saveRDS(res, paste0(downstream_dir, "/", cur_tp, "_", next_tp, "_Knn_umap.rds"))
    saveRDS(res, paste0(downstream_dir, "/", cur_tp, "_", next_tp, "_Knn_embedding.rds"))
}

replication_times=500
res_median_umap = list()
for (tp_i in 1:(length(timepoints)-1)) {
    cur_tp = timepoints[tp_i]
    next_tp = timepoints[tp_i+1]
    # dat = readRDS(paste0(downstream_dir, "/", cur_tp, "_", next_tp, "_Knn_umap.rds"))
    dat = readRDS(paste0(downstream_dir, "/", cur_tp, "_", next_tp, "_Knn_embedding.rds"))
    state_1 = row.names(dat[[1]])
    state_1 = paste0(next_tp, ":", state_1)
    state_2 = names(dat[[1]])
    state_2 = paste0(cur_tp, ":", state_2)
    tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))

    for(i in 1:nrow(dat[[1]])){
        for(j in 1:ncol(dat[[1]])){
        xx = NULL
        for(k in 1:replication_times){
            xx = c(xx, dat[[k]][i,j])
        }
        tmp_1[i,j] = median(xx[!is.na(xx)])
        }
    }
    tmp_1 = data.frame(tmp_1)
    row.names(tmp_1) = state_1
    names(tmp_1) = state_2
    res_median_umap[[cur_tp]] = tmp_1
}

dat = NULL
for(i in 1:length(res_median_umap)){
  print(timepoints[i])
  dat = rbind(dat, melt(as.matrix(res_median_umap[[i]])))
}

dat = data.frame(dat)
names(dat) = c("nex", "pre", "prob")

dat$pre_time = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][1]))
dat$pre_cell = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat$nex_time = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][1]))
dat$nex_cell = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][2]))
saveRDS(dat, paste0(downstream_dir, "/edge_all.rds"))

x = dat[dat$prob>=0.3,]
x = x[,c("pre","nex","prob")]
print(paste0("how many nodes now: ", length(unique(c(as.vector(x$pre), as.vector(x$nex))))))

dat_sub = res = x
dat_sub$pre_cell = unlist(lapply(as.vector(dat_sub$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat_sub$nex_cell = unlist(lapply(as.vector(dat_sub$nex), function(x) strsplit(x,"[:]")[[1]][2]))
sum(dat_sub$pre_cell == dat_sub$nex_cell)
sum(dat_sub$pre_cell == dat_sub$nex_cell)/nrow(dat_sub)
sum(dat_sub$pre_cell != dat_sub$nex_cell)
sum(dat_sub$pre_cell != dat_sub$nex_cell)/nrow(dat_sub)

print(paste0("how many edges: ", nrow(res)))
print(paste0("how many nodes: ", length(unique(c(as.vector(res$pre), as.vector(res$nex))))))

nex_list = as.vector(unique(res$nex))
tree = NULL
for(i in 1:length(nex_list)){
  res_sub = res[res$nex==nex_list[i],]
  if(nrow(res_sub)==1){
    tree = rbind(tree, res_sub)
  } else {
    res_sub = res_sub[order(res_sub$prob, decreasing = TRUE),]
    tree = rbind(tree, res_sub[1,])
  }
}
tree = data.frame(tree)
tree = tree[,c("pre","nex")]

write.table(res, paste0(downstream_dir, "/edge_prob.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(tree, paste0(downstream_dir, "/mouse_edge.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
