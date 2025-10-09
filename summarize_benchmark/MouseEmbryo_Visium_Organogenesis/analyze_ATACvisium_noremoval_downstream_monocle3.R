.libPaths("/net/mulan/home/wenjinma/Rlib")

Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(monocle3)
library(dplyr)
library(reshape2)


base_dir = "/net/mulan/home/wenjinma"
result_dir = file.path(base_dir, "projects/spatialATAC-explore")
spatial_layour_dir = file.path(result_dir, "summarize_benchmark/ATAC-Visium_noremoval_figures/spatial_layout")
dir.create(spatial_layour_dir, showWarnings = FALSE, recursive = TRUE)
downstream_dir = file.path(result_dir, "summarize_benchmark/ATAC-Visium_noremoval_figures/ATAC-Visium_noremoval_downstream")
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)
latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$pixel_x_bak = latent_adata$obs$pixel_x
latent_adata$obs$pixel_y_bak = latent_adata$obs$pixel_y
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


# --- explore trajectory analysis on t0 cluster 0, t1 cluster3, t2 cluster2 and t2 cluster 0
t0_cluster0_idx = which(latent_adata$obs$timepoint == 't0' & latent_adata$obs$pred_labels == 0)
t1_cluster3_idx = which(latent_adata$obs$timepoint == 't1' & latent_adata$obs$pred_labels == 3)
t2_cluster0_idx = which(latent_adata$obs$timepoint == 't2' & latent_adata$obs$pred_labels == 0)
t2_cluster2_idx = which(latent_adata$obs$timepoint == 't2' & latent_adata$obs$pred_labels == 2)
total_idx = c(t0_cluster0_idx, t1_cluster3_idx, t2_cluster0_idx, t2_cluster2_idx)
adata = latent_adata[total_idx, ]

# --- obtain marker genes
t0_markers = read.csv(file.path(downstream_dir, 't0_Seurat_FindAllMarkers_original_cluster.csv'),
                          header=TRUE, row.names=1)
sub_t0_markers = t0_markers[t0_markers$avg_log2FC > 1 & t0_markers$p_val_adj < 0.05,]
t0_marker_list = list()
for (cluster in unique(sub_t0_markers$cluster)) {
    selected_markers = sub_t0_markers[sub_t0_markers$cluster == cluster,]$gene
    t0_marker_list[[paste0('t0_cluster', cluster)]] = selected_markers[1:100]
}
t1_markers = read.csv(file.path(downstream_dir, 't1_Seurat_FindAllMarkers_original_cluster.csv'),
                            header=TRUE, row.names=1)
sub_t1_markers = t1_markers[t1_markers$avg_log2FC > 1 & t1_markers$p_val_adj < 0.05,]
t1_marker_list = list()
for (cluster in unique(sub_t1_markers$cluster)) {
    selected_markers = sub_t1_markers[sub_t1_markers$cluster == cluster,]$gene
    t1_marker_list[[paste0('t1_cluster', cluster)]] = selected_markers[1:100]
}
t2_markers = read.csv(file.path(downstream_dir, 't2_Seurat_FindAllMarkers_original_cluster.csv'),
                            header=TRUE, row.names=1)
sub_t2_markers = t2_markers[t2_markers$avg_log2FC > 1 & t2_markers$p_val_adj < 0.05,]
t2_marker_list = list()
for (cluster in unique(sub_t2_markers$cluster)) {
    selected_markers = sub_t2_markers[sub_t2_markers$cluster == cluster,]$gene
    t2_marker_list[[paste0('t2_cluster', cluster)]] = selected_markers[1:100]
}

t0_t1_overlap_markers <- sapply(t1_marker_list, function(cluster2) {
  sapply(t0_marker_list, function(cluster1) {
    length(intersect(cluster1, cluster2)) # Find overlapping genes
  })
})
t0_t1_overlap_markers_df = melt(t0_t1_overlap_markers)
t0_t1_overlap_markers_df$Var1 = substr(t0_t1_overlap_markers_df$Var1, nchar(as.character(t0_t1_overlap_markers_df$Var1)), nchar(as.character(t0_t1_overlap_markers_df$Var1)))
t0_t1_overlap_markers_df$Var2 = substr(t0_t1_overlap_markers_df$Var2, nchar(as.character(t0_t1_overlap_markers_df$Var2)), nchar(as.character(t0_t1_overlap_markers_df$Var2)))
ggplot(t0_t1_overlap_markers_df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = value), color = "white") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 100)) +
    labs(x = "E13.5 Clusters", y = "E12.5 Clusters", title = "Overlapping Marker Genes (top 100)") +
    theme_pubr() + 
    theme(legend.position="right")
ggsave(file.path(downstream_dir, 't0_t1_overlap_markers_heatmap.png'), width=4.5, height=4)

t1_t2_overlap_markers = sapply(t2_marker_list, function(cluster2) {
  sapply(t1_marker_list, function(cluster1) {
    length(intersect(cluster1, cluster2)) # Find overlapping genes
  })
})
t1_t2_overlap_markers_df = melt(t1_t2_overlap_markers)
t1_t2_overlap_markers_df$Var1 = substr(t1_t2_overlap_markers_df$Var1, nchar(as.character(t1_t2_overlap_markers_df$Var1)), nchar(as.character(t1_t2_overlap_markers_df$Var1)))
t1_t2_overlap_markers_df$Var2 = substr(t1_t2_overlap_markers_df$Var2, nchar(as.character(t1_t2_overlap_markers_df$Var2)), nchar(as.character(t1_t2_overlap_markers_df$Var2)))
ggplot(t1_t2_overlap_markers_df, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = value), color = "white") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 100)) +
    labs(x = "E15.5 Clusters", y = "E13.5 Clusters", title = "Overlapping Marker Genes (top 100)") +
    theme_pubr() + 
    theme(legend.position="right")
ggsave(file.path(downstream_dir, 't1_t2_overlap_markers_heatmap.png'), width=4.5, height=4)

# --- run enrichGO for each marker gene list
for (cluster in names(t0_marker_list)) {
    t0_cluster_markers = t0_marker_list[[cluster]]
    t0_cluster_GO = enrichGO(gene = t0_cluster_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
    g = dotplot(t0_cluster_GO, showCategory=10)
    ggsave(file.path(downstream_dir, paste0('t0_', cluster, '_top100_enrichment_clusterprofiler.png')), width=8, height=6)
}

for (cluster in names(t1_marker_list)) {
    t1_cluster_markers = t1_marker_list[[cluster]]
    t1_cluster_GO = enrichGO(gene = t1_cluster_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
    g = dotplot(t1_cluster_GO, showCategory=10)
    ggsave(file.path(downstream_dir, paste0('t1_', cluster, '_top100_enrichment_clusterprofiler.png')), width=8, height=6)
}


for (cluster in names(t2_marker_list)) {
    t2_cluster_markers = t2_marker_list[[cluster]]
    t2_cluster_GO = enrichGO(gene = t2_cluster_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
    g = dotplot(t2_cluster_GO, showCategory=10)
    ggsave(file.path(downstream_dir, paste0('t2_', cluster, '_top100_enrichment_clusterprofiler.png')), width=8, height=6)
}

# --- venn diagram for top 100 marker genes
t0_cluster0_markers = sub_t0_markers[sub_t0_markers$cluster == 0,]$gene
t1_cluster3_markers = sub_t1_markers[sub_t1_markers$cluster == 3,]$gene
t2_cluster0_markers = sub_t2_markers[sub_t2_markers$cluster == 0,]$gene
t2_cluster2_markers = sub_t2_markers[sub_t2_markers$cluster == 2,]$gene

markers_list = list(
    t0_cluster0 = t0_cluster0_markers,
    t1_cluster3 = t1_cluster3_markers,
    t2_cluster0 = t2_cluster0_markers,
    t2_cluster2 = t2_cluster2_markers
)
venn.plot = venn.diagram(
  x = markers_list,
  category.names = names(markers_list),
  filename = file.path(downstream_dir, 't0_cluster0_VennDiagram.png'),  # Keep output in R instead of saving to a file
  output = TRUE
)


# --- marker genes overlap between timepoints in each cluster
library(clusterProfiler)
library(org.Mm.eg.db)


t0_t1_common_markers = intersect(t0_cluster0_markers,t1_cluster3_markers)
t0_cluster0_specific_markers = setdiff(t0_cluster0_markers, t0_t1_common_markers)
t1_cluster3_specific_markers = setdiff(t1_cluster3_markers, t0_t1_common_markers)

t0_t1_common_markers_GO = enrichGO(gene = t0_t1_common_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
g = dotplot(t0_t1_common_markers_GO, showCategory=10)
ggsave(file.path(downstream_dir, 't0_t1_common_markers_enrichment_clusterprofiler.png'), width=8, height=6)

t0_specific_GO = enrichGO(gene = t0_cluster0_specific_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
g = dotplot(t0_specific_GO, showCategory=10)
ggsave(file.path(downstream_dir, 't0_specific_markers_enrichment_clusterprofiler.png'), width=8, height=6)

t1_specific_GO = enrichGO(gene = t1_cluster3_specific_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
g = dotplot(t1_specific_GO, showCategory=10)
ggsave(file.path(downstream_dir, 't1_specific_markers_enrichment_clusterprofiler.png'), width=8, height=6)


t2_cluster_genes = c(t2_cluster0_markers, t2_cluster2_markers)
t1_t2_common_markers = intersect(t1_cluster3_markers,t2_cluster_genes)
t1_specific_markers = setdiff(t1_cluster3_markers, t1_t2_common_markers)
t2_cluster0_specific_markers = setdiff(t2_cluster0_markers, t1_t2_common_markers)
t2_cluster2_specific_markers = setdiff(t2_cluster2_markers, t1_t2_common_markers)

t1_t2_common_markers_GO = enrichGO(gene = t1_t2_common_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
g = dotplot(t1_t2_common_markers_GO, showCategory=10)
ggsave(file.path(downstream_dir, 't1_t2_common_markers_enrichment_clusterprofiler.png'), width=8, height=6)

t2_cluster0_specific_GO = enrichGO(gene = t2_cluster0_specific_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
g = dotplot(t2_cluster0_specific_GO, showCategory=10)
ggsave(file.path(downstream_dir, 't2_cluster0_specific_markers_enrichment_clusterprofiler.png'), width=8, height=6)

t2_cluster2_specific_GO = enrichGO(gene = t2_cluster2_specific_markers, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
g = dotplot(t2_cluster2_specific_GO, showCategory=10)
ggsave(file.path(downstream_dir, 't2_cluster2_specific_markers_enrichment_clusterprofiler.png'), width=8, height=6)


# --- run monocle3 on the selected clusters
t0_cluster0_idx = which(latent_adata$obs$timepoint == 't0' & latent_adata$obs$pred_labels == 0)
t1_cluster2_idx = which(latent_adata$obs$timepoint == 't1' & latent_adata$obs$pred_labels == 2)
t1_cluster3_idx = which(latent_adata$obs$timepoint == 't1' & latent_adata$obs$pred_labels == 3)
t1_cluster4_idx = which(latent_adata$obs$timepoint == 't1' & latent_adata$obs$pred_labels == 4)
t2_cluster0_idx = which(latent_adata$obs$timepoint == 't2' & latent_adata$obs$pred_labels == 0)
t2_cluster2_idx = which(latent_adata$obs$timepoint == 't2' & latent_adata$obs$pred_labels == 2)

# first focus on t0_cluster0 -> t1_cluster2/3/4
total_idx = c(t0_cluster0_idx, t1_cluster2_idx, t1_cluster3_idx, t1_cluster4_idx)
adata = latent_adata[total_idx, ]
cluster_markers = c()
for (tp in c('t0', 't1')) {
    tp_markers = read.csv(file.path(downstream_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')),
                          header=TRUE, row.names=1)
    subset_tp_markers = tp_markers[tp_markers$avg_log2FC > 1 & tp_markers$p_val_adj < 0.05,]
    if (tp == 't0') {
        selected_markers = subset_tp_markers[subset_tp_markers$cluster == 0,]
        cluster_markers = c(cluster_markers, selected_markers$gene[1:100])
    }
    if (tp == 't1') {
        selected_markers = subset_tp_markers[subset_tp_markers$cluster == 2,]
        cluster_markers = c(cluster_markers, selected_markers$gene[1:100])
        selected_markers = subset_tp_markers[subset_tp_markers$cluster == 3,]
        cluster_markers = c(cluster_markers, selected_markers$gene[1:100])
        selected_markers = subset_tp_markers[subset_tp_markers$cluster == 4,]
        cluster_markers = c(cluster_markers, selected_markers$gene[1:100])
    }
}
cluster_markers = unique(cluster_markers)
# load gene expression
data_dir = file.path(base_dir, 'data/spatial_scRNAseq_GSE214989')
h5_file = file.path(data_dir, "GSE214989_counts_embryo_visium.h5")
data = Read10X_h5(h5_file)
common_cells = intersect(colnames(data), rownames(adata$obs))
common_genes = intersect(rownames(data), cluster_markers)
data = data[common_genes, common_cells]
adata = adata[common_cells, ]
gene_metadata = data.frame(rownames(data))
colnames(gene_metadata) = c('gene_short_name')
rownames(gene_metadata) = gene_metadata$gene_short_name
cds = new_cell_data_set(data, cell_metadata=adata$obs, gene_metadata=gene_metadata)
cds <- cds[,Matrix::colSums(exprs(cds)) != 0]
cds <- estimate_size_factors(cds)
cds = preprocess_cds(cds)
cds = reduce_dimension(cds, reduction_method='PCA', preprocess_method='PCA')
# cds = reduce_dimension(cds, reduction_method='UMAP', preprocess_method='PCA')
cds = cluster_cells(cds, reduction_method='PCA')
# cds = cluster_cells(cds, reduction_method='UMAP')
cds = learn_graph(cds)
plot_cells(cds, reduction_method='PCA', color_cells_by = "timepoint_original", cell_size=1)
ggsave(file.path(downstream_dir, 'monocle3_t0t1_latent_timepoint.png'), width=5, height=5)
plot_cells(cds, reduction_method='PCA', color_cells_by = "pred_labels", cell_size=1)
ggsave(file.path(downstream_dir, 'monocle3_t0t1_latent_predlabels.png'), width=5, height=5)



total_idx = c(t1_cluster3_idx, t2_cluster0_idx, t2_cluster2_idx)
adata = latent_adata[total_idx, ]
cluster_markers = c()
for (tp in c('t1', 't2')) {
    tp_markers = read.csv(file.path(downstream_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')),
                          header=TRUE, row.names=1)
    subset_tp_markers = tp_markers[tp_markers$avg_log2FC > 1 & tp_markers$p_val_adj < 0.05,]
    if (tp == 't1') {
        selected_markers = subset_tp_markers[subset_tp_markers$cluster == 3,]
        cluster_markers = c(cluster_markers, selected_markers$gene[1:100])
    }
    if (tp == 't2') {
        selected_markers = subset_tp_markers[subset_tp_markers$cluster == 0,]
        cluster_markers = c(cluster_markers, selected_markers$gene[1:100])
        selected_markers = subset_tp_markers[subset_tp_markers$cluster == 2,]
        cluster_markers = c(cluster_markers, selected_markers$gene[1:100])
    }
}
cluster_markers = unique(cluster_markers)
data_dir = file.path(base_dir, 'data/spatial_scRNAseq_GSE214989')
h5_file = file.path(data_dir, "GSE214989_counts_embryo_visium.h5")
data = Read10X_h5(h5_file)
common_cells = intersect(colnames(data), rownames(adata$obs))
common_genes = intersect(rownames(data), cluster_markers)
data = data[common_genes, common_cells]
adata = adata[common_cells, ]
gene_metadata = data.frame(rownames(data))
colnames(gene_metadata) = c('gene_short_name')
rownames(gene_metadata) = gene_metadata$gene_short_name
cds = new_cell_data_set(data, cell_metadata=adata$obs, gene_metadata=gene_metadata)
cds <- cds[,Matrix::colSums(exprs(cds)) != 0]
cds <- estimate_size_factors(cds)
cds = preprocess_cds(cds)
cds = reduce_dimension(cds, reduction_method='PCA', preprocess_method='PCA')
cds = reduce_dimension(cds, reduction_method='UMAP', preprocess_method='PCA')
cds = cluster_cells(cds, reduction_method='PCA')
cds = cluster_cells(cds, reduction_method='UMAP')
cds = learn_graph(cds)
plot_cells(cds, reduction_method='PCA', color_cells_by = "timepoint_original", cell_size=1)
ggsave(file.path(downstream_dir, 'monocle3_t1t2_latent_timepoint.png'), width=5, height=5)
plot_cells(cds, reduction_method='PCA', color_cells_by = "pred_labels", cell_size=1)
ggsave(file.path(downstream_dir, 'monocle3_t1t2_latent_predlabels.png'), width=5, height=5)
plot_cells(cds, reduction_method='UMAP', color_cells_by = "timepoint_original", cell_size=1)
ggsave(file.path(downstream_dir, 'monocle3_t1t2_latent_timepoint_UMAP.png'), width=5, height=5)
plot_cells(cds, reduction_method='UMAP', color_cells_by = "pred_labels", cell_size=1)
ggsave(file.path(downstream_dir, 'monocle3_t1t2_latent_predlabels_UMAP.png'), width=5, height=5)
