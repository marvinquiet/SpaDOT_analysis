Sys.setenv(RETICULATE_PYTHON = "/media/sda/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(dplyr)
library(stringr)
library(NbClust)
library(factoextra)
source('MOSTA_utils.R')

base_dir = "/media/sda/wenjinma"
result_dir = file.path(base_dir, "projects/spatialATAC-explore")
spatial_layout_dir = file.path(result_dir, "summarize_benchmark/MOSTA_brain_figures/spatial_layout")
dir.create(spatial_layout_dir, showWarnings = FALSE, recursive = TRUE)
downstream_dir = file.path(result_dir, "summarize_benchmark/MOSTA_brain_figures/MOSTA_brain_downstream")
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)
latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'MOSTA_brain_results', 'final_latent.h5ad'))

# --- explore number of clusters
# library(NbClust)
# for (tp in unique(latent_adata$obs$timepoint)) {
#     tp_latent_adata = latent_adata[latent_adata$obs$timepoint == tp, ]
#     set.seed(42)
#     g = fviz_nbclust(tp_latent_adata$X, kmeans, method = "gap_stat", k.max=20)
#     ggsave(file.path(spatial_layout_dir, paste0(tp, '_KMeans_elbow_gapstat.png')), plot=g, height=10, width=10)
#     g = fviz_nbclust(tp_latent_adata$X, kmeans, method = "silhouette", k.max=20)
#     ggsave(file.path(spatial_layout_dir, paste0(tp, '_KMeans_elbow_silhouette.png')), plot=g, height=10, width=10)
#     g = fviz_nbclust(tp_latent_adata$X, kmeans, method = "wss", k.max=20)
#     ggsave(file.path(spatial_layout_dir, paste0(tp, '_KMeans_elbow_wss.png')), plot=g, height=10, width=10)
# }

# --- explore wss
# Create a data frame for plotting
for (tp in unique(latent_adata$obs$timepoint)) {
    tp_latent_adata = latent_adata[latent_adata$obs$timepoint == tp, ]
    wss = sapply(5:20, function(k) {
        set.seed(42)
        kmeans(tp_latent_adata$X, centers = k, nstart = 10)$tot.withinss
    })
    wss_df = data.frame(
        k = 5:20,
        wss = wss
    )
    wss_df$wss_diff = c(-(diff(wss_df$wss)), NA)
    wss_df$wss_ratio = c(NA, wss_df$wss_diff[-nrow(wss_df)] / wss_df$wss_diff[-1])
    ggplot(wss_df, aes(x = k, y = wss)) +
        geom_point() +
        geom_line() +
        ggtitle("Elbow Method for K Selection") +
        xlab("Number of Clusters (K)") +
        ylab("Within-Cluster Sum of Squares") + 
        theme_pubr()
    ggsave(file.path(spatial_layout_dir, paste0(tp, '_KMeans_elbow_wss_myself.pdf')), height=5, width=5)  
}
# --- in the end, there are 7,8,10,9,10

data_dir = file.path('/media/sda/wenjinma/data/spatial_MOSTA')
original_adata = read_h5ad(file.path(data_dir, 'adata_brain_clean.h5ad'))
common_cells = intersect(rownames(original_adata), rownames(latent_adata))
adata = original_adata[common_cells]
latent_adata = latent_adata[common_cells, ]
# --- create Seurat object for analysis
Seurat_obj = CreateSeuratObject(counts=t(adata$X), meta.data=latent_adata$obs)
Seurat_obj[['spaDOT_latent']] = CreateAssayObject(data = t(latent_adata$X))
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# --- plot gene expression 
marker_genes = c('Fabp7', 'Lhx2', 'Tfdp1', 'Adgrv1', 'Dbx1', 'Lmo1', 'Slc1a3', 'Ttyh1', 'Neurod6', 'Bhlhe22', 'Tbr1', 'Islr2', 'Dlx1', 'Dlx5', 'Zdbf2', 'Hap1', 'Tcf7l2', 'Lhx9', 'Tfap2b', 'Pou4f1', 'Ina', 'Snhg11', 'Ebf2', 'Ebf3', 'Hoxa5', 'Hoxb5', 'Prph', 'Sncg', 'Ttr', 'Clic6', 'Gad2', 'Hba-a1', 'Slc17a6', 'Foxp2', 'Snap25', 'Slc32a1', 'Aldh1a1', 'Rora', 'Rorb', 'Gbx2', 'Egr2', 'Lmx1b', 'Phox2b', 'Dbh', 'Mafb', 'Hoxa1', 'Hoxa2', 'Hoxb1')
for (gene in marker_genes) {
    gene_idx = which(rownames(Seurat_obj) == gene)
    gene_exprs = Seurat_obj@assays$RNA@layers$data[gene_idx, ]
    gene_exprs_df = Seurat_obj@meta.data
    gene_exprs_df[, gene] = gene_exprs   
    ggplot(gene_exprs_df, aes(x=pixel_x, y=pixel_y, color=.data[[gene]])) + geom_point(size=0.1) +
        facet_wrap(~timepoint, scale='free', nrow=1) + 
        theme_pubr() + 
        theme(axis.text =element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) + 
        scale_color_gradient(low = "grey", high = "red")
    ggsave(file.path(downstream_dir, paste0(gene, '_geneexprs.png')),height=3, width=12)    
}

# perform KMeans clustering 
# concatenate tp_latent_adata
# cluster_num = c('t0'=8, 't1'=11, 't2'=11, 't3'=9, 't4'=10)
cluster_num = c('t0'=7, 't1'=8, 't2'=10, 't3'=9, 't4'=10)
tp_latent_adata_list = list()
for (tp in unique(latent_adata$obs$timepoint)) {
    tp_latent_adata = latent_adata[latent_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_latent_adata$obs$pred_labels = kmeans(tp_latent_adata$X, centers=unname(cluster_num[tp]), nstart = 10)$cluster
    tp_latent_adata$obs$pred_labels = as.factor(tp_latent_adata$obs$pred_labels)
    tp_latent_adata_list[[tp]] = tp_latent_adata
}
pred_adata = concat(tp_latent_adata_list)
write_h5ad(pred_adata, file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'MOSTA_brain_results', 'prediction_adata.h5ad'))
write.csv(pred_adata$obs, file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'MOSTA_brain_results', 'spaDOT_pred.csv'))


# --- generate marker gene enrichment
marker_genes = c('Fabp7', 'Lhx2', 'Tfdp1', 'Adgrv1', 'Dbx1', 'Lmo1', 'Slc1a3', 'Ttyh1', 'Neurod6', 'Bhlhe22', 'Tbr1', 'Islr2', 'Dlx1', 'Dlx5', 'Zdbf2', 'Hap1', 'Tcf7l2', 'Lhx9', 'Tfap2b', 'Pou4f1', 'Ina', 'Snhg11', 'Ebf2', 'Ebf3', 'Hoxa5', 'Hoxb5', 'Prph', 'Sncg', 'Ttr', 'Clic6')
names(marker_genes) = rep('Brain', length(marker_genes))
Seurat_obj$spaDOT_pred_labels = pred_adata$obs[colnames(Seurat_obj), 'pred_labels']
# --- plot gene expression in each domain
for (tp in unique(pred_adata$obs$timepoint)) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    Idents(tp_Seurat_obj) = tp_Seurat_obj$spaDOT_pred_labels
    g = DotPlot(tp_Seurat_obj, features=unname(marker_genes), dot.scale = 6) + RotatedAxis() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank())
    ggsave(file.path(downstream_dir, paste0(tp, '_marker_genes.pdf')), plot=g, height=5, width=10)
}
# --- calculate enrichment
# enrichment_df = generate_marker_enrichment(Seurat_obj, 'spaDOT_pred_labels', marker_genes)
# write.csv(enrichment_df, file.path(downstream_dir, 'marker_enrichment.csv'), row.names = FALSE)

# --- plot enrichment results
enrichment_df = read.csv(file.path(downstream_dir, 'marker_enrichment.csv'))
enrichment_df$enrichment = as.numeric(enrichment_df$enrichment)
enrichment_df$gene = factor(enrichment_df$gene, levels=unname(marker_genes))
enrichment_df$cluster = as.character(enrichment_df$cluster)
id_to_tp_vec = c('t0'='E12.5', 't1'='E13.5', 't2'='E14.5', 't3'='E15.5', 't4'='E16.5')
for (tp in unique(enrichment_df$timepoint)) {
    tp_enrichment_df = enrichment_df[enrichment_df$timepoint == tp, ]
    tp_enrichment_df = tp_enrichment_df[tp_enrichment_df$enrichment > 0, ]
    ggplot(data=tp_enrichment_df, aes(x=gene, y=log10(enrichment), group=cluster)) +
        geom_point(aes(color=cluster)) +
        geom_hline(yintercept=0, linetype="dashed", color = "red") +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    ggsave(file.path(downstream_dir, paste0(tp, '_marker_enrichment.pdf')), height=4, width=6)
}

# --- plot original labels
for (tp in unique(pred_adata$obs$timepoint)) {
    tp_latent_adata = pred_adata[pred_adata$obs$timepoint == tp, ]
    p = ggplot(tp_latent_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(spatial_layout_dir, paste0(tp, '_KMeans_original.pdf')), plot=p, height=5, width=5)
}
# --- plot original transport matrix
timepoints = as.vector(unique(pred_adata$obs$timepoint))
timepoints = sort(timepoints)
plot_OT_heatmap(timepoints, 'spaDOT', file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'MOSTA_brain_results'))


# --- change color palette
tp_vector = c('t0'='E12.5', 't1'='E13.5', 't2'='E14.5', 't3'='E15.5', 't4'='E16.5')
pred_adata$obs$original_timepoint = unname(tp_vector[pred_adata$obs$timepoint])
# E12.5: 1(Mb VZ+Hb VZ), 2(SPall VZ), 3 (Pall VZ+Pall), 4(Thalamus), 5(FMN), 6(Mb+Hb), 7(SPall)
# E13.5: 1(Pall+OB), 2(Thalamus), 3(Mb), 4(Pall VZ), 5(HB), 6(Hb VZ), 7(SPall VZ+MB VZ), 8(SPall)
# E14.5: 1(Mb VZ), 2(SPall+OB), 3(Pall+OB), 4(spC), 5(Hypo), 6 (Mb), 7(Hb), 8 (Hb VZ), 9(chP), 10(Pall VZ)
# E15.5: 1(Hb+Mb), 2(Thalamus), 3(Hb VZ), 4(Pall VZ), 5(Chp+Mb VZ), 6 (Cere + Mb), 7 (spC+FMN), 8(SPall), 9(Pall)
# E16.5: 1(Mb), 2(Hb), 3(Hb Vz), 4(Spall1), 5 (Cere), 6 (SPall 2), 7(FMN+spC), 8(Thalamus), 9(Pall), 10(Mb VZ)

# --- correspondance
# Mb+Hb=1, Mb=1-1, Hb=1-2, MbVZ+HbVZ=2, Mb VZ=2-1, Hb VZ=2-2, PallVZ+Pall=3, PallVZ=3-1, Pall=3-2, SPall=4, 
# SPall VZ=5, spC=6, chP=7, FMN=8, Hypo=9, Cere=10, Thalamus=11

cluster_mapping = list(
    'E12.5'=c("1"="2", "2"="5", "3"="3", "4"="11", "5"="8", "6"="1", "7"="4"),
    'E13.5'=c("1"="3-2", "2"="11", "3"="1-1", "4"="3-1", "5"="1-2", "6"="2-2", "7"="2-1", "8"="4"),
    'E14.5'=c("1"="2-1", "2"="4", "3"="3-2", "4"="6", "5"="9", "6"="1-1", "7"="1-2", "8"="2-2", "9"="7", "10"="3-1"),
    'E15.5'=c("1"="1", "2"="11", "3"="2-2", "4"="3-1", "5"="2-1", "6"="10", "7"="6", "8"="4", "9"="3-2"),
    'E16.5'=c("1"="1-1", "2"="1-2", "3"="2-2", "4"="4-1", "5"="10", "6"="4-2", "7"="6", "8"="11", "9"="3-2", "10"="2-1")
)
pred_adata$obs = pred_adata$obs %>% mutate(
    spaDOT_pred_labels = case_when(
        timepoint_original == 'E12.5' ~ cluster_mapping[['E12.5']][as.character(pred_labels)],
        timepoint_original == 'E13.5' ~ cluster_mapping[['E13.5']][as.character(pred_labels)],
        timepoint_original == 'E14.5' ~ cluster_mapping[['E14.5']][as.character(pred_labels)],
        timepoint_original == 'E15.5' ~ cluster_mapping[['E15.5']][as.character(pred_labels)],
        timepoint_original == 'E16.5' ~ cluster_mapping[['E16.5']][as.character(pred_labels)]
    )
)
pred_adata$obs$spaDOT_pred_labels = factor(pred_adata$obs$spaDOT_pred_labels,
levels=c("1", "1-1", "1-2", "2", "2-1", "2-2", "3", "3-1", "3-2", "4", "4-1", "4-2",
"5", "6", "7", "8", "9", "10", "11"))

tp_color_palette = c("1"="#cdb4db", "1-1"="#ffafcc", "1-2"="#a2d2ff",
"2"="#8338ec", "2-1"="#ff006e", "2-2"="#3a86ff",
"3"="#38b000", "3-1"="#004b23", "3-2"="#ccff33",
"4"="#e85d04", "5"="#ffba08", "6"="#d9d9d9", "7"="#7f5539", "8"="#a89f68", "9"="#227c9d", "11"="#17c3b2", "10"="#41521f",
"4-1"="#fcf300", "4-2"="#ff0000")
plot_pred_label(pred_adata, color_col='spaDOT_pred_labels', file_name=file.path(spatial_layout_dir, 'spaDOT_prediction.pdf'), tp_color_palette=tp_color_palette, legend_ncol=4, height=3.5, width=15)


# E12.5: 1 (SPall), 2 (Mb), 3 (big region: HB??), 4 (Pall), 5 (FMN), 6 (SPall VZ), 7 (Thalamus), 8 (Mb Vz)
# E13.5: 1 (Hb VZ), 2 (spC+FMN), 3 (Cere), 4 (Mb VZ), 5 (SPall), 6 (SPall VZ), 7 (Mb), 8 (Pall), 9 (Thalamus), 10 (big region: HB??), 11 (Pall VZ)
# E14.5: 1 (Mb VZ), 2 (Pall VZ), 3 (Pall), 4 (Hb VZ), 5 (Mb), 6 (Chor), 7 (big region: HB??), 8 (Cavity), 9 (spC), 10 (Hy), 11(SPall)
# E15.5: 1 (big region: HB??), 2 (Thalamus), 3 (Hb VZ), 4 (Pall VZ), 5 (Mb VZ+FMN), 6 (Cere), 7 (spC), 8 (SPall), 9 (Pall)
# E16.5 : 1 (Mb), 2 (HB), 3 (Hb VZ), 4 (SPall), 5 (Cere+Chor), 6 (HY), 7 (sPC+FMN), 8 (Thalamus), 9 (Pall), 10 (Mb VZ)

# --- assign colors
# SPall: #b67c3b, Spall VZ: #AF4E24, Pall: #C6CDF7, Pall VZ: #7294D4, Mb: #A2A475, Mb Vz: #02401B, Thalamus: #5B1A18, HY: #FD6467, FMN: #0F0D0E, Cavity: #B0AFA2, spC: #D9D0D3, sPC+FMN: #8D8680, Cere: #6C8645, HB: #F1BB7B, Hb VZ: #D67236, Chor: #D8B70A


# tp_color_palette = list(
#     'E12.5'=c('1'='#b67c3b', '2'='#A2A475', '3'='#F1BB7B', '4'='#C6CDF7', '5'='#0F0D0E', '6'='#AF4E24', '7'='#5B1A18', '8'='#02401B'),
#     'E13.5'=c('1'='#D67236', '2'='#8D8680', '3'='#6C8645', '4'='#02401B', '5'='#b67c3b', '6'='#AF4E24', '7'='#A2A475', '8'='#C6CDF7', '9'='#5B1A18', '10'='#F1BB7B', '11'='#7294D4'),
#     'E14.5'=c('1'='#02401B', '2'='#7294D4', '3'='#C6CDF7', '4'='#D67236', '5'='#A2A475', '6'='#D8B70A', '7'='#F1BB7B', '8'='#B0AFA2', '9'='#D9D0D3', '10'='#FD6467', '11'='#b67c3b'),
#     'E15.5'=c('1'='#F1BB7B', '2'='#5B1A18', '3'='#D67236', '4'='#7294D4', '5'='#02401B', '6'='#6C8645', '7'='#D9D0D3', '8'='#b67c3b', '9'='#C6CDF7'),
#     'E16.5'=c('1'='#A2A475', '2'='#F1BB7B', '3'='#D67236', '4'='#b67c3b', '5'='#6C8645', '6'='#FD6467', '7'='#8D8680', '8'='#5B1A18', '9'='#C6CDF7', '10'='#02401B')
# )

# for (tp in unique(pred_adata$obs$original_timepoint)) {
#     tp_latent_adata = pred_adata[pred_adata$obs$original_timepoint == tp, ]
#     p = ggplot(tp_latent_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
#         geom_point(size=0.5) +
#         guides(color = guide_legend(ncol=4, title=NULL)) + 
#         theme_void() + 
#         ggtitle(tp) + 
#         theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
#     p = p + scale_color_manual(values=tp_color_palette[[tp]])
#     ggsave(file.path(spatial_layout_dir, paste0(tp, '_KMeans.png')), plot=p, height=5, width=5)
# }

timepoints = as.vector(unique(pred_adata$obs$timepoint))
timepoints = sort(timepoints)
plot_OT_heatmap(timepoints, 'spaDOT', file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'MOSTA_brain_results'), tp_color_palette=tp_color_palette)

# --- plot gene enrichment analysis
enrichment_df = read.csv(file.path(downstream_dir, 'marker_enrichment.csv'))
enrichment_df = enrichment_df %>% mutate(
    spaDOT_pred_labels = case_when(
        timepoint == 't0' ~ cluster_mapping[['E12.5']][as.character(cluster)],
        timepoint == 't1' ~ cluster_mapping[['E13.5']][as.character(cluster)],
        timepoint == 't2' ~ cluster_mapping[['E14.5']][as.character(cluster)],
        timepoint == 't3' ~ cluster_mapping[['E15.5']][as.character(cluster)],
        timepoint == 't4' ~ cluster_mapping[['E16.5']][as.character(cluster)]
    )
)
enrichment_df$spaDOT_pred_labels = factor(enrichment_df$spaDOT_pred_labels, levels=c("1", "1-1", "1-2", "2", "2-1", "2-2", "3", "3-1", "3-2", "4", "4-1", "4-2",
    "5", "6", "7", "8", "9", "10", "11"))
selected_markers = c('Tcf7l2', 'Lhx9', 'Tfap2b', 'Dbx1', 'Slc1a3', 'Ttyh1',
'Lhx2', 'Neurod6', 'Tbr1', 'Dlx1', 'Dlx5', 'Adgrv1', 'Hoxa5', 'Hoxb5', 'Ttr', 'Clic6',
'Prph', 'Sncg', 'Zdbf2', 'Hap1', 'Ebf2')
sub_enrichment_df = enrichment_df[enrichment_df$gene %in% selected_markers, ]
sub_enrichment_df$gene = factor(sub_enrichment_df$gene, levels=selected_markers)
# # Mb+Hb=1, Mb=1-1, Hb=1-2, MbVZ+HbVZ=2, Mb VZ=2-1, Hb VZ=2-2, PallVZ+Pall=3, PallVZ=3-1, Pall=3-2, SPall=4, 
# SPall VZ=5, spC=6, chP=7, FMN=8, Hypo=9, Cere=10, Thalamus=11
for (tp in unique(sub_enrichment_df$timepoint)) {
    tp_enrichment_df = sub_enrichment_df[sub_enrichment_df$timepoint == tp, ]
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    tp_Seurat_counts = tp_Seurat_obj@assays$RNA@layers$counts
    rownames(tp_Seurat_counts) = rownames(tp_Seurat_obj)
    colnames(tp_Seurat_counts) = colnames(tp_Seurat_obj)
    tp_Seurat_counts = tp_Seurat_counts[selected_markers,]
    # --- obtain gene expresson percentage per cluster for enriched genes
    tp_Seurat_counts = as.data.frame(t(tp_Seurat_counts))
    tp_Seurat_counts$barcode = rownames(tp_Seurat_counts)
    tp_Seurat_counts$spaDOT_pred_labels = unname(cluster_mapping[[unname(id_to_tp_vec[tp])]][as.character(tp_Seurat_obj$spaDOT_pred_labels)])
    tp_expr_pct = tp_Seurat_counts %>%
        group_by(spaDOT_pred_labels) %>%
        summarise(across(
            .cols = where(is.numeric),
            .fns = ~ mean(. > 0) * 100
        )) %>% as.data.frame()
    tp_expr_pct_long = reshape2::melt(tp_expr_pct, id.vars = "spaDOT_pred_labels")
    colnames(tp_expr_pct_long) = c('spaDOT_pred_labels', 'gene', 'pct')
    tp_enrichment_df = merge(tp_enrichment_df, tp_expr_pct_long, by=c('spaDOT_pred_labels', 'gene'))
    # --- preprocess enrichment score, keep those larger than 0
    tp_enrichment_df$processed_enrichment = as.numeric(tp_enrichment_df$enrichment)
    tp_enrichment_df$processed_enrichment[tp_enrichment_df$processed_enrichment < 0] = 0
    tp_enrichment_df$processed_enrichment = log10(tp_enrichment_df$processed_enrichment)
    # --- plot dotplot
    ggplot(data=tp_enrichment_df, aes(x=gene, y=spaDOT_pred_labels, size=pct, color=processed_enrichment)) +
        geom_point() +
        scale_color_gradient(low = "grey", high = "red", breaks=c(-2, -1, 0, 1, 2)) +
        scale_size_continuous(range = c(1, 6), breaks=c(25, 50, 75)) +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        labs(color='Enrichment', size='Ratio (%)', x="", y='SpaDOT Domains')
    ggsave(file.path(downstream_dir, paste0(tp, '_marker_enrichment_dotplot.pdf')), height=5, width=6)

    # tp_enrichment_df = tp_enrichment_df[tp_enrichment_df$enrichment > 0, ]
    # ggplot(data=tp_enrichment_df, aes(x=gene, y=log(enrichment), group=spaDOT_pred_labels)) +
    #     geom_point(aes(color=spaDOT_pred_labels)) +
    #     geom_hline(yintercept=0, linetype="dashed", color = "red") +
    #     scale_color_manual(values=tp_color_palette) +
    #     theme_pubr() +
    #     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    # ggsave(file.path(downstream_dir, paste0(tp, '_marker_enrichment.pdf')), height=5, width=6)
}

# --- plot results for SEDR
SEDR_adata = read_h5ad(file.path(result_dir, 'explore_SEDR', 'MOSTA_brain_results', 'SEDR_adata.h5ad'))
SEDR_adata$obs$pixel_y = -SEDR_adata$obs$pixel_y
# SEDR_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
SEDR_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
SEDR_adata_list = list()
for (tp in unique(SEDR_adata$obs$timepoint)) {
    tp_SEDR_adata = SEDR_adata[SEDR_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_SEDR_adata$obs$pred_labels = kmeans(tp_SEDR_adata$obsm$'SEDR.Harmony', centers=unname(SEDR_cluster_num[tp]), nstart = 10)$cluster
    tp_SEDR_adata$obs$pred_labels = as.factor(tp_SEDR_adata$obs$pred_labels)
    SEDR_adata_list[[tp]] = tp_SEDR_adata
    # plot
    p = ggplot(tp_SEDR_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_SEDR', 'MOSTA_brain_results', paste0(tp, '_SEDR_KMeans.png')), plot=p, height=5, width=5)
}
SEDR_pred_adata = concat(SEDR_adata_list)
write.csv(SEDR_pred_adata$obs, file.path(result_dir, 'explore_SEDR', 'MOSTA_brain_results', 'SEDR_pred.csv'))

# --- plot results for STAligner
STAligner_adata = read_h5ad(file.path(result_dir, 'explore_STAligner', 'MOSTA_brain_results', 'STAligner_adata.h5ad'))
STAligner_adata$obs$pixel_x = STAligner_adata$obsm$spatial[, 1]
STAligner_adata$obs$pixel_y = -STAligner_adata$obsm$spatial[, 2]
# STAligner_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
STAligner_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
STAligner_adata_list = list()
for (tp in unique(STAligner_adata$obs$timepoint)) {
    tp_STAligner_adata = STAligner_adata[STAligner_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_STAligner_adata$obs$pred_labels = kmeans(tp_STAligner_adata$obsm$'STAligner', centers=unname(STAligner_cluster_num[tp]), nstart = 10)$cluster
    tp_STAligner_adata$obs$pred_labels = as.factor(tp_STAligner_adata$obs$pred_labels)
    STAligner_adata_list[[tp]] = tp_STAligner_adata
    # plot
    p = ggplot(tp_STAligner_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_STAligner', 'MOSTA_brain_results', paste0(tp, '_STAligner_KMeans.png')), plot=p, height=5, width=5)
}
STAligner_pred_adata = concat(STAligner_adata_list)
write.csv(STAligner_pred_adata$obs, file.path(result_dir, 'explore_STAligner', 'MOSTA_brain_results', 'STAligner_pred.csv'))

# --- plot results for spatialPCA
spatialPCA_merged_obj = readRDS(file.path(result_dir, 'explore_spatialPCA', 'MOSTA_brain_results', 'spatial_merge_obj.RDS'))
spatialPCA_res = spatialPCA_merged_obj@SpatialPCs
common_cells = intersect(colnames(spatialPCA_res), original_adata$obs_names)
spatialPCA_metadata = original_adata$obs[common_cells, ]
spatialPCA_metadata$pixel_x = original_adata[common_cells, ]$obsm$spatial[, 1]
spatialPCA_metadata$pixel_y = -original_adata[common_cells, ]$obsm$spatial[, 2]
spatialPCA_res = spatialPCA_res[, common_cells]
# spatialPCA_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
spatialPCA_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
spatialPCA_pred_df = NULL
for (tp in unique(spatialPCA_metadata$timepoint)) {
    tp_idx = which(spatialPCA_metadata$timepoint == tp)
    tp_spatialPCA_res = spatialPCA_res[, tp_idx]
    set.seed(42)
    tp_pred_labels = kmeans(t(tp_spatialPCA_res), centers=unname(spatialPCA_cluster_num[tp]), nstart = 10)$cluster
    tp_spatialPCA_res_metadata = spatialPCA_metadata[tp_idx, ]
    tp_spatialPCA_res_metadata$pred_labels = as.factor(tp_pred_labels)
    if (is.null(spatialPCA_pred_df)) {
        spatialPCA_pred_df = tp_spatialPCA_res_metadata
    } else {
        spatialPCA_pred_df = rbind(spatialPCA_pred_df, tp_spatialPCA_res_metadata)
    }
    # plot
    p = ggplot(tp_spatialPCA_res_metadata, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_spatialPCA', 'MOSTA_brain_results', paste0(tp, '_spatialPCA_KMeans.png')), plot=p, height=5, width=5)
}
write.csv(spatialPCA_pred_df, file.path(result_dir, 'explore_spatialPCA', 'MOSTA_brain_results', 'spatialPCA_pred.csv'))

marker_genes = c('Fabp7', 'Lhx2', 'Tfdp1', 'Adgrv1', 'Dbx1', 'Lmo1', 'Slc1a3', 'Ttyh1', 'Neurod6', 'Bhlhe22', 'Tbr1', 'Islr2', 'Dlx1', 'Dlx5', 'Zdbf2', 'Hap1', 'Tcf7l2', 'Lhx9', 'Tfap2b', 'Pou4f1', 'Ina', 'Snhg11', 'Ebf2', 'Ebf3', 'Hoxa5', 'Hoxb5', 'Prph', 'Sncg', 'Ttr', 'Clic6')
names(marker_genes) = rep('Brain', length(marker_genes))
Seurat_obj$SpatialPCA_pred_labels = spatialPCA_pred_df[colnames(Seurat_obj), 'pred_labels']
enrichment_df = generate_marker_enrichment(Seurat_obj, 'SpatialPCA_pred_labels', marker_genes)
write.csv(enrichment_df, file.path(downstream_dir, 'SpatialPCA_marker_enrichment.csv'), row.names = FALSE)

enrichment_df$enrichment = as.numeric(enrichment_df$enrichment)
enrichment_df$gene = factor(enrichment_df$gene, levels=unname(marker_genes))
for (tp in unique(enrichment_df$timepoint)) {
    tp_enrichment_df = enrichment_df[enrichment_df$timepoint == tp, ]
    tp_enrichment_df = tp_enrichment_df[tp_enrichment_df$enrichment > 0, ]
    ggplot(data=tp_enrichment_df, aes(x=gene, y=log10(enrichment), group=cluster)) +
        geom_point(aes(color=cluster)) +
        geom_hline(yintercept=0, linetype="dashed", color = "red") +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    ggsave(file.path(downstream_dir, paste0(tp, '_marker_enrichment_SpatialPCA.png')), height=4, width=4)
}

# --- plot results for spatiAlign
spatiAlign_adata = read_h5ad(file.path(result_dir, 'explore_spatiAlign', 'MOSTA_brain_results', 'spatiAlign_adata.h5ad'))
spatiAlign_adata$obs$pixel_x = spatiAlign_adata$obsm$spatial[, 1]
spatiAlign_adata$obs$pixel_y = -spatiAlign_adata$obsm$spatial[, 2]
# spatiAlign_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
spatiAlign_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
spatiAlign_adata_list = list()
for (tp in unique(spatiAlign_adata$obs$timepoint)) {
    tp_spatiAlign_adata = spatiAlign_adata[spatiAlign_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_spatiAlign_adata$obs$pred_labels = kmeans(tp_spatiAlign_adata$obsm$'correct', centers=unname(spatiAlign_cluster_num[tp]), nstart = 10)$cluster
    tp_spatiAlign_adata$obs$pred_labels = as.factor(tp_spatiAlign_adata$obs$pred_labels)
    spatiAlign_adata_list[[tp]] = tp_spatiAlign_adata
    # plot
    p = ggplot(tp_spatiAlign_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_spatiAlign', 'MOSTA_brain_results', paste0(tp, '_spatiAlign_KMeans.png')), plot=p, height=5, width=5)
}
spatiAlign_pred_adata = concat(spatiAlign_adata_list)
write.csv(spatiAlign_pred_adata$obs, file.path(result_dir, 'explore_spatiAlign', 'MOSTA_brain_results', 'spatiAlign_pred.csv'))

# --- spiral 
SPIRAL_adata = read_h5ad(file.path(result_dir, 'explore_SPIRAL', 'MOSTA_brain_results', 'SPIRAL_corrected.h5ad'))
SPIRAL_adata$obs$pixel_x = original_adata[SPIRAL_adata$obs_names, ]$obsm$spatial[, 1]
SPIRAL_adata$obs$pixel_y = -original_adata[SPIRAL_adata$obs_names, ]$obsm$spatial[, 2]
# SPIRAL_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
SPIRAL_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
SPIRAL_adata_list = list()
for (tp in unique(SPIRAL_adata$obs$timepoint)) {
    tp_SPIRAL_adata = SPIRAL_adata[SPIRAL_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_SPIRAL_adata$obs$pred_labels = kmeans(tp_SPIRAL_adata$obsm$'spiral', centers=unname(SPIRAL_cluster_num[tp]), nstart = 10)$cluster
    tp_SPIRAL_adata$obs$pred_labels = as.factor(tp_SPIRAL_adata$obs$pred_labels)
    SPIRAL_adata_list[[tp]] = tp_SPIRAL_adata
    # plot
    p = ggplot(tp_SPIRAL_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_SPIRAL', 'MOSTA_brain_results', paste0(tp, '_SPIRAL_KMeans.png')), plot=p, height=5, width=5)
}
SPIRAL_pred_adata = concat(SPIRAL_adata_list)
write.csv(SPIRAL_pred_adata$obs, file.path(result_dir, 'explore_SPIRAL', 'MOSTA_brain_results', 'SPIRAL_pred.csv'))

# --- analyze graphST
GraphST_adata = read_h5ad(file.path(result_dir, 'explore_GraphST', 'MOSTA_brain_results', 'GraphST_adata.h5ad'))
GraphST_adata$obs$pixel_x = GraphST_adata$obsm$spatial[, 1]
GraphST_adata$obs$pixel_y = -GraphST_adata$obsm$spatial[, 2]
# GraphST_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
GraphST_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
GraphST_adata_list = list()
for (tp in unique(GraphST_adata$obs$timepoint)) {
    tp_GraphST_adata = GraphST_adata[GraphST_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_GraphST_adata$obs$pred_labels = kmeans(tp_GraphST_adata$obsm$'emb', centers=unname(GraphST_cluster_num[tp]), nstart = 10)$cluster
    tp_GraphST_adata$obs$pred_labels = as.factor(tp_GraphST_adata$obs$pred_labels)
    GraphST_adata_list[[tp]] = tp_GraphST_adata
    # plot
    p = ggplot(tp_GraphST_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_GraphST', 'MOSTA_brain_results', paste0(tp, '_GraphST_KMeans.png')), plot=p, height=5, width=5)
}
GraphST_pred_adata = concat(GraphST_adata_list)
write.csv(GraphST_pred_adata$obs, file.path(result_dir, 'explore_GraphST', 'MOSTA_brain_results', 'GraphST_pred.csv'))

# --- analyze spaGCN
spaGCN_adata = read_h5ad(file.path(result_dir, 'explore_spaGCN', 'MOSTA_brain_results', 'spaGCN.h5ad'))
spaGCN_adata$obs$pixel_x = spaGCN_adata$obsm$spatial[, 1]
spaGCN_adata$obs$pixel_y = -spaGCN_adata$obsm$spatial[, 2]
# spaGCN_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
spaGCN_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
spaGCN_adata_list = list()
for (tp in unique(spaGCN_adata$obs$timepoint)) {
    tp_spaGCN_adata = spaGCN_adata[spaGCN_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_spaGCN_adata$obs$pred_labels = kmeans(tp_spaGCN_adata$obsm$'spaGCN_embed', centers=unname(spaGCN_cluster_num[tp]), nstart = 10)$cluster
    tp_spaGCN_adata$obs$pred_labels = as.factor(tp_spaGCN_adata$obs$pred_labels)
    spaGCN_adata_list[[tp]] = tp_spaGCN_adata
    # plot
    p = ggplot(tp_spaGCN_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_spaGCN', 'MOSTA_brain_results', paste0(tp, '_spaGCN_KMeans.png')), plot=p, height=5, width=5)
}
spaGCN_pred_adata = concat(spaGCN_adata_list)
write.csv(spaGCN_pred_adata$obs, file.path(result_dir, 'explore_spaGCN', 'MOSTA_brain_results', 'spaGCN_pred.csv'))


# --- analyze PRECAST
PRECAST_seuratobj = readRDS(file.path(result_dir, 'explore_PRECAST', 'MOSTA_brain_results', 'PRECAST_res.RDS'))
# PRECAST_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
PRECAST_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
PRECAST_seuratobj$pixel_y = -PRECAST_seuratobj$pixel_y
PRECAST_pred_df = NULL
for (tp in unique(PRECAST_seuratobj$timepoint)) {
    tp_PRECAST_res = subset(PRECAST_seuratobj, subset=timepoint == tp)
    set.seed(42)
    tp_pred_labels = kmeans(Embeddings(tp_PRECAST_res, reduction='PRECAST'), centers=unname(PRECAST_cluster_num[tp]), nstart = 10)$cluster
    tp_PRECAST_res$pred_labels = as.factor(tp_pred_labels)
    # plot
    if (is.null(PRECAST_pred_df)) {
        PRECAST_pred_df = tp_PRECAST_res@meta.data
    } else {
        PRECAST_pred_df = rbind(PRECAST_pred_df, tp_PRECAST_res@meta.data)
    }
    p = ggplot(tp_PRECAST_res@meta.data, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_PRECAST', 'MOSTA_brain_results', paste0(tp, '_PRECAST_KMeans.png')), plot=p, height=5, width=5)
}
write.csv(PRECAST_pred_df, file.path(result_dir, 'explore_PRECAST', 'MOSTA_brain_results', 'PRECAST_pred.csv'))

# --- analyze STADIA
STADIA_res = readRDS(file.path(result_dir, 'explore_STADIA', 'MOSTA_brain_results', 'stadia.RDS'))
STADIA_adata = original_adata[colnames(STADIA_res$factors), ]
STADIA_factors = t(STADIA_res$factors)
rownames(STADIA_factors) = NULL
STADIA_adata$obsm$stadia = t(STADIA_res$factors)
STADIA_adata$obs$pixel_x = STADIA_adata$obsm$spatial[, 1]
STADIA_adata$obs$pixel_y = -STADIA_adata$obsm$spatial[, 2]
# STADIA_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
STADIA_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
STADIA_adata_list = list()
for (tp in unique(STADIA_adata$obs$timepoint)) {
    tp_STADIA_adata = STADIA_adata[STADIA_adata$obs$timepoint == tp, ]
    set.seed(42)
    tp_STADIA_adata$obs$pred_labels = kmeans(tp_STADIA_adata$obsm$'stadia', centers=unname(STADIA_cluster_num[tp]), nstart = 10)$cluster
    tp_STADIA_adata$obs$pred_labels = as.factor(tp_STADIA_adata$obs$pred_labels)
    STADIA_adata_list[[tp]] = tp_STADIA_adata
    # plot
    p = ggplot(tp_STADIA_adata$obs, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_STADIA', 'MOSTA_brain_results', paste0(tp, '_STADIA_KMeans.png')), plot=p, height=5, width=5)
}
STADIA_pred_adata = concat(STADIA_adata_list)
write.csv(STADIA_pred_adata$obs, file.path(result_dir, 'explore_STADIA', 'MOSTA_brain_results', 'STADIA_pred.csv'))

# --- analyze spaVAE
spaVAE_latent = read.csv(file.path(result_dir, 'explore_spaVAE', 'MOSTA_brain_results', 'final_latent.txt'), header=FALSE)
spaVAE_metadata = read.csv(file.path(result_dir, 'explore_spaVAE', 'MOSTA_brain_results', 'combined_metadata.csv'), header=TRUE, row.names=1)
spaVAE_metadata$pixel_x = original_adata[rownames(spaVAE_metadata), ]$obsm$spatial[, 1]
spaVAE_metadata$pixel_y = -original_adata[rownames(spaVAE_metadata), ]$obsm$spatial[, 2]
# spaVAE_cluster_num = c('E12.5'=8, 'E13.5'=11, 'E14.5'=11, 'E15.5'=9, 'E16.5'=10)
spaVAE_cluster_num = c('E12.5'=7, 'E13.5'=8, 'E14.5'=10, 'E15.5'=9, 'E16.5'=10)
spaVAE_pred_df = NULL
for (tp in unique(spaVAE_metadata$timepoint)) {
    tp_idx = which(spaVAE_metadata$timepoint == tp)
    tp_spaVAE_latent = spaVAE_latent[tp_idx, ]
    set.seed(42)
    tp_pred_labels = kmeans(tp_spaVAE_latent, centers=unname(spaVAE_cluster_num[tp]), nstart = 10)$cluster
    tp_spaVAE_metadata = spaVAE_metadata[tp_idx, ]
    tp_spaVAE_metadata$pred_labels = as.factor(tp_pred_labels)
    if (is.null(spaVAE_pred_df)) {
        spaVAE_pred_df = tp_spaVAE_metadata
    } else {
        spaVAE_pred_df = rbind(spaVAE_pred_df, tp_spaVAE_metadata)
    }
    # plot
    p = ggplot(tp_spaVAE_metadata, aes(x=pixel_x, y=pixel_y, color=pred_labels)) + 
        geom_point(size=0.5) +
        guides(color = guide_legend(ncol=4, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
    ggsave(file.path(result_dir, 'explore_spaVAE', 'MOSTA_brain_results', paste0(tp, '_spaVAE_KMeans.png')), plot=p, height=5, width=5)
}
write.csv(spaVAE_pred_df, file.path(result_dir, 'explore_spaVAE', 'MOSTA_brain_results', 'spaVAE_pred.csv'))

# --- compute lisi
domain_metrics_df = read.csv(file.path(dirname(spatial_layout_dir), "final_MOSTA_brain_domain_metrics.csv"), header=T)
domain_metrics_df$X = NULL
domain_metrics_df = domain_metrics_df %>% mutate(method = case_when (
    method == 'GATVAE-SVGP-VAE' ~ 'SpaDOT',
    TRUE ~ method
))
methods = unique(domain_metrics_df$method)
# method_color_palette = c("#F21A00", "#F3DF6C", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#446455")
method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae")
names(method_color_palette) = c('SpaDOT', 'spatialPCA', 'GraphST', 'SEDR', 'spaGCN', 'spaVAE', 'STAligner', 'PRECAST', 'spatiAlign', 'SPIRAL', 'STADIA')

library(lisi)
for (method in methods) {
    if (method == 'SpaDOT') {
        method_pred_file = file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'MOSTA_brain_results', 'spaDOT_pred.csv')
    } else {
        method_pred_file = file.path(result_dir, paste0('explore_', method), 'MOSTA_brain_results', paste0(method, '_pred.csv'))
    }
    method_pred_df = read.csv(method_pred_file, header=T, row.names=1)
    for (tp in unique(method_pred_df$timepoint)) {
        tp_method_pred_df = method_pred_df[method_pred_df$timepoint == tp, ]
        locations = data.frame(x=tp_method_pred_df$pixel_x, y=tp_method_pred_df$pixel_y)
        lisi = compute_lisi(locations, tp_method_pred_df, c('pred_labels'))
        median_lisi = median(lisi[, 1])
        domain_metrics_df[nrow(domain_metrics_df)+1,] = c('median-LISI', tp, median_lisi, method)
    }
}

# --- rename time point
domain_metrics_df = domain_metrics_df %>% mutate(timepoint = case_when(
    timepoint == 't0' ~ 'E12.5',
    timepoint == 't1' ~ 'E13.5',
    timepoint == 't2' ~ 'E14.5',
    timepoint == 't3' ~ 'E15.5',
    timepoint == 't4' ~ 'E16.5',
    TRUE ~ timepoint
))
# domain_metrics_df$method = factor(domain_metrics_df$method, levels=unique(domain_metrics_df$method))
domain_metrics_df$method = factor(domain_metrics_df$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST'))
metric_df = domain_metrics_df[domain_metrics_df$index == 'median-LISI', ]
metric_df$value = as.numeric(metric_df$value)
# compute average
avg_metric_df = metric_df %>% group_by(method) %>% summarize(value=mean(value), .groups="drop") %>% as.data.frame()
avg_metric_df$timepoint = 'avg'
avg_metric_df$index = 'median-LISI'
avg_metric_df = avg_metric_df[, colnames(domain_metrics_df)]
domain_metrics_df = rbind(domain_metrics_df, avg_metric_df)
domain_metrics_df$value = as.numeric(domain_metrics_df$value)
write.csv(domain_metrics_df, file.path(dirname(spatial_layout_dir), "final_MOSTA_brain_domain_metrics_medianLISIadded.csv"), row.names=FALSE)

for (metric in unique(domain_metrics_df$index)) {
    metric_df = domain_metrics_df[domain_metrics_df$index == metric & domain_metrics_df$timepoint != 'avg', ]
    avg_metric_df = domain_metrics_df[domain_metrics_df$index == metric & domain_metrics_df$timepoint == 'avg', ]
    p = ggplot() +
        geom_point(data=avg_metric_df, aes(x=method, y=value, color=method, group=method), size=6) +
        scale_color_manual(values=method_color_palette) + 
        geom_point(data=metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.5) +
        theme_pubr() +
        theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric)
    p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
    legend = get_plot_component(p, "guide-box", return_all=TRUE)[[4]]
    ggsave(file.path(dirname(spatial_layout_dir), 'legend.pdf'), cowplot::plot_grid(legend), width = 8, height = 2)
    p = p + theme(legend.position="none")
    ggsave(file.path(dirname(spatial_layout_dir), paste0(metric, ".pdf")), plot=p, width=4, height=3.5)
}
