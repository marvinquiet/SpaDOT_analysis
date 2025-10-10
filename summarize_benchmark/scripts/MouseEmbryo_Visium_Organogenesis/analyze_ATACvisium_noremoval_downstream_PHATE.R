Sys.setenv(RETICULATE_PYTHON = "/media/sda/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)

library(anndata)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(dplyr)
library(reshape2)
library(phateR)


base_dir = "/media/sda/wenjinma"
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

data_phate = phate(adata$X)
data_phate_labels = paste0(adata$obs$timepoint_original, "_", adata$obs$pred_labels)
data_phate$labels = data_phate_labels
ggplot(data_phate, aes(x=PHATE1, y=PHATE2, color=data_phate_labels)) + geom_point(size=0.5) + theme_pubr() + theme(legend.position="right")
ggsave(file.path(downstream_dir, 'phate_mesenchymal_clusters.png'), width=5, height=5)


data_phate = phate(adata$X, gamma=0, t=10)
ggplot(data_phate, aes(x=PHATE1, y=PHATE2, color=data_phate_labels)) + geom_point(size=0.5) + theme_pubr() + theme(legend.position="right")
ggsave(file.path(downstream_dir, 'phate_mesenchymal_clusters_gamma0_t10.png'), width=5, height=5)


data_pca = prcomp(adata$X)$x
ggplot(data_pca, aes(x=PC1, y=PC2, color=data_phate_labels)) + geom_point(size=0.5) + theme_pubr() + theme(legend.position="right")
ggsave(file.path(downstream_dir, 'PCA_mesenchymal_clusters.png'), width=5, height=5)


png(file.path(downstream_dir, 'PCA_mesenchymal_clusters.png'))


plot(prcomp(adata$X)$x, col=data_phate_labels)
dev.off()