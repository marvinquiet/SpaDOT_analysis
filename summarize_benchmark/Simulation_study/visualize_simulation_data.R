Sys.setenv(RETICULATE_PYTHON = "/media/sda/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)

library(anndata)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(dplyr)

# load simulated data
simulated_data_dir = '/media/sda/wenjinma/data/spatial_simulation_brain/666/20_50'
adata = read_h5ad(file.path(simulated_data_dir, 'adata.h5ad'))
adata$obs$label = factor(adata$obs$label)

spaDOT_pred_dir = "/media/sda/wenjinma/projects/spatialATAC-explore/explore_GATVAE-SVGP-VAE/simulation_results/666/20_50"
spaDOT_latent_adata = read_h5ad(file.path(spaDOT_pred_dir, 'prediction_adata.h5ad'))
common_cells = intersect(adata$obs_names, spaDOT_latent_adata$obs_names)
adata = adata[common_cells,]
spaDOT_latent_adata = spaDOT_latent_adata[common_cells,]
adata$obs$spaDOT_pred = spaDOT_latent_adata$obs[adata$obs_names, 'pred_labels']

# another method spatiAlign
spatiAlign_pred_dir = "/media/sda/wenjinma/projects/spatialATAC-explore/explore_spatiAlign/simulation_results/666/20_50"
spatiAlign_latent_adata = read_h5ad(file.path(spatiAlign_pred_dir, 'prediction_adata.h5ad'))
adata$obs$spatiAlign_pred = spatiAlign_latent_adata$obs[adata$obs_names, 'pred_labels']

spaVAE_pred_dir = "/media/sda/wenjinma/projects/spatialATAC-explore/explore_spaVAE/simulation_results/666/20_50"
spaVAE_latent_adata = read_h5ad(file.path(spaVAE_pred_dir, 'prediction_adata.h5ad'))
adata$obs$spaVAE_pred = spaVAE_latent_adata$obs[adata$obs_names, 'pred_labels']

spaGCN_pred_dir = "/media/sda/wenjinma/projects/spatialATAC-explore/explore_spaGCN/simulation_results/666/20_50"
spaGCN_latent_adata = read_h5ad(file.path(spaGCN_pred_dir, 'prediction_adata.h5ad'))
adata$obs$spaGCN_pred = spaGCN_latent_adata$obs[adata$obs_names, 'pred_labels']

# mutate labels
adata$obs = adata$obs %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "t1",
        timepoint == 't1' ~ "t2",
        timepoint == 't2' ~ "t3"
    ))

adata$obs = adata$obs %>%
    mutate(label_new = case_when(
        timepoint == 't0' & label == 1 ~ "1",
        timepoint == 't0' & label == 2 ~ "2",
        timepoint == 't0' & label == 3 ~ "3",
        timepoint == 't1' & label == 3 ~ "1-1",
        timepoint == 't1' & label == 2 ~ "1-2",
        timepoint == 't1' & label == 1 ~ "2",
        timepoint == 't1' & label == 4 ~ "3",
        timepoint == 't2' & label == 1 ~ "2",
        timepoint == 't2' & label == 2 ~ "3-1",
        timepoint == 't2' & label == 3 ~ "1-2",
        timepoint == 't2' & label == 4 ~ "1-1",
        timepoint == 't2' & label == 5 ~ "3-2"
    ))
adata$obs$label_new = factor(adata$obs$label_new)

color_palettes = c("1"="red", "1-1"="purple", "1-2"="pink", "2"="green", "3"="yellow", "3-1"="brown", "3-2"="gold")
# draw spatial layout of region labels
ggplot(adata$obs, aes(x=x_coor, y=y_coor, color=label_new)) + geom_point(size=0.5) +
    facet_wrap(~timepoint_original, nrow=1, scale='free') +
    theme_pubr() + 
    theme(axis.text =element_blank(), axis.title=element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12)) + 
    scale_color_manual(name='Region Annotation', values=color_palettes)
ggsave(file.path(simulated_data_dir, 'simulated_region_annotation.png'), height=4, width=10)

color_palettes = c("GABA"="red", "Neurons"="green", "Astrocytes"="yellow", "OPC"="purple")
ggplot(adata$obs, aes(x=x_coor, y=y_coor, color=celltype)) + geom_point(size=0.5, alpha=0.5) +
    facet_wrap(~timepoint_original, nrow=1, scale='free') +
    theme_pubr() + 
    theme(axis.text =element_blank(), axis.title=element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12)) + 
    scale_color_manual(name='Cell Type Annotation', values=color_palettes)
ggsave(file.path(simulated_data_dir, 'simulated_celltype_annotation.png'), height=4, width=10)

ggplot(adata$obs, aes(x=x_coor, y=y_coor, color=spaDOT_pred)) + geom_point(size=0.5, alpha=0.5) +
    facet_wrap(~timepoint_original, nrow=1, scale='free') +
    theme_pubr() + 
    theme(axis.text =element_blank(), axis.title=element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12)) + 
    scale_color_discrete(name='spaDOT Prediction')
ggsave(file.path(spaDOT_pred_dir, 'spaDOT_prediction.png'), height=4, width=10)

ggplot(adata$obs, aes(x=x_coor, y=y_coor, color=spatiAlign_pred)) + geom_point(size=0.5, alpha=0.5) +
    facet_wrap(~timepoint_original, nrow=1, scale='free') +
    theme_pubr() + 
    theme(axis.text =element_blank(), axis.title=element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12)) + 
    scale_color_discrete(name='spatiAlign Prediction')
ggsave(file.path(spatiAlign_pred_dir, 'spatiAlign_prediction.png'), height=4, width=10)

ggplot(adata$obs, aes(x=x_coor, y=y_coor, color=spaVAE_pred)) + geom_point(size=0.5, alpha=0.5) +
    facet_wrap(~timepoint_original, nrow=1, scale='free') +
    theme_pubr() + 
    theme(axis.text =element_blank(), axis.title=element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12)) + 
    scale_color_discrete(name='spaVAE Prediction')
ggsave(file.path(spaVAE_pred_dir, 'spaVAE_prediction.png'), height=4, width=10)

ggplot(adata$obs, aes(x=x_coor, y=y_coor, color=spaGCN_pred)) + geom_point(size=0.5, alpha=0.5) +
    facet_wrap(~timepoint_original, nrow=1, scale='free') +
    theme_pubr() + 
    theme(axis.text =element_blank(), axis.title=element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12)) + 
    scale_color_discrete(name='spaGCN Prediction')
ggsave(file.path(spaGCN_pred_dir, 'spaGCN_prediction.png'), height=4, width=10)