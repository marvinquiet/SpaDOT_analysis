Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggpubr)
library(rhdf5)
library(dplyr)
library(stringr)
library(tidyr)
source('Axolotl_utils.R')

base_dir = "/net/mulan/home/wenjinma"
result_dir = file.path(base_dir, "projects/SpaDOT_analysis")
spatial_layout_dir = file.path(result_dir, "summarize_benchmark/figures/Axolotl_figures/clustering")
dir.create(spatial_layout_dir, showWarnings = FALSE, recursive = TRUE)


# -- here, pred_labels is the results of kmeans
plot_spatial = function(latent_adata, method) {
    latent_adata$obs = latent_adata$obs %>%
        mutate(timepoint_original = case_when(
            timepoint == 't0' ~ "St.44",
            timepoint == 't1' ~ "St.54",
            timepoint == 't2' ~ "St.57"
        )
    )

    clustering_methods = c('kmeans', 'mclust', 'hc', 'leiden', 'louvain')
    for (clustering in clustering_methods) {
        latent_adata$obs[[clustering]] = as.factor(latent_adata$obs[[clustering]])
        plot_pred_label(latent_adata$obs, color_col=clustering, 
            color_legend_title=paste(method, clustering, 'Predicted Labels', sep=' '),
            file_name=file.path(spatial_layout_dir, paste(method, clustering, 'prediction.png', sep='_')))
    }
}

# SpaDOT
latent_adata = read_h5ad(file.path(result_dir, 'explore_SpaDOT', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(latent_adata, method='SpaDOT')

# --- load all methods' results
# spatialPCA
spatialPCA_latent = read_h5ad(file.path(result_dir, 'explore_spatialPCA', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(spatialPCA_latent, method='SpatialPCA')

# GraphST
GraphST_latent = read_h5ad(file.path(result_dir, 'explore_GraphST', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(GraphST_latent, method='GraphST')

# SEDR
SEDR_latent = read_h5ad(file.path(result_dir, 'explore_SEDR', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(GraphST_latent, method='SEDR')

# spaGCN
spaGCN_latent = read_h5ad(file.path(result_dir, 'explore_spaGCN', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(spaGCN_latent, method='spaGCN')

# spaVAE
spaVAE_latent = read_h5ad(file.path(result_dir, 'explore_spaVAE', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(spaVAE_latent, method='spaVAE')

# STAligner
STAligner_latent = read_h5ad(file.path(result_dir, 'explore_STAligner', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(STAligner_latent, method='STAligner')

# PRECAST
PRECAST_latent = read_h5ad(file.path(result_dir, 'explore_PRECAST', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(PRECAST_latent, method='PRECAST')

# spatiAlign
spatiAlign_latent = read_h5ad(file.path(result_dir, 'explore_spatiAlign', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(spatiAlign_latent, method='spatiAlign')

# SPIRAL
SPIRAL_latent = read_h5ad(file.path(result_dir, 'explore_SPIRAL', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(SPIRAL_latent, method='SPIRAL')

# STADIA
STADIA_latent = read_h5ad(file.path(result_dir, 'explore_STADIA', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(STADIA_latent, method='STADIA')

# STAMP
STAMP_latent = read_h5ad(file.path(result_dir, 'explore_STAMP', 'Axolotl_results', 'prediction_adata.h5ad'))
plot_spatial(STAMP_latent, method='STAMP')

methods = c("SpaDOT", "spatialPCA", "GraphST", "SEDR", "spaGCN", "spaVAE", "STAligner", "PRECAST", "spatiAlign", "SPIRAL", "STADIA", 'STAMP')
clusterings = c('kmeans', 'hc', 'louvain', 'leiden', 'mclust')
clustering_df = NULL
for (method in methods) {
    for (clustering in clusterings) {
        method_clustering_df = read.csv(file.path(result_dir, paste0('explore_', method), 'Axolotl_results', paste(method, clustering, 'domain_metrics.csv', sep='_')),
                header=TRUE, row.names=1)
        method_clustering_df_long = method_clustering_df %>%
                tibble::rownames_to_column("metric") %>%   # convert rownames to a column
                pivot_longer(cols = c(t0, t1, t2, avg),
                            names_to = "timepoint",
                            values_to = "value") %>%
                as.data.frame()
        method_clustering_df_long$method = method
        method_clustering_df_long$clustering = clustering
        if (is.null(clustering_df)){
            clustering_df = method_clustering_df_long
        } else {
            clustering_df = rbind(clustering_df, method_clustering_df_long)
        }
    }
}
write.csv(clustering_df, file.path(spatial_layout_dir, 'all_clustering_results.csv'))

avg_clustering_df = clustering_df[clustering_df$timepoint == 'avg', ]
methods = unique(avg_clustering_df$method)
method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae", "#ABDDDE")
customized_method_colors = method_color_palette[1:length(methods)]
names(customized_method_colors) = methods

for (metric in c('ARI', 'AMI', 'H')) {
    avg_clustering_metric_df = avg_clustering_df[avg_clustering_df$metric == metric, ]
    SpaDOT_res = avg_clustering_metric_df$value[
        avg_clustering_metric_df$method == "SpaDOT" &
        avg_clustering_metric_df$clustering == "kmeans"
    ]
    cat('---', metric, '---\n')
    if (metric in c('ARI', 'AMI', 'H')) {
        print(avg_clustering_metric_df[avg_clustering_metric_df$value >= SpaDOT_res, ])
    }


    # plot per metric
    clustering_metric_df = clustering_df[clustering_df$timepoint != 'avg' & clustering_df$metric == metric, ]
    clustering_metric_df$method = factor(clustering_metric_df$method, levels=c('SpaDOT', 'spatiAlign', 'SpatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST', 'STAMP'))
    clustering_metric_df$clustering = factor(clustering_metric_df$clustering, levels=c('kmeans', 'leiden', 'louvain', 'mclust', 'hc'))
    p = ggplot() +
        geom_point(data=avg_clustering_metric_df, aes(x=method, y=value, color=method, group=method), size=6) +
        scale_color_manual(values=customized_method_colors) + 
        geom_point(data=clustering_metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.5) +
        facet_wrap(~clustering, scale='free_x', nrow=1) + 
        theme_pubr() +
        theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric)
    p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
    legend = get_legend(p)
    ggsave(file.path(spatial_layout_dir, 'legend.pdf'), plot=legend, width=6, height=2)
    p = p + theme(legend.position="none")
    ggsave(file.path(spatial_layout_dir, paste0(metric, "_clustering.pdf")), plot=p, width=20, height=4)
}
