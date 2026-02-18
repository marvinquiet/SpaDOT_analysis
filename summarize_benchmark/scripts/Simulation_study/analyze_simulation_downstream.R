# conda activate ~/envs/spatialATAC
Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(dplyr)
library(Rtsne)
library(umap)
source('simulation_utils.R')

result_dir = '/net/zootopia/disk1/wenjinma/minicluster/minicluster_backup_06102025/wenjinma/projects/spatialATAC-explore'
downstream_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/summarize_benchmark/figures/simulation_figures_updated"
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)

# --- summarize running time
simulation_seeds = c(1993, 2024, 1105, 327, 113, 905, 789, 1314, 888, 666)
simulation_combs = c('20_30', '20_50', '30_30', '30_50')
methods = c('spaDOT', 'spatialPCA', 'GraphST', 'SEDR', 'spaGCN', 'spaVAE', 'STAligner', 'PRECAST', 'spatiAlign', 'SPIRAL', 'STADIA', 'STAMP')
simu_df = data.frame(matrix(nrow=0, ncol=4))
colnames(simu_df) = c('seed', 'comb', 'method', 'time')
for (seed in simulation_seeds) {
    for (comb in simulation_combs) {
        for (method in methods) {
            if (method == 'spaDOT') {
                time = scan(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'simulation_results', seed, comb, 'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1', 'spaDOT_running_time.txt'), what=character(), n=1)
            } else {
                time = scan(file.path(result_dir, paste0('explore_', method), 'simulation_results', seed, comb, paste0(method, '_running_time.txt')), what=character(), n=1)
            }
            simu_df[nrow(simu_df)+1, ] = c(seed, comb, method, time)
        }
    }
}
write.csv(simu_df, file.path(downstream_dir, 'runtime.csv'), quote=F)
# --- spaDOT
spaDOT_simu_df = simu_df[simu_df$method == 'spaDOT',]
spaDOT_simu_df %>% group_by(comb) %>% summarize(avg_time=mean(as.numeric(time)), median_time=median(as.numeric(time)))
# --- rename spaDOT to SpaDOT
simu_df$method[simu_df$method == 'spaDOT'] = 'SpaDOT'
method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae", "#ABDDDE")
customized_method_colors = method_color_palette[1:length(unique(simu_df$method))]
names(customized_method_colors) = unique(simu_df$method)

# plot boxplot for each method
simu_df$time = as.numeric(simu_df$time)
# simu_df$method = factor(simu_df$method, levels=unique(simu_df$method))
simu_df$method = factor(simu_df$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST', 'STAMP'))

# log10 transform
# simu_df$time = log10(simu_df$time)

p = ggplot() +
    geom_boxplot(data=simu_df, aes(x=method, y=time, color=method, group=method)) +
    scale_color_manual(values=customized_method_colors) + 
    theme_pubr() +
    theme(axis.text.x = element_blank()) + xlab(NULL) + ylab('Running time(s)')
    # facet_wrap(~comb, ncol=4)
p = p + guides(color = guide_legend(ncol=1, title=NULL), shape=guide_legend(ncol=1, title=NULL))
legend = get_legend(p)
ggsave(file.path(downstream_dir, 'legend.pdf'), plot=legend, width=2, height=6)

p = p + theme(legend.position="none")
ggsave(file.path(downstream_dir, "running_time_nolog.pdf"), plot=p, width=4, height=3)

simu_df$time = log10(simu_df$time)

p = ggplot() +
    geom_boxplot(data=simu_df, aes(x=method, y=time, color=method, group=method)) +
    scale_color_manual(values=customized_method_colors) + 
    theme_pubr() +
    theme(axis.text.x = element_blank()) + xlab(NULL) + ylab('Running time(s) - log10')
    # facet_wrap(~comb, ncol=4)
p = p + guides(color = guide_legend(ncol=1, title=NULL), shape=guide_legend(ncol=1, title=NULL))
legend = get_legend(p)
ggsave(file.path(downstream_dir, 'legend.pdf'), plot=legend, width=2, height=6)

p = p + theme(legend.position="none")
ggsave(file.path(downstream_dir, "running_time_log.pdf"), plot=p, width=4, height=3)


# --- load data for seed 1993
# sim_combinations = c('1993/20_30', '1993/20_50', '1993/30_30', '1993/30_50')
sim_comb = '1993/20_30'
latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'simulation_results', sim_comb, 'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1', 'prediction_adata.h5ad'))
latent_adata$obs$spaDOT_pred_labels = paste(latent_adata$obs$timepoint, latent_adata$obs$pred_labels, sep='_')

latent_adata$obs$renamed_GT_labels = rename_GT_labels(latent_adata)
timepoints = as.vector(unique(latent_adata$obs$timepoint))
timepoints = sort(timepoints)
g_total = NULL
for (tp in timepoints) {
    g = ggplot(subset(latent_adata$obs, timepoint == tp), aes(x=pixel_x, y=pixel_y, color=renamed_GT_labels)) + 
    geom_point() + theme_void() + 
    scale_color_manual(values=GT_custom_palette()) + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (is.null(g_total)) {
        g_total = g
    } else {
        g_total = g_total + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_ground_truth.png')), height=4, width=10)

#custom_palette = c(
#    't0_0'='#EA3423',
#    't0_1'='#76FC4C',
#    't0_2'='#FEFC53',
#    't1_0'='#EA33F4',
#    't1_1'='#FEFC53',
#    't1_2'='#891AF5',
#    't1_3'='#76FC4C',
#    't2_0'='#76FC4C',
#    't2_4'='#F5C242',
#    't2_2'='#891AF5',
#    't2_1'='#716A1E',
#    't2_3'='#EA33F4'
#)

# ---  assign domain colors based on domain overlap
assign_colors = function(GT_domain_label, pred_domain_label) {
    library(clue)
    conf_matrix = table(GT_domain_label, pred_domain_label)
    assignment = solve_LSAP(conf_matrix, maximum = TRUE)
    mapping = data.frame(
      pred_label = colnames(conf_matrix)[assignment],
      true_label = rownames(conf_matrix)
    )
    mapping$true_color_label = sapply(strsplit(mapping$true_label, '_'), '[[', 2)
    mapping$true_color = unname(GT_custom_palette()[mapping$true_color_label])
    custome_palette = mapping$true_color
    names(custome_palette) = mapping$pred_label
    return(custome_palette)
}


g_total = NULL
for (tp in timepoints) {
    g = ggplot(subset(latent_adata$obs, timepoint == tp), aes(x=pixel_x, y=pixel_y, color=spaDOT_pred_labels)) + 
    geom_point() + theme_void() + 
    scale_color_manual(values=assign_colors(paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$renamed_GT_labels), latent_adata$obs$spaDOT_pred_labels)) + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (is.null(g_total)) {
        g_total = g
    } else {
        g_total = g_total + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spaDOT_prediction.png')), height=4, width=10)


spatiAlign_latent = read_h5ad(file.path(result_dir, 'explore_spatiAlign', 'simulation_results', sim_comb, 'prediction_adata.h5ad'))
# -- rename spatiAlign_latent
spatiAlign_latent$obs_names = sub("-.*$", "", spatiAlign_latent$obs_names)
latent_adata$obs$spatiAlign_pred_labels = spatiAlign_latent$obs[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$spatiAlign_pred_labels = paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$spatiAlign_pred_labels)
g_total = NULL
for (tp in timepoints) {
    g = ggplot(subset(latent_adata$obs, timepoint == tp), aes(x=pixel_x, y=pixel_y, color=spatiAlign_pred_labels)) + 
    geom_point() + theme_void() + 
    scale_color_manual(values=assign_colors(paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$renamed_GT_labels), 
                                            latent_adata$obs$spatiAlign_pred_labels)) + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (is.null(g_total)) {
        g_total = g
    } else {
        g_total = g_total + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spatiAlign_prediction.png')), height=4, width=10)

# graphST results
graphST_latent = read_h5ad(file.path(result_dir, 'explore_GraphST', 'simulation_results', sim_comb, 'prediction_adata.h5ad'))
latent_adata$obs$graphST_pred_labels = graphST_latent$obs[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$graphST_pred_labels = paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$graphST_pred_labels)
g_total = NULL
for (tp in timepoints) {
    g = ggplot(subset(latent_adata$obs, timepoint == tp), aes(x=pixel_x, y=pixel_y, color=graphST_pred_labels)) + 
    geom_point() + theme_void() + 
    scale_color_manual(values=assign_colors(paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$renamed_GT_labels),
                                                    latent_adata$obs$graphST_pred_labels)) +
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (is.null(g_total)) {
        g_total = g
    } else {
        g_total = g_total + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_graphST_prediction.png')), height=4, width=10)

spaVAE_latent = read_h5ad(file.path(result_dir, 'explore_spaVAE', 'simulation_results', sim_comb, 'prediction_adata.h5ad'))
latent_adata$obs$spaVAE_pred_labels = spaVAE_latent$obs[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$spaVAE_pred_labels = paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$spaVAE_pred_labels)
g_total = NULL
for (tp in timepoints) {
    g = ggplot(subset(latent_adata$obs, timepoint == tp), aes(x=pixel_x, y=pixel_y, color=spaVAE_pred_labels)) + 
    geom_point() + theme_void() + 
    scale_color_manual(values=assign_colors(paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$renamed_GT_labels),
                                                   latent_adata$obs$spaVAE_pred_labels)) +
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (is.null(g_total)) {
        g_total = g
    } else {
        g_total = g_total + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spaVAE_prediction.png')), height=4, width=10)

#spatialPCA_latent = read_h5ad(file.path(result_dir, 'explore_spatialPCA', 'simulation_results', sim_comb, 'prediction_adata.h5ad'))
#latent_adata$obs$spatialPCA_pred_labels = spatialPCA_latent$obs[latent_adata$obs_names, 'pred_labels']
#g_total = NULL
#for (tp in timepoints) {
#    g = ggplot(subset(latent_adata$obs, timepoint == tp), aes(x=pixel_x, y=pixel_y, color=spatialPCA_pred_labels)) + 
#    geom_point() + theme_void() + 
#    theme(legend.position = "bottom") + 
#    guides(color = guide_legend(ncol = 3, title=NULL))
#    if (is.null(g_total)) {
#        g_total = g
#    } else {
#        g_total = g_total + g
#    }
#}

# create UMAP plots
umap_result = umap(latent_adata$X)
latent_adata$obs$spaDOT_umap1 = umap_result$layout[,1]
latent_adata$obs$spaDOT_umap2 = umap_result$layout[,2]
g_list = NULL
# for (variable in c('timepoint', 'celltype', 'renamed_GT_labels', 'spaDOT_pred_labels')) {
for (variable in c('timepoint', 'celltype', 'renamed_GT_labels')) {
    g = ggplot(latent_adata$obs, aes(x=spaDOT_umap1, y=spaDOT_umap2, color=!!sym(variable))) + geom_point(size=0.1, shape=19) + theme_pubr() + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL)) + 
    labs(x = 'SpaDOT_UMAP1', y = 'SpaDOT_UMAP2')
    if (variable == 'renamed_GT_labels') {
        g = g + scale_color_manual(values=GT_custom_palette())
    }
    if (variable == 'spaDOT_pred_labels') {
        g = g + scale_color_manual(values=custom_palette)
    }
    if (is.null(g_list)) {
        g_list = g
    } else {
        g_list = g_list + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spaDOT_umap_all.png')), g_list, height=4, width=8)

umap_result = umap(spatiAlign_latent[latent_adata$obs_names]$obsm$correct)
latent_adata$obs$spatiAlign_umap1 = umap_result$layout[,1]
latent_adata$obs$spatiAlign_umap2 = umap_result$layout[,2]
latent_adata$obs$spatiAlign_final_pred_labels = paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$spatiAlign_pred_labels)
g_list = NULL
for (variable in c('timepoint', 'celltype', 'renamed_GT_labels')) {
    g = ggplot(latent_adata$obs, aes(x=spatiAlign_umap1, y=spatiAlign_umap2, color=!!sym(variable))) + geom_point(size=0.1, shape=19, alpha=0.3) + theme_pubr() + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (variable == 'renamed_GT_labels') {
        g = g + scale_color_manual(values=GT_custom_palette())
    }
    if (is.null(g_list)) {
        g_list = g
    } else {
        g_list = g_list + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spatiAlign_umap_all.png')), g_list, height=4, width=8)

umap_result = umap(graphST_latent[latent_adata$obs_names]$obsm$emb)
latent_adata$obs$graphST_umap1 = umap_result$layout[,1]
latent_adata$obs$graphST_umap2 = umap_result$layout[,2]
latent_adata$obs$graphST_final_pred_labels = paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$graphST_pred_labels)
g_list = NULL
for (variable in c('timepoint', 'celltype', 'renamed_GT_labels')) {
    g = ggplot(latent_adata$obs, aes(x=graphST_umap1, y=graphST_umap2, color=!!sym(variable))) + geom_point(size=0.1, shape=19) + theme_pubr() + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (variable == 'renamed_GT_labels') {
        g = g + scale_color_manual(values=GT_custom_palette())
    }
    if (is.null(g_list)) {
        g_list = g
    } else {
        g_list = g_list + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_graphST_umap_all.png')), g_list, height=4, width=8)

umap_result = umap(spaVAE_latent[latent_adata$obs_names]$X)
latent_adata$obs$spaVAE_umap1 = umap_result$layout[,1]
latent_adata$obs$spaVAE_umap2 = umap_result$layout[,2]
latent_adata$obs$spaVAE_final_pred_labels = paste0(latent_adata$obs$timepoint, '_', latent_adata$obs$spaVAE_pred_labels)
g_list = NULL
for (variable in c('timepoint', 'celltype', 'renamed_GT_labels')) {
    g = ggplot(latent_adata$obs, aes(x=spaVAE_umap1, y=spaVAE_umap2, color=!!sym(variable))) + geom_point(size=0.1, shape=19) + theme_pubr() + 
    theme(legend.position = "bottom") + 
    guides(color = guide_legend(ncol = 3, title=NULL))
    if (variable == 'renamed_GT_labels') {
        g = g + scale_color_manual(values=GT_custom_palette())
    }
    if (is.null(g_list)) {
        g_list = g
    } else {
        g_list = g_list + g
    }
}
ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spaVAE_umap_all.png')), g_list, height=4, width=8)


# tsne_result = Rtsne(latent_adata$X)
# latent_adata$obs$spaDOT_tsne1 = tsne_result$Y[,1]
# latent_adata$obs$spaDOT_tsne2 = tsne_result$Y[,2]
# ggplot(latent_adata$obs, aes(x=spaDOT_tsne1, y=spaDOT_tsne2, color=spaDOT_pred_labels)) + geom_point() + theme_pubr()
# ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spaDOT_tsne.png')), height=4, width=4)

# --- analyze results
# method_color_palette =c("#F21A00", "#F3DF6C", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#446455")
# method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae")
simulation_seeds = c(1993, 2024, 1105, 327, 113, 905, 789, 1314, 888, 666)
simulation_combs = c('20_30', '20_50', '30_30', '30_50')
simu_df = data.frame()
for (seed in simulation_seeds) {
    for (comb in simulation_combs) {
        simu_comb = paste0(seed, '-', comb)
        metric_df = read.csv(file.path(downstream_dir, simu_comb, "final_simulation_domain_metrics.csv"), header=T)
        metric_df$simu_seed = seed
        metric_df$simu_comb = comb
        simu_df = rbind(simu_df, metric_df)
    }
}
simu_df$X = NULL
simu_df$method[simu_df$method == 'GATVAE-SVGP-VAE'] = 'SpaDOT'
#customized_method_colors = method_color_palette[1:length(unique(simu_df$method))]
#names(customized_method_colors) = unique(simu_df$method)
# simu_df$method = factor(simu_df$method, levels=unique(simu_df$method))
simu_df$method = factor(simu_df$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST', 'STAMP'))


# --- calculate LISI for all methods
methods = unique(simu_df$method)
library(lisi)
timepoints = c('t0', 't1', 't2')
for (seed in simulation_seeds) {
    for (comb in simulation_combs) {
        for (method in methods) {
            if (method == 'SpaDOT') {
                pred_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'simulation_results', paste0(seed, '/', comb), 'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1', 'prediction_adata.h5ad'))
            } else {
                pred_adata = read_h5ad(file.path(result_dir, paste0('explore_', method), 'simulation_results', paste0(seed, '/', comb), 'prediction_adata.h5ad'))
            }
            for (tp in unique(pred_adata$obs$timepoint)) {
                tp_pred_adata = pred_adata[pred_adata$obs$timepoint == tp, ]
                locations = data.frame(x=tp_pred_adata$obs$pixel_x, y=tp_pred_adata$obs$pixel_y)
                lisi = compute_lisi(locations, tp_pred_adata$obs, c('pred_labels'))
                median_lisi = median(lisi[, 1])
                simu_df[nrow(simu_df)+1,] = c('median-LISI', tp, median_lisi, method, seed, comb)
            }
        }
    }
}

simu_df$value = as.numeric(simu_df$value)
simu_df$simu_comb = factor(simu_df$simu_comb, levels=c('30_50', '30_30', '20_50', '20_30'))
metrics = unique(simu_df$index)
for (metric in metrics) {
    metric_df = simu_df[simu_df$index == metric,]
    for (tp in unique(metric_df$timepoint)) {
        tp_metric_df = metric_df[metric_df$timepoint == tp,]
        # calculate average
        avg_metric_df = tp_metric_df %>% group_by(method, simu_comb) %>% summarize(avg_value=mean(value), median_value=median(value))

        p = ggplot() +
            geom_boxplot(data=tp_metric_df, aes(x=method, y=value, color=method, group=method)) +
            scale_color_manual(values=customized_method_colors) + 
            theme_pubr() +
            theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric) + 
            facet_wrap(~simu_comb, ncol=4)
        p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
        legend = get_legend(p)
        p = p + theme(legend.position="none")
        ggsave(file.path(downstream_dir, paste0(tp, '_', metric, ".pdf")), plot=p, width=12, height=3)
    }

    timepoint_metric_df = metric_df[metric_df$timepoint != 'avg',]
    if (metric == 'median-LISI') {
        # for median-LISI, we need to calculate the average across timepoints
        avg_metric_df = timepoint_metric_df %>% group_by(method, simu_seed, simu_comb) %>% summarize(value=mean(value))
        metric_df %>% group_by(method, simu_comb) %>% summarize(avg=mean(value), median=median(value), .groups="drop") %>% as.data.frame()
    } else {
        avg_metric_df = metric_df[metric_df$timepoint == 'avg',]
    }
    p = ggplot() +
        geom_boxplot(data=avg_metric_df, aes(x=method, y=value, color=method, group=method)) +
        scale_color_manual(values=customized_method_colors) + 
        # geom_point(data=timepoint_metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.5) +
        theme_pubr() +
        theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric) + 
        facet_wrap(~simu_comb, ncol=4)
    p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
    legend = get_legend(p)
    ggsave(file.path(downstream_dir, 'legend.pdf'), plot=legend, width=6, height=2)
    p = p + theme(legend.position="none")
    ggsave(file.path(downstream_dir, paste0(metric, ".pdf")), plot=p, width=12, height=3)
}
write.csv(simu_df, file.path(downstream_dir, "final_simulation_domain_metrics.csv"), row.names=FALSE)




# --- plot for seed 1993
for (sim_comb in c('1993/20_30', '1993/30_30', '1993/20_50', '1993/30_50')) {
    latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'simulation_results', sim_comb, 'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1', 'prediction_adata.h5ad'))
    latent_adata$obs$spaDOT_pred_labels = paste(latent_adata$obs$timepoint, latent_adata$obs$pred_labels, sep='_')
    # --- plot prediction
    g_total = NULL
    for (tp in timepoints) {
        g = ggplot(subset(latent_adata$obs, timepoint == tp), aes(x=pixel_x, y=pixel_y, color=spaDOT_pred_labels)) + 
        geom_point() + theme_void() + 
        scale_color_manual(values=custome_palette_list[[sim_comb]]) + 
        theme(legend.position = "bottom") + 
        guides(color = guide_legend(ncol = 3, title=NULL))
        if (is.null(g_total)) {
            g_total = g
        } else {
            g_total = g_total + g
        }
    }
    ggsave(file.path(downstream_dir, paste0(gsub('/', '-', sim_comb), '_spaDOT_prediction.png')), height=4, width=10)

    # plot OT heatmap
    timepoints = as.vector(unique(latent_adata$obs$timepoint))
    timepoints = sort(timepoints)
    plot_OT_heatmap(timepoints, file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'simulation_results', sim_comb, 'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1'), custome_palette_list[[sim_comb]])

    # --- perform Harmony
    library(harmony)
    meta_data = latent_adata$obs
    harmony_latent = RunHarmony(latent_adata$X, meta_data, 'timepoint')
    write.csv(harmony_latent, file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'simulation_results', sim_comb, 'lambda1_0.1_lambda2_0_beta1_1_beta2_0.0001_omiga1_0.1_omiga2_0.1_omiga3_1', 'harmony_latent.csv'))
}
