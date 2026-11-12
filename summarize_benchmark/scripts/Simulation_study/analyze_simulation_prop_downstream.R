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
library(glue)
library(optparse)
source('simulation_utils.R')

option_list = list(
    make_option(c("--DEG"), type="integer", help="number of DE genes for each cell type"),
    make_option(c("--mixedNPCs"), action="store_true", default=FALSE, help="Whether the NPCs is a mixture")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)

result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis"
if (opt$mixedNPCs) {
    downstream_dir = file.path(result_dir, glue("summarize_benchmark/figures/simulation_prop_figures_mixedNPCs_DEG{opt$DEG}/simulation_downstream"))
} else {
    downstream_dir = file.path(result_dir, glue("summarize_benchmark/figures/simulation_prop_figures_DEG{opt$DEG}/simulation_downstream"))
}
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)

simulation_seeds = c(1993, 2024, 1105)
simulation_combs = c('prop_0.9', 'prop_0.8', 'prop_0.7', 'prop_0.6', 'prop_0.5')
methods = c('SpaDOT', 'spatialPCA', 'GraphST', 'SEDR', 'spaGCN', 'spaVAE', 'STAligner', 'PRECAST', 'spatiAlign', 'SPIRAL', 'STADIA', 'STAMP')
simu_df = data.frame(matrix(nrow=0, ncol=3))
colnames(simu_df) = c('seed', 'comb', 'method')

method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae", "#ABDDDE")
customized_method_colors = method_color_palette[1:length(methods)]
names(customized_method_colors) = methods

# plot boxplot for each method
simu_df$method = factor(simu_df$method, levels=methods)
# --- analyze results
simu_df = data.frame()
for (seed in simulation_seeds) {
    for (comb in simulation_combs) {
        metric_df = read.csv(file.path(dirname(downstream_dir), seed, comb, "final_simulation_domain_metrics.csv"), header=T)
        metric_df$simu_seed = seed
        metric_df$comb = comb
        simu_df = rbind(simu_df, metric_df)
    }
}
simu_df$X = NULL
simu_df$method = factor(simu_df$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST', 'STAMP'))

# --- calculate LISI for all methods
methods = unique(simu_df$method)
library(lisi)
timepoints = c('t0', 't1', 't2')
for (seed in simulation_seeds) {
    for (comb in simulation_combs) {
        for (method in methods) {
            if (opt$mixedNPCs) {
                pred_adata = read_h5ad(file.path(result_dir, "benchmark_study", paste0('explore_', method), paste0('simulation_prop_results_mixedNPCs_DEG', opt$DEG), paste0(seed, '/', comb), 'prediction_adata.h5ad'))
            } else {
                pred_adata = read_h5ad(file.path(result_dir, "benchmark_study", paste0('explore_', method), paste0('simulation_prop_results_DEG', opt$DEG), paste0(seed, '/', comb), 'prediction_adata.h5ad'))
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
write.csv(simu_df, file.path(downstream_dir, "final_simulation_prop_domain_metrics.csv"), row.names=FALSE)

simu_df$value = as.numeric(simu_df$value)
simu_df$comb = factor(simu_df$comb, levels=simulation_combs)
metrics = unique(simu_df$index)
for (metric in metrics) {
    metric_df = simu_df[simu_df$index == metric,]
    for (tp in unique(metric_df$timepoint)) {
        tp_metric_df = metric_df[metric_df$timepoint == tp,]
        # calculate average
        avg_metric_df = tp_metric_df %>% group_by(method, comb) %>% summarize(avg_value=mean(value), median_value=median(value))

        p = ggplot() +
            geom_boxplot(data=tp_metric_df, aes(x=method, y=value, color=method, group=method)) +
            scale_color_manual(values=customized_method_colors) + 
            theme_pubr() +
            theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric) + 
            facet_wrap(~comb, nrow=1, ncol=length(simulation_combs))
        p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
        legend = get_legend(p)
        p = p + theme(legend.position="none")
        ggsave(file.path(downstream_dir, paste0(tp, '_', metric, ".pdf")), plot=p, width=3*length(simulation_combs), height=3)
    }

    timepoint_metric_df = metric_df[metric_df$timepoint != 'avg',]
    if (metric == 'median-LISI') {
        # for median-LISI, we need to calculate the average across timepoints
        avg_metric_df = timepoint_metric_df %>% group_by(method, simu_seed, comb) %>% summarize(value=mean(value))
        metric_df %>% group_by(method, comb) %>% summarize(avg=mean(value), median=median(value), .groups="drop") %>% as.data.frame()
    } else {
        avg_metric_df = metric_df[metric_df$timepoint == 'avg',]
    }
    p = ggplot() +
        geom_boxplot(data=avg_metric_df, aes(x=method, y=value, color=method, group=method)) +
        scale_color_manual(values=customized_method_colors) + 
        # geom_point(data=timepoint_metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.5) +
        theme_pubr() +
        theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric) + 
        facet_wrap(~comb, nrow=1, ncol=length(simulation_combs))
    p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
    legend = get_legend(p)
    ggsave(file.path(downstream_dir, 'legend.pdf'), plot=legend, width=6, height=2)
    p = p + theme(legend.position="none")
    ggsave(file.path(downstream_dir, paste0(metric, ".pdf")), plot=p, width=3*length(simulation_combs), height=3)
}
