library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(optparse)

ggplot2::theme_set(theme_pubr())

option_list = list(
    make_option(c("--dataset"), type='character', help="Dataset.")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)
dataset = opt$dataset
figure_dir = paste0(dataset, "_figures")

method_color_palette = c("#F3DF6C", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#000000", "#446455")
fill_color_palette = c("gray90", "gray80", "gray70", "gray60", "gray50", "gray40", "gray30", "gray20", "gray10")
plot_domain_metrics = function(domain_metrics_df, figure_dir, dataset) {
    # set colors
    customized_method_colors = method_color_palette[1:length(unique(domain_metrics_df$method))]
    names(customized_method_colors) = unique(domain_metrics_df$method)
    domain_metrics_df$X = NULL
    domain_metrics_df$method = factor(domain_metrics_df$method, levels=unique(domain_metrics_df$method))

    metrics = unique(domain_metrics_df$index)
    for (metric in metrics) {
        metric_df = domain_metrics_df[domain_metrics_df$index == metric,]
        # for (tp in unique(metric_df$timepoint)) {
        #     tp_meric_df = metric_df[metric_df$timepoint == tp,]
        #     ggplot(tp_meric_df, aes(x=method, y=value, color=method, group=method)) + 
        #     geom_point(size=6) +
        #     scale_color_manual(values=customized_method_colors) + 
        #     theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        #     xlab("Methods") + ylab(metric)
        #     ggsave(file.path(figure_dir, paste0(metric, '_', tp, ".png")))
        # }
        avg_metric_df = metric_df[metric_df$timepoint == 'avg',]
        timepoint_metric_df = metric_df[metric_df$timepoint != 'avg',]

        p = ggplot() +
            geom_point(data=avg_metric_df, aes(x=method, y=value, color=method, group=method), size=6) +
            scale_color_manual(values=customized_method_colors) + 
            geom_point(data=timepoint_metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.5) +
            theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric)
        if (dataset == 'ChickenHeart') {
            p = p + guides(color = guide_legend(nrow = 4, title=NULL), shape=guide_legend(nrow=4, title=NULL))

        } else {
            p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
        }
        legend = get_legend(p)
        ggsave(file.path(figure_dir, 'legend.png'), plot=legend, width=6, height=2)
        p = p + theme(legend.position="none")
        ggsave(file.path(figure_dir, paste0(metric, ".png")), plot=p, width=4, height=3.5)
    }
}

if (dataset == 'ATAC-Visium' || dataset == 'ATAC-Visium_noremoval' || dataset == 'Axolotl' || dataset == 'ChickenHeart') {
    domain_metrics_df = read.csv(file.path(figure_dir, paste0("final_", dataset, "_domain_metrics.csv")), header=T)
    # rename one of the methods name
    domain_metrics_df$method[domain_metrics_df$method == 'GATVAE-SVGP-VAE'] = 'spaDOT'
    plot_domain_metrics(domain_metrics_df, figure_dir, dataset)

} else if (dataset == 'simulation') {
    simulation_combinations = c('666-10_30', '666-20_30', '666-30_30',
                                '666-10_50', '666-20_50', '666-30_50')
    for (simulation_comb in simulation_combinations) {
        simulation_figure_dir = file.path(figure_dir, simulation_comb)
        domain_metrics_df = read.csv(file.path(simulation_figure_dir, "final_simulation_domain_metrics.csv"), header=T)
            # rename one of the methods name
        domain_metrics_df$method[domain_metrics_df$method == 'GATVAE-SVGP-VAE'] = 'spaDOT'
        plot_domain_metrics(domain_metrics_df, simulation_figure_dir, dataset)
    }
} else {
    print("Invalid dataset.")
}


# Usage
# Rscript plot_domain_metrics.R --dataset ATAC-Visium
# Rscript plot_domain_metrics.R --dataset Axolotl
# Rscript plot_domain_metrics.R --dataset simulation
