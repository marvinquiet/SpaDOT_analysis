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
source('ChickenHeart_utils.R')

base_dir = "/media/sda/wenjinma"
result_dir = file.path(base_dir, "projects/spatialATAC-explore")
spatial_layout_dir = file.path(result_dir, "summarize_benchmark/ChickenHeart_figures/spatial_layout")
dir.create(spatial_layout_dir, showWarnings = FALSE, recursive = TRUE)
downstream_dir = file.path(result_dir, "summarize_benchmark/ChickenHeart_figures/ChickenHeart_downstream")
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)
latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ChickenHeart_results', 'prediction_adata.h5ad'))


annotations = unique(latent_adata$obs$true_labels)
customized_color_palette = c('#000000',  # Epicardium- like  
'#00b4d8', # Ventricle
'#ffb3c6', # Valves
'#fca311', # Outflow tract      
'#a3b18a', # Atria
'#000000', # Epicardium
'#344e41', # Endothelium
'#caf0f8', # Trabecular LV and \nendocardium          
'#0077b6', # Right ventricle
'#90e0ef', # Compact LV and \ninter-ventricular septum
'#48cae4' # Endocardium
)

# customized_color_palette = color_palette[1:length(annotations)]
names(customized_color_palette) = annotations
# customized_color_palette['Epicardium- like'] = '#F98400'
# customized_color_palette['Endocardium'] = "#5BBCD6"
# --- edit the coordinates
latent_adata$obs$pixel_x_bak = latent_adata$obs$pixel_x
latent_adata$obs$pixel_y_bak = latent_adata$obs$pixel_y
latent_adata$obs = latent_adata$obs %>%
    mutate(pixel_x = case_when(
        timepoint == 't0' ~ pixel_y_bak,
        timepoint == 't1' ~ pixel_y_bak,
        timepoint == 't2' ~ pixel_y_bak,
        timepoint == 't3' ~ pixel_y_bak
    )) %>%
    mutate(pixel_y = pixel_x_bak) 

# --- draw ground truth
latent_adata$obs = latent_adata$obs %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "Day 4",
        timepoint == 't1' ~ "Day 7",
        timepoint == 't2' ~ "Day 10",
        timepoint == 't3' ~ "Day 14"
    ))
latent_adata$obs$timepoint_original = factor(latent_adata$obs$timepoint_original, levels=c('Day 4', 'Day 7', 'Day 10', 'Day 14'))
# --- change labels in Day 14 Trabecular LV and \nendocardium  to be endocardium
latent_adata$obs$true_labels = as.character(latent_adata$obs$true_labels)
latent_adata$obs$true_labels[latent_adata$obs$timepoint_original == 'Day 14' & latent_adata$obs$true_labels == 'Trabecular LV and \nendocardium'] = 'Endocardium'
latent_adata$obs$true_labels = factor(latent_adata$obs$true_labels, levels=unique(latent_adata$obs$true_labels))

plot_pred_label(latent_adata, color_col='true_labels', file_name=file.path(spatial_layout_dir, 'Ground_Truth.pdf'), tp_color_palette=customized_color_palette, height=4, width=13)


# --- customize color palette for spaDOT prediction
# set color palette for different timepoint and plot spaDOT prediction
# --- re-order prediced labels
cluster_mapping = list(
    'Day 4'=c("0"="1", "1"="3", "2"="5", "3"="4", "4"="2"),
    'Day 7'=c("0"="4", "1"="3-1", "2"="3-3", "3"="6", "4"="3-2", "5"="5", "6"="1-1"),
    'Day 10'=c("0"="3-2", "1"="1-1", "2"="2", "3"="5", "4"="3-3", "5"="4", "6"="3-1"),
    'Day 14'=c("0"="3-1", "1"="5-1", "2"="4", "3"="1-1", "4"="5-2", "5"="3-3")
)
latent_adata$obs = latent_adata$obs %>% mutate(
    spaDOT_pred_labels = case_when(
        timepoint_original == 'Day 4' ~ cluster_mapping[['Day 4']][as.character(pred_labels)],
        timepoint_original == 'Day 7' ~ cluster_mapping[['Day 7']][as.character(pred_labels)],
        timepoint_original == 'Day 10' ~ cluster_mapping[['Day 10']][as.character(pred_labels)],
        timepoint_original == 'Day 14' ~ cluster_mapping[['Day 14']][as.character(pred_labels)]
    )
)

tp_color_palette = c("1"=unname(customized_color_palette['Epicardium- like']),
"1-1"=unname(customized_color_palette['Epicardium- like']),
"2"=unname(customized_color_palette['Outflow tract']),
"3"=unname(customized_color_palette['Ventricle']),
"3-1"=unname(customized_color_palette['Compact LV and \ninter-ventricular septum']),
"3-2"=unname(customized_color_palette['Trabecular LV and \nendocardium']),
"3-3"=unname(customized_color_palette['Right ventricle']),
"4"=unname(customized_color_palette['Atria']),
"5"=unname(customized_color_palette['Valves']),
"5-1"="#ffe5ec",
# "5-2"=unname(customized_color_palette['Valves']),
"5-2"="#fb6f92",
"6"=unname(customized_color_palette['Endothelium'])
)
plot_pred_label(latent_adata, color_col='spaDOT_pred_labels', file_name=file.path(spatial_layout_dir, 'spaDOT_prediction.pdf'), tp_color_palette=tp_color_palette, legend_ncol=4, height=3.5, width=13)

# latent_adata$obs$spaDOT_pred_labels = paste(latent_adata$obs$timepoint_original, latent_adata$obs$pred_labels, sep='_')
# tp_color_palette = list(
#     "Day 4"=c("0"=unname(customized_color_palette['Epicardium- like']), 
#               "1"=unname(customized_color_palette['Ventricle']),
#               "2"=unname(customized_color_palette['Valves']),
#               "3"=unname(customized_color_palette['Atria']),
#               "4"=unname(customized_color_palette['Outflow tract'])),
#     "Day 7"=c("0"=unname(customized_color_palette['Atria']), 
#               "1"=unname(customized_color_palette['Compact LV and \ninter-ventricular septum']),
#               "2"=unname(customized_color_palette['Right ventricle']), 
#               "3"=unname(customized_color_palette['Endothelium']), 
#               "4"=unname(customized_color_palette['Trabecular LV and \nendocardium']), 
#               "5"=unname(customized_color_palette['Valves']), 
#               "6"=unname(customized_color_palette['Epicardium'])), 
#     "Day 10"=c("0"=unname(customized_color_palette['Trabecular LV and \nendocardium']), 
#                 "1"=unname(customized_color_palette['Epicardium']), 
#                 "2"=unname(customized_color_palette['Outflow tract']), 
#                 "3"=unname(customized_color_palette['Valves']), 
#                 "4"=unname(customized_color_palette['Right ventricle']), 
#                 "5"=unname(customized_color_palette['Atria']), 
#                 "6"=unname(customized_color_palette['Compact LV and \ninter-ventricular septum'])),
#     "Day 14"=c("0"=unname(customized_color_palette['Compact LV and \ninter-ventricular septum']), 
#                "1"="#E6A0C4",
#                "2"=unname(customized_color_palette['Atria']), 
#                "3"=unname(customized_color_palette['Epicardium']), 
#                "4"=unname(customized_color_palette['Valves']), 
#                "5"=unname(customized_color_palette['Right ventricle']))
# )
# plot_pred_label(latent_adata, color_col='pred_labels', file_name=file.path(spatial_layout_dir, 'spaDOT_prediction.png'), tp_color_palette=tp_color_palette, legend_ncol=4, height=3.5, width=13)

# --- plot OT heatmap
timepoints = as.vector(unique(latent_adata$obs$timepoint))
timepoints = sort(timepoints)
plot_OT_heatmap(timepoints, 'spaDOT', file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ChickenHeart_results'), tp_color_palette=tp_color_palette)

# --- load all methods's results
# spatialPCA
spatialPCA_latent = read_h5ad(file.path(result_dir, 'explore_spatialPCA', 'ChickenHeart_results', 'prediction_adata.h5ad'))
spatialPCA_latent$obs_names = paste0(spatialPCA_latent$obs_names, '_', spatialPCA_latent$obs$timepoint)
latent_adata$obs$spatialPCA_pred_labels = spatialPCA_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spatialPCA_pred_labels', file_name=file.path(spatial_layout_dir, 'spatialPCA_prediction.png'))

# GraphST
GraphST_latent = read_h5ad(file.path(result_dir, 'explore_GraphST', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$GraphST_pred_labels = GraphST_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='GraphST_pred_labels', file_name=file.path(spatial_layout_dir, 'GraphST_prediction.png'))

# SEDR
SEDR_latent = read_h5ad(file.path(result_dir, 'explore_SEDR', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$SEDR_pred_labels = SEDR_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='SEDR_pred_labels', file_name=file.path(spatial_layout_dir, 'SEDR_prediction.png'))

# spaGCN
spaGCN_latent = read_h5ad(file.path(result_dir, 'explore_spaGCN', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$spaGCN_pred_labels = spaGCN_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spaGCN_pred_labels', file_name=file.path(spatial_layout_dir, 'spaGCN_prediction.png'))

# spaVAE
spaVAE_latent = read_h5ad(file.path(result_dir, 'explore_spaVAE', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$spaVAE_pred_labels = spaVAE_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spaVAE_pred_labels', file_name=file.path(spatial_layout_dir, 'spaVAE_prediction.png'))

# STAligner
STAligner_latent = read_h5ad(file.path(result_dir, 'explore_STAligner', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$STAligner_pred_labels = STAligner_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='STAligner_pred_labels', file_name=file.path(spatial_layout_dir, 'STAligner_prediction.png'))

# PRECAST
PRECAST_latent = read_h5ad(file.path(result_dir, 'explore_PRECAST', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$PRECAST_pred_labels = PRECAST_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='PRECAST_pred_labels', file_name=file.path(spatial_layout_dir, 'PRECAST_prediction.png'))

# spatiAlign
spatiAlign_latent = read_h5ad(file.path(result_dir, 'explore_spatiAlign', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$spatiAlign_pred_labels = spatiAlign_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spatiAlign_pred_labels', file_name=file.path(spatial_layout_dir, 'spatiAlign_prediction.png'))

# SPIRAL
SPIRAL_latent = read_h5ad(file.path(result_dir, 'explore_SPIRAL', 'ChickenHeart_results', 'prediction_adata.h5ad'))
latent_adata$obs$SPIRAL_pred_labels = SPIRAL_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='SPIRAL_pred_labels', file_name=file.path(spatial_layout_dir, 'SPIRAL_prediction.png'))

# STADIA
STADIA_latent = read_h5ad(file.path(result_dir, 'explore_STADIA', 'ChickenHeart_results', 'prediction_adata.h5ad'))
STADIA_latent$obs_names = sub("^[^.]*\\.", "", STADIA_latent$obs_names)
latent_adata$obs$STADIA_pred_labels = STADIA_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='STADIA_pred_labels', file_name=file.path(spatial_layout_dir, 'STADIA_prediction.png'))

# load SVGP results
SVGP_latent = read.csv(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ChickenHeart_results', 'SVGPVAE_predictions.csv'), header=T, row.names=1)
latent_adata$obs$SVGP_pred_labels = SVGP_latent[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$SVGP_pred_labels = as.factor(latent_adata$obs$SVGP_pred_labels)
plot_pred_label(latent_adata, color_col='SVGP_pred_labels', file_name=file.path(spatial_layout_dir, 'spaDOT_SVGP_prediction.png'))

GAT_pred_labels = read.csv(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ChickenHeart_results', 'GATVAE_predictions.csv'), header=T, row.names=1)
latent_adata$obs$GAT_pred_labels = GAT_pred_labels[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$GAT_pred_labels = as.factor(latent_adata$obs$GAT_pred_labels)
plot_pred_label(latent_adata, color_col='GAT_pred_labels', file_name=file.path(spatial_layout_dir, 'spaDOT_GAT_prediction.png'))

# --- compute lisi
domain_metrics_df = read.csv(file.path(dirname(spatial_layout_dir), "final_ChickenHeart_domain_metrics.csv"), header=T)
domain_metrics_df$X = NULL
domain_metrics_df$method[domain_metrics_df$method == 'GATVAE-SVGP-VAE'] = 'SpaDOT'
methods = unique(domain_metrics_df$method)
# method_color_palette = c("#F21A00", "#F3DF6C", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#446455")
method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae")
customized_method_colors = method_color_palette[1:length(unique(domain_metrics_df$method))]
names(customized_method_colors) = unique(domain_metrics_df$method)

library(lisi)
timepoints = as.vector(unique(latent_adata$obs$timepoint))
timepoints = sort(timepoints)
for (tp in timepoints) {
    tp_latent_adata = latent_adata[latent_adata$obs$timepoint == tp, ]
    locations = data.frame(x=tp_latent_adata$obs$pixel_x, y=tp_latent_adata$obs$pixel_y)
    for (method in methods) {
        if (method == 'SpaDOT') {
            lisi = compute_lisi(locations, tp_latent_adata$obs, c('pred_labels'))
            median_lisi = median(lisi[, 1])
            domain_metrics_df[nrow(domain_metrics_df)+1,] = c('median-LISI', tp, median_lisi, method)
        }
        else {
            lisi = compute_lisi(locations, tp_latent_adata$obs, c(paste0(method,'_pred_labels')))
            median_lisi = median(lisi[, 1])
            domain_metrics_df[nrow(domain_metrics_df)+1,] = c('median-LISI', tp, median_lisi, method)
        }
    }
}
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
write.csv(domain_metrics_df, file=file.path(dirname(spatial_layout_dir), 'SupplementaryTable5_metric.csv'), row.names=F)

# --- draw radar plot
library(fmsb)
library(tidyr)
library(dplyr)
for (tp in unique(domain_metrics_df$timepoint)) {
    tp_domain_metrics_df = domain_metrics_df[domain_metrics_df$timepoint == tp,]
    tp_domain_metrics_df_wide = tp_domain_metrics_df %>% pivot_wider(names_from=index, values_from=value)
    tp_domain_metrics_df_wide[['1-PAS']] = 1 - as.numeric(tp_domain_metrics_df_wide[['PAS']])
    # tp_domain_metrics_df_wide[['1-CHAOS']] = 1 - as.numeric(tp_domain_metrics_df_wide[['CHAOS']])
    # --- edit to normalized CHAOS
    tp_domain_metrics_df_wide[['norm_CHAOS']] = as.numeric(tp_domain_metrics_df_wide[['CHAOS']])/max(as.numeric(tp_domain_metrics_df_wide[['CHAOS']]))
    tp_domain_metrics_df_wide[['1-norm_CHAOS']] = 1 - tp_domain_metrics_df_wide[['norm_CHAOS']]
    tp_domain_metrics_df_wide[['norm_median-LISI']] = as.numeric(tp_domain_metrics_df_wide[['median-LISI']])/max(as.numeric(tp_domain_metrics_df_wide[['median-LISI']]))
    tp_domain_metrics_df_wide[['1-norm_LISI']] = 1 - tp_domain_metrics_df_wide[['norm_median-LISI']]
    tp_domain_metrics_df_wide = as.data.frame(tp_domain_metrics_df_wide[, c('method', 'ARI', 'AMI', 'H', '1-norm_LISI', '1-PAS', '1-norm_CHAOS')])
    tp_domain_metrics_df_wide$ARI = as.numeric(tp_domain_metrics_df_wide$ARI)
    tp_domain_metrics_df_wide$AMI = as.numeric(tp_domain_metrics_df_wide$AMI)
    tp_domain_metrics_df_wide$H = as.numeric(tp_domain_metrics_df_wide$H)
    rownames(tp_domain_metrics_df_wide) = tp_domain_metrics_df_wide$method
    tp_domain_metrics_df_wide$method = NULL
    # normalize each column
    tp_domain_metrics_df_wide = apply(tp_domain_metrics_df_wide, 2, function(x) (x-min(x))/(max(x)-min(x)))
    tp_domain_metrics_df_wide = rbind(rep(1,6) , rep(0,6) , tp_domain_metrics_df_wide)
    colors_in = alpha(unname(customized_method_colors),0.8)
    if (tp == 'avg') {
        png(file.path(downstream_dir, paste0(tp, '_radar.png')), width=800, height=1200)
    } else {
        png(file.path(downstream_dir, paste0(tp, '_radar.png')), width=800, height=800)
    }
    par(mar = c(0.1, 0.1, 0.1, 0.1)) 
    radarchart(as.data.frame(tp_domain_metrics_df_wide), axistype=1, pcol=colors_in, #pfcol=unname(customized_method_colors),
    plwd=4, plty=1, cglcol="grey", cglty=1, axislabcol="grey", caxislabels=c(0, 0.25, 0.5, 0.75, 1.0),
    cglwd=0.8, vlcex=1.8)
    if (tp == 'avg') {
        methods = as.character(levels(domain_metrics_df$method))
        legend(x=-1, y=-1.2, legend = methods, bty = "n", pch=20 , col=unname(customized_method_colors[methods]), text.col = "black", cex=2, pt.cex=3, ncol=4)
    }
    dev.off()
}


domain_metrics_df$value = as.numeric(domain_metrics_df$value)
for (metric in unique(domain_metrics_df$index)) {
    metric_df = domain_metrics_df[domain_metrics_df$index == metric, ]
    avg_metric_df = metric_df[metric_df$timepoint == 'avg', ]
    timepoint_metric_df = metric_df[metric_df$timepoint != 'avg', ]
    p = ggplot() +
        geom_point(data=avg_metric_df, aes(x=method, y=value, color=method, group=method), size=6) +
        scale_color_manual(values=customized_method_colors) + 
        geom_point(data=timepoint_metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.5) +
        theme_pubr() +
        theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric)
    p = p + guides(color = guide_legend(nrow = 4, title=NULL), shape=guide_legend(nrow=4, title=NULL))
    legend = get_plot_component(p, "guide-box", return_all=TRUE)[[4]]
    ggsave(file.path(dirname(spatial_layout_dir), 'legend.pdf'), plot=legend, width=6, height=2)
    p = p + theme(legend.position="none")
    ggsave(file.path(dirname(spatial_layout_dir), paste0(metric, ".pdf")), plot=p, width=4, height=3.5)
}




## derive cluster specific genes
data_dir = file.path('/media/sda/wenjinma/data/spatial_ChickenHeart')
adata = read_h5ad(file.path(data_dir, 'processed_ChickenHeart.h5ad'))
common_cells = intersect(rownames(adata), rownames(latent_adata))
adata = adata[common_cells]
latent_adata = latent_adata[common_cells, ]
# --- create Seurat object for analysis
Seurat_obj = CreateSeuratObject(counts=t(adata$X), meta.data=latent_adata$obs)
Seurat_obj[['spaDOT_latent']] = CreateAssayObject(data = t(latent_adata$X))
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
# Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
# Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
# Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)


# --- derive cluster specific genes
# cluster_marker_genes(Seurat_obj, 'spaDOT', method_result_dir=file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ChickenHeart_results'))

# t0_Seurat_obj = subset(Seurat_obj, subset=timepoint=='t0')
# Idents(t0_Seurat_obj) = "spaDOT_pred_labels"
# DotPlot(t0_Seurat_obj, features=c('NPPA', 'FHL2', 'GJA5', 'HBZ', 'Wpkci-7', 'SVIP', 'TGFB2', 'MRPL41', 'ACADS', 'PSMD7', 'COL1A1', 'CSRP2', 'CNMD', 'CCDC80', 'CRIP1', 'MYH7', 'MYH1D', 'PTN', 'APOA1', 'SERPINE2','MSX1', 'CHODL', 'TMSB4X', 'PENK', 'COL9A2')) + coord_flip()
# ggsave(file.path(downstream_dir, 't0_dotplot.png'), height=15, width=8)

# t1_Seurat_obj = subset(Seurat_obj, subset=timepoint=='t1')
# Idents(t1_Seurat_obj) = "spaDOT_pred_labels"
# DotPlot(t1_Seurat_obj, features=c('MYH15', 'HBM', 'HBBA', 'HBA1', 'HBE1', 'FN1', 'TBX5', 'MYOZ2', 'RARRES1', 'SLN', 'ACTG2', 'ACTA1', 'CST3', 'MYL1', 'MYH7', 'MYH1D', 'ATP2A2', 'SERPINE2', 'BMP10', 'FBLN1', 'KRT7', 'COL3A1', 'HBBR', 'HBBR.1', 'CHODL', 'CSRP2', 'S100A6', 'COL1A1', 'CRIP1')) + coord_flip()
# ggsave(file.path(downstream_dir, 't1_dotplot.png'), height=15, width=8)

# find overlap between cluster 1 and cluster 4 at Day 14
method_result_dir = file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ChickenHeart_results')
Day14_marker = read.csv(file.path(method_result_dir, paste0('t3', '_Seurat_FindAllMarkers_original_cluster.csv')), row.names=1)
Day14_sig_marker = Day14_marker[Day14_marker$p_val_adj < 0.05 & Day14_marker$avg_log2FC > 0, ]
Day14_cluster1_marker = Day14_sig_marker[Day14_sig_marker$cluster == 1, 'gene']
Day14_cluster4_marker = Day14_sig_marker[Day14_sig_marker$cluster == 4, 'gene']
go_list = list(Day14_cluster1_marker, Day14_cluster4_marker)
names(go_list) = c('Domain 1', 'Domain 4')
library(ggvenn)
g = ggvenn(go_list)
ggsave(file.path(downstream_dir, 'Day14_cluster1_vs_cluster4.png'), plot=g, height=6, width=6)
Day14_cluster1_exclusive_genes = setdiff(Day14_cluster1_marker, Day14_cluster4_marker)
write(Day14_cluster1_exclusive_genes, file=file.path(downstream_dir, 'Day14_cluster1_exclusive_genes.txt'))




Day14_cluster4_exclusive_genes = setdiff(Day14_cluster4_marker, Day14_cluster1_marker)
write(Day14_cluster4_exclusive_genes, file=file.path(downstream_dir, 'Day14_cluster4_exclusive_genes.txt'))

# --- plot certain marker genes
Day14_Seurat_obj = subset(Seurat_obj, subset=timepoint=='t3')
Idents(Day14_Seurat_obj) = "spaDOT_pred_labels"
features = c('LYZ', 'GZMA', 'TPSB2', 'KIT', 'TPSAB1', 'ACTA2', 'MYH11', 'MYL9', 'AREG', 'CPR3', 'CTSG', 'PECAM1', 'EMCN', 'COL3A1', 'TAGLN', 'RGS5', 'POSTN', 'TWIST1', 'TBX20', 'DCN', 'CDH5', 'PODXL', 'ENG', 'RARRES1', 'NPR3', 'COL1C11', 'CYP1B1', 'PLAT', 'FABP4', 'FABP5', 'CLDN5', 'HES4', 'RARRES2', 'GUCY1A1', 'CNN1', 'CRABP2', 'PENK', 'LEPR', 'GATA5', 'ID2')
DotPlot(Day14_Seurat_obj, features=features) + coord_flip()
# DotPlot(Day14_Seurat_obj, features=c('NFATC1', 'NOTCH1', 'GATA4', 'SOX9', 'COL1A1', 'COL3A1', 'POSTN', 'ACTA2', 'TAGLN', 'TBX20', 'BMP2', 'BMP4', 'HAS2', 'COL1A2', 'VCAN', 'MYH6', 'FN1', 'MGP', 'COL8A1', 'EFNA1', 'HEY1', 'DLL4', 'HMGB1', 'STAT1', 'CDH5', 'PECAM1', 'HEY2', 'HOTAIR', 'EIF2A', 'GNRH', 'NPR3', 'CDH11', 'NTRK2', 'EDN1', 'GATA5', 'MEF2C', 'SOX17', 'THSD1', 'EMCN', 'TGFB2', 'ATF3', 'AGTR1', 'MYL9', 'MYH11', 'TPM2', 'TPSB2', 'TPSAB1', 'AREG', 'CPA3', 'CTSG', 'KIT', 'LAMP1', 'LYZ', 'PLAT', 'GJA1', 'HES4', 'CRLF1', 'PLA2G2A', 'C1R', 'CXCR4', 'MSX1')) + coord_flip()
ggsave(file.path(downstream_dir, 'Day14_dotplot.png'), height=15, width=8)

# --- plot marker genes on these two domains only
# Day14_cluster1_cluster4_Seurat_obj = subset(Day14_Seurat_obj, subset=spaDOT_pred_labels %in% c('Day 14_1', 'Day 14_4'))
Day14_cluster1_cluster4_Seurat_obj = subset(Day14_Seurat_obj, subset=spaDOT_pred_labels %in% c('5-1', '5-2'))
# marker_genes = c('MYL9', 'CCL21', 'FABP4', 'TFF3', 'ACTA2', 'MYH11', 'DES', 'DCN', 'CDH19', 'COLEC11', 'HES4', 'ERBB3', 'SERPINE2', 'S100B', 'S100A1')
# marker_genes = c('LYZ', 'GZMA', 'TPSB2', 'KIT', 'TPSAB1', 'ACTA2', 'MYH11', 'MYL9', 'AREG', 'CPR3', 'CTSG', 'PECAM1', 'EMCN', 'COL3A1', 'TAGLN', 'RGS5', 'POSTN', 'CNN1', 'CRABP2', 'LEPR', 'PENK', 'GATA5', 'ATP2A2', 'NKX2-5', 'MYH1E', 'TBX5', 'TBX20', 'SMAD6', 'BMP2', 'COL1A1', 'COL1A2', 'COL2A1', 'TNC', 'BMP4', 'BMP5', 'BMP6', 'BMP7')

marker_genes = c('GATA5', 'CNN1', 'ATP2A2', 'COL1A1', 'COL1A2', 'COL3A1')
for (gene in marker_genes) {
    gene_idx = which(rownames(Day14_cluster1_cluster4_Seurat_obj) == gene)
    if (length(gene_idx) == 0) {
        cat(gene, 'is not in the Seurat object..\n')
        next
    }
    gene_exprs = Day14_cluster1_cluster4_Seurat_obj@assays$RNA@layers$data[gene_idx, ]
    gene_exprs_df = Day14_cluster1_cluster4_Seurat_obj@meta.data
    gene_exprs_df[, gene] = gene_exprs
    # explore gene expression percentage in each group
    summary_df = gene_exprs_df %>%
            group_by(spaDOT_pred_labels) %>%
            summarize(non_zero_pct = mean(!!sym(gene) > 0))
    # if (all(summary_df$non_zero_pct < 0.5)) {
    #     cat(gene, 'does not have enough expression..\n')
    #     next
    # }
    ggplot(gene_exprs_df, aes(x=spaDOT_pred_labels, y=!!sym(gene), color=spaDOT_pred_labels)) + geom_boxplot() + theme_pubr() + labs(color = "SpaDOT Domains", x="")
    ggsave(file.path(downstream_dir, paste0(gene, '_geneexprs_boxplot.pdf')),height=4, width=4)    

    ggplot(gene_exprs_df, aes(x=pixel_x, y=pixel_y, color=.data[[gene]])) + geom_point(size=0.5) +
        theme_pubr() + 
        theme(axis.text =element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) + 
        scale_color_gradient(low = "grey", high = "red")
    ggsave(file.path(downstream_dir, paste0(gene, '_geneexprs.pdf')),height=4, width=2.5)    
}


# --- perform GO
GO_analysis(timepoints, 'spaDOT', file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ChickenHeart_results'), tp_color_palette=tp_color_palette)

# --- overlap with cluster 4 at Day14
Day14_Seurat_obj = subset(Seurat_obj, subset=timepoint=='t3')
overlap_df = data.frame(valve.cells=Day14_Seurat_obj$Valve.cells)
overlap_df$valve_scores = ifelse(overlap_df$valve.cells > 0.5, overlap_df$valve.cells, 0)
overlap_df$valve_labels = ifelse(overlap_df$valve.cells > 0.5, TRUE, FALSE)
overlap_df$pixel_x = Day14_Seurat_obj$pixel_x
overlap_df$pixel_y = Day14_Seurat_obj$pixel_y
ggplot(overlap_df, aes(x=pixel_x, y=pixel_y)) + 
    geom_point(aes(color=valve.cells), size=0.1, alpha=0.6) +
    scale_color_gradientn(colours=c("grey","red", "purple"),
    values=c(0, 0.5, 1), limits=c(0, 1)) + 
    # breaks=c(0, 0.5, 1), limits=c(0, 1)) + 
    # scale_color_gradientn(colors = c("yellow", "red"),
    # values = c(0, 0.5, 1),
    # limits = c(0.5, 1), na.value='grey') + 
    theme_void() + 
    theme(plot.background = element_rect(fill = "white"),
    legend.position = "bottom", legend.title=element_blank())
ggsave(file.path(downstream_dir, 'Day14_valve_cells.png'), height=3, width=2)

spaDOT_valve_labels = ifelse(Day14_Seurat_obj$spaDOT_pred_labels == 'Day 14_4', TRUE, FALSE)
STADIA_valve_labels = ifelse(Day14_Seurat_obj$STADIA_pred_labels == 2, TRUE, FALSE)
SEDR_valve_labels = ifelse(Day14_Seurat_obj$SEDR_pred_labels == 3, TRUE, FALSE)
SPIRAL_valve_labels = ifelse(Day14_Seurat_obj$SPIRAL_pred_labels == 2, TRUE, FALSE)
spaGCN_valve_labels = ifelse(Day14_Seurat_obj$spaGCN_pred_labels == 1, TRUE, FALSE)
spaVAE_valve_labels = ifelse(Day14_Seurat_obj$spaVAE_pred_labels == 2, TRUE, FALSE)
STAligner_valve_labels = ifelse(Day14_Seurat_obj$STAligner_pred_labels == 2, TRUE, FALSE)
spatialPCA_valve_labels = ifelse(Day14_Seurat_obj$spatialPCA_pred_labels == 3, TRUE, FALSE)
PRECAST_valve_labels = ifelse(Day14_Seurat_obj$PRECAST_pred_labels == 3, TRUE, FALSE)
GraphST_valve_labels = ifelse(Day14_Seurat_obj$GraphST_pred_labels == 2, TRUE, FALSE)
spatiAlign_valve_labels = ifelse(Day14_Seurat_obj$spatiAlign_pred_labels == 4, TRUE, FALSE)


overlap_df$spaDOT = spaDOT_valve_labels
overlap_df$STADIA = STADIA_valve_labels
overlap_df$SEDR = SEDR_valve_labels
overlap_df$SPIRAL = SPIRAL_valve_labels
overlap_df$spaGCN = spaGCN_valve_labels
overlap_df$spaVAE = spaVAE_valve_labels
overlap_df$STAligner = STAligner_valve_labels
overlap_df$spatialPCA = spatialPCA_valve_labels
overlap_df$PRECAST = PRECAST_valve_labels
overlap_df$GraphST = GraphST_valve_labels
overlap_df$spatiAlign = spatiAlign_valve_labels

library(caret)
spaDOT_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$spaDOT))$byClass[7])
STADIA_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$STADIA))$byClass[7])
SEDR_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$SEDR))$byClass[7])
SPIRAL_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$SPIRAL))$byClass[7])
spaGCN_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$spaGCN))$byClass[7])
spaVAE_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$spaVAE))$byClass[7])
STAligner_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$STAligner))$byClass[7])
spatialPCA_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$spatialPCA))$byClass[7])
PRECAST_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$PRECAST))$byClass[7])
GraphST_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$GraphST))$byClass[7])
spatiAlign_pred_f1 = unname(confusionMatrix(as.factor(overlap_df$valve_labels), as.factor(overlap_df$spatiAlign))$byClass[7])

f1_df = data.frame(method=c('SpaDOT', 'STADIA', 'SEDR', 'SPIRAL', 'spaGCN', 'spaVAE', 'STAligner', 'spatialPCA', 'PRECAST', 'GraphST', 'spatiAlign'),
                    f1=c(spaDOT_pred_f1, STADIA_pred_f1, SEDR_pred_f1, SPIRAL_pred_f1, spaGCN_pred_f1, spaVAE_pred_f1, STAligner_pred_f1, spatialPCA_pred_f1, PRECAST_pred_f1, GraphST_pred_f1, spatiAlign_pred_f1))
f1_df$f1 = as.numeric(f1_df$f1)
f1_df$method = factor(f1_df$method, levels=unique(f1_df$method))

customized_method_colors['SpaDOT'] = '#F3DF6C'
ggplot(f1_df, aes(x=method, y=f1, fill=method)) + geom_bar(stat='identity') + coord_cartesian(ylim=c(0.95,1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab('F1 score') + xlab('Method') + scale_fill_manual(values=customized_method_colors)
ggsave(file.path(downstream_dir, 'valve_cell_F1.png'), height=4, width=6)




# library(PRROC)
# spaDOT_pr_curve = pr.curve(scores.class0 = overlap_df$valve.cells[overlap_df$spaDOT == TRUE], 
#                      scores.class1 = overlap_df$valve.cells[overlap_df$spaDOT == FALSE], curve = TRUE)


# STADIA_pr_curve = pr.curve(scores.class0 = overlap_df$valve.cells[overlap_df$STADIA == TRUE], 
#                      scores.class1 = overlap_df$valve.cells[overlap_df$STADIA == FALSE], curve = TRUE)
# png(file.path(downstream_dir, 'spaDOT_AUPR.png'))
# plot(spaDOT_pr_curve)
# dev.off()
# png(file.path(downstream_dir, 'STADIA_AUPR.png'))
# plot(STADIA_pr_curve)
# dev.off()









# Day14_Seurat_obj = subset(Seurat_obj, subset=timepoint=='t3')

# # Day14_valve_scores = ifelse(Day14_Seurat_obj$Valve.cells > 0.5, TRUE, FALSE)
# # Day14_valve_labels = ifelse(Day14_Seurat_obj$spaDOT_pred_labels == 'Day 14_4', TRUE, FALSE)
# Day14_spaDOT_valve_labels = ifelse(Day14_Seurat_obj$spaDOT_pred_labels == 'Day 14_4', TRUE, FALSE)
# # sum(Day14_valve_scores & Day14_spaDOT_valve_labels) / sum(Day14_valve_scores | Day14_spaDOT_valve_labels) # 0.112

# Day14_STADIA_valve_labels = ifelse(Day14_Seurat_obj$STADIA_pred_labels == 2, TRUE, FALSE)
# # sum(Day14_valve_scores & Day14_STADIA_valve_labels) / sum(Day14_valve_scores | Day14_STADIA_valve_labels) # 0.084

# table(Day14_Seurat_obj$Valve.cells >= 0.5, Day14_spaDOT_valve_labels)
# table(Day14_Seurat_obj$Valve.cells >= 0.5, Day14_STADIA_valve_labels)


# roc(Day14_valve_labels, Day14_Seurat_obj$Valve.cells) # 0.959
# roc(Day14_valve_labels, Day14_valve_scores) # 0.556
# STADIA_valve_labels = ifelse(Day14_Seurat_obj$STADIA_pred_labels == 2, TRUE, FALSE)
# roc(Day14_STADIA_valve_labels, Day14_Seurat_obj$Valve.cells) # 0.975
# roc(Day14_STADIA_valve_labels, Day14_valve_scores) # 0.542




# plot Day 14 valve.cells

# --- perform DEGs between t3 1 and 4
Idents(Seurat_obj) = "spaDOT_pred_labels"
# Day14_1_vs_4_DE_genes = FindMarkers(Seurat_obj, ident.1 = "Day 14_1", ident.2 = "Day 14_4") # here, we want to consider 4 as baseline
Day14_1_vs_4_DE_genes = FindMarkers(Seurat_obj, ident.1 = "5-1", ident.2 = "5-2")

# volcano plot
Day14_1_vs_4_DE_genes$negLog10adjP = -log10(Day14_1_vs_4_DE_genes$p_val_adj)
Day14_1_vs_4_DE_genes$significant = ifelse(Day14_1_vs_4_DE_genes$p_val_adj < 0.05 & abs(Day14_1_vs_4_DE_genes$avg_log2FC) > 1, 'Significant', 'Not Significant')
Day14_1_vs_4_DE_genes = Day14_1_vs_4_DE_genes[order(Day14_1_vs_4_DE_genes$p_val_adj),]
# top 20 genes with smallest p-value and avg_log2FC > 1
# labeled_genes = rownames(Day14_1_vs_4_DE_genes)[which(abs(Day14_1_vs_4_DE_genes$avg_log2FC) > 1)[1:20]]
labeled_genes = c('COL1A1', 'COL1A2', 'ACTC1', 'MYL2', 'MYL3', 'S100A6', 'S100A10', 'SPARC')

Day14_1_vs_4_DE_genes$label = ifelse(rownames(Day14_1_vs_4_DE_genes) %in% labeled_genes, rownames(Day14_1_vs_4_DE_genes), NA)
ggplot(Day14_1_vs_4_DE_genes, aes(x = avg_log2FC, y = negLog10adjP)) +
    geom_point(aes(color = significant), size = 2) +
    geom_point(
        data = Day14_1_vs_4_DE_genes[!is.na(Day14_1_vs_4_DE_genes$label), ], # Only for labeled genes
        aes(fill = significant), 
        shape = 21, # Circle with outline
        size = 3,   # Slightly larger size for outlined points
        color = "black" # Outline color
    ) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.title = element_blank()
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
    geom_text_repel(aes(label = label), 
                  box.padding = 0.35, 
                  point.padding = 0.3, 
                  segment.color = "grey50", # Color of the connecting line
                  segment.size = 0.5)  # Adjust `max.overlaps` as needed
ggsave(file.path(downstream_dir, 'Day14_1_vs_4_DEGs_volcano.pdf'), height=5, width=6)

# obtain upregulated genes
Day14_1_vs_4_up_genes = rownames(Day14_1_vs_4_DE_genes[Day14_1_vs_4_DE_genes$avg_log2FC > 1 & Day14_1_vs_4_DE_genes$p_val_adj < 0.05,])
library(org.Gg.eg.db)
Day14_1_vs_4_up_genes_GO = enrichGO(Day14_1_vs_4_up_genes, keyType = 'SYMBOL',
                                OrgDb = org.Gg.eg.db, 
                                ont = "BP", # Biological process
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)
g = dotplot(Day14_1_vs_4_up_genes_GO, showCategory=10)
ggsave(file.path(downstream_dir,'Day14_1_vs_4_UP_enrichment_clusterprofiler.png'), width=8, height=6)
write(Day14_1_vs_4_up_genes, file=file.path(downstream_dir, 'Day14_1_vs_4_upregulated_genes.txt'))

# down regulated genes
Day14_1_vs_4_down_genes = rownames(Day14_1_vs_4_DE_genes[Day14_1_vs_4_DE_genes$avg_log2FC < -1 & Day14_1_vs_4_DE_genes$p_val_adj < 0.05,])
Day14_1_vs_4_down_genes_GO = enrichGO(Day14_1_vs_4_down_genes, keyType = 'SYMBOL',
                                OrgDb = org.Gg.eg.db, 
                                ont = "BP", # Biological process
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)
g = dotplot(Day14_1_vs_4_down_genes_GO, showCategory=10)
ggsave(file.path(downstream_dir,'Day14_1_vs_4_DOWN_enrichment_clusterprofiler.png'), width=8, height=6)
write(Day14_1_vs_4_down_genes, file=file.path(downstream_dir, 'Day14_1_vs_4_downregulated_genes.txt'))

# plot overall gene counts
ggplot(Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=nCount_RNA)) + geom_point(size=0.5) +
        facet_wrap(~timepoint_original, scale='free') + 
        theme_pubr() + 
        theme(axis.text =element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) + 
        scale_color_gradient(low = "grey", high = "red")
ggsave(file.path(downstream_dir, 'nCount_RNA_geneexprs.png'),height=4, width=10) 

ggplot(Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=nFeature_RNA)) + geom_point(size=0.5) +
        facet_wrap(~timepoint_original, scale='free') + 
        theme_pubr() + 
        theme(axis.text =element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) + 
        scale_color_gradient(low = "grey", high = "red")
ggsave(file.path(downstream_dir, 'nFeature_RNA_geneexprs.png'),height=4, width=10) 


# --- plot marker gene enrichment
ventricle_genes = c('TBX5', 'ACTG2', 'CHGB', 'TBX20', 'MB', 'NPPA')
names(ventricle_genes) = rep('Ventricle', length(ventricle_genes))
marker_genes = c(ventricle_genes)
df = data.frame(matrix(nrow=0, ncol=5))
    colnames(df) = c('timepoint', 'cluster', 'gene', 'region', 'enrichment')
# subset Seurat object to ventricles only
Day7_Seurat_obj = subset(Seurat_obj, subset=spaDOT_pred_labels %in% c('Day 7_1', 'Day 7_2', 'Day 7_4'))
Day7_2_cluster_cells_idx = Idents(Day7_Seurat_obj) == 'Day 7_2'
Day7_1_cluster_cells_idx = Idents(Day7_Seurat_obj) == 'Day 7_1'
Day7_4_cluster_cells_idx = Idents(Day7_Seurat_obj) == 'Day 7_4'
for (gene_idx in 1:length(marker_genes)) {
    gene = marker_genes[gene_idx]
    Day7_2_avg_exprs = mean(Day7_Seurat_obj[unname(gene), Day7_2_cluster_cells_idx]$RNA@layers$counts)
    Day7_1_avg_exprs = mean(Day7_Seurat_obj[unname(gene), Day7_1_cluster_cells_idx]$RNA@layers$counts)
    Day7_4_avg_exprs = mean(Day7_Seurat_obj[unname(gene), Day7_4_cluster_cells_idx]$RNA@layers$counts)
    Day7_1_vs_2_enrichment = (Day7_1_avg_exprs-Day7_2_avg_exprs)/Day7_2_avg_exprs
    Day7_4_vs_2_enrichment = (Day7_4_avg_exprs-Day7_2_avg_exprs)/Day7_2_avg_exprs
    df[nrow(df) + 1,] = c('Day 7', 'Compact', unname(gene), names(gene), Day7_1_vs_2_enrichment)
    df[nrow(df) + 1,] = c('Day 7', 'Trabecular', unname(gene), names(gene), Day7_4_vs_2_enrichment)
}

Day10_Seurat_obj = subset(Seurat_obj, subset=spaDOT_pred_labels %in% c('Day 10_0', 'Day 10_4', 'Day 10_6'))
Day10_0_cluster_cells_idx = Idents(Day10_Seurat_obj) == 'Day 10_0'
Day10_4_cluster_cells_idx = Idents(Day10_Seurat_obj) == 'Day 10_4'
Day10_6_cluster_cells_idx = Idents(Day10_Seurat_obj) == 'Day 10_6'
for (gene_idx in 1:length(marker_genes)) {
    gene = marker_genes[gene_idx]
    Day10_0_avg_exprs = mean(Day10_Seurat_obj[unname(gene), Day10_0_cluster_cells_idx]$RNA@layers$counts)
    Day10_4_avg_exprs = mean(Day10_Seurat_obj[unname(gene), Day10_4_cluster_cells_idx]$RNA@layers$counts)
    Day10_6_avg_exprs = mean(Day10_Seurat_obj[unname(gene), Day10_6_cluster_cells_idx]$RNA@layers$counts)
    Day10_0_vs_4_enrichment = (Day10_0_avg_exprs-Day10_4_avg_exprs)/Day10_4_avg_exprs
    Day10_6_vs_4_enrichment = (Day10_6_avg_exprs-Day10_4_avg_exprs)/Day10_4_avg_exprs
    df[nrow(df) + 1,] = c('Day 10', 'Trabecular', unname(gene), names(gene), Day10_0_vs_4_enrichment)
    df[nrow(df) + 1,] = c('Day 10', 'Compact', unname(gene), names(gene), Day10_6_vs_4_enrichment)
}

Day14_Seurat_obj = subset(Seurat_obj, subset=spaDOT_pred_labels %in% c('Day 14_0', 'Day 14_1', 'Day 14_5'))
Day14_0_cluster_cells_idx = Idents(Day14_Seurat_obj) == 'Day 14_0'
Day14_1_cluster_cells_idx = Idents(Day14_Seurat_obj) == 'Day 14_1'
Day14_5_cluster_cells_idx = Idents(Day14_Seurat_obj) == 'Day 14_5'
for (gene_idx in 1:length(marker_genes)) {
    gene = marker_genes[gene_idx]
    Day14_0_avg_exprs = mean(Day14_Seurat_obj[unname(gene), Day14_0_cluster_cells_idx]$RNA@layers$counts)
    Day14_1_avg_exprs = mean(Day14_Seurat_obj[unname(gene), Day14_1_cluster_cells_idx]$RNA@layers$counts)
    Day14_5_avg_exprs = mean(Day14_Seurat_obj[unname(gene), Day14_5_cluster_cells_idx]$RNA@layers$counts)
    Day14_0_vs_5_enrichment = (Day14_0_avg_exprs-Day14_5_avg_exprs)/Day14_5_avg_exprs
    df[nrow(df) + 1,] = c('Day 14', 'Compact', unname(gene), names(gene), Day14_0_vs_5_enrichment)
}

# write.csv(df, file.path(downstream_dir, 'cluster_gene_enrichment_compared_to_right_ventricle.csv'), quote=F)


df = read.csv(file.path(downstream_dir, 'cluster_gene_enrichment_compared_to_right_ventricle.csv'), row.names=1)
df$enrichment = as.numeric(df$enrichment)
df$cluster = as.factor(df$cluster)
df = df %>% mutate(domain=ifelse(cluster == 'Compact', 'Compact LV', 'Trabecular LV'))
df$timepoint = factor(df$timepoint, levels=c('Day 7', 'Day 10', 'Day 14'))
df$gene = factor(df$gene, levels=c('TBX5', 'ACTG2','NPPA' , 'TBX20', 'MB', 'CHGB'))

ggplot(data=df, aes(x=timepoint, y=enrichment, group=domain)) +
geom_line(aes(color=domain, linetype=domain)) +
geom_hline(yintercept = 0, linetype = "dotted", color = "red") +  # Add red dotted line
facet_wrap(~gene) +
scale_color_manual(values=c('Compact LV'='#90e0ef', 'Trabecular LV'='#caf0f8')) +  
theme_pubr() + 
theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave(file.path(downstream_dir, 'gene_enrichment_comparedto_rightventricle.pdf'), height=4, width=6)