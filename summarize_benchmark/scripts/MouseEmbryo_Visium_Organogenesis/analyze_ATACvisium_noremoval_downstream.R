Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(ClusterR)
# library(fgsea)
# library(msigdbr)
library(dplyr)
# C5_gene_sets = msigdbr(species = "mouse", category = "C5", subcategory = 'GO:BP') # subcategory=
# msigdbr_list = split(x = C5_gene_sets$gene_symbol, f = C5_gene_sets$gs_name)
library(clusterProfiler)
library(org.Mm.eg.db)
source('ATACVisium_utils.R')

base_dir = "/net/mulan/home/wenjinma"
result_dir = file.path(base_dir, "projects/SpaDOT_analysis")
spatial_layout_dir = file.path(result_dir, "summarize_benchmark/figures/ATAC-Visium_noremoval_figures/spatial_layout")
dir.create(spatial_layout_dir, showWarnings = FALSE, recursive = TRUE)
downstream_dir = file.path(result_dir, "summarize_benchmark/figures/ATAC-Visium_noremoval_figures/ATAC-Visium_noremoval_downstream")
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)
latent_adata = read_h5ad(file.path(result_dir, 'explore_SpaDOT', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
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

# ground truth
plot_pred_label(latent_adata$obs, color_col='true_labels', color_legend_title='Region Annotation', file_name=file.path(spatial_layout_dir, 'ground_truth.pdf'), width=9, height=3)


# set color palette for different timepoint and plot spaDOT prediction
# latent_adata$obs$spaDOT_pred_labels = paste(latent_adata$obs$timepoint_original, latent_adata$obs$pred_labels, sep='_')
# --- edit predicted labels
E125_cluster_mapping = c("0"="5", "1"="1", "2"="4", "3"="3", "4"="2")
E135_cluster_mapping = c("0"="4", "1"="1", "2"="3", "3"="5", "4"="2")
E155_cluster_mapping = c("0"="5-1", "1"="1", "2"="5-2", "3"="4", "4"="2")
latent_adata$obs = latent_adata$obs %>% mutate(
    spaDOT_pred_labels = case_when(
        timepoint_original == 'E12.5' ~ E125_cluster_mapping[as.character(pred_labels)],
        timepoint_original == 'E13.5' ~ E135_cluster_mapping[as.character(pred_labels)],
        timepoint_original == 'E15.5' ~ E155_cluster_mapping[as.character(pred_labels)]
    )
)
tp_color_palette = c("1"="#F8766D", "2"="#B79F00", "3"="#00BA38", "4"="#4EAEF0", "5"="#D772EC", "5-1"="#a06cd5", "5-2"="#e0aaff")


# tp_color_palette = c("E12.5_1"="#F8766D", "E12.5_4"="#B79F00", "E12.5_3"="#00BA38", "E12.5_2"="#4EAEF0", "E12.5_0"="#D772EC", "E13.5_0"="#4EAEF0", "E13.5_1"="#F8766D", "E13.5_2"="#00BA38", "E13.5_3"="#D772EC", "E13.5_4"="#B79F00", "E15.5_1"="#F8766D", "E15.5_3"="#4EAEF0", "E15.5_4"="#B79F00", "E15.5_0"="#990099", "E15.5_2"="#FF99FF")
plot_pred_label(latent_adata$obs, color_col='spaDOT_pred_labels', color_legend_title='SpaDOT Domains', file_name=file.path(spatial_layout_dir, 'spaDOT_prediction.pdf'), tp_color_palette=tp_color_palette, height=4.5)

# --- plot OT
timepoints = as.vector(unique(latent_adata$obs$timepoint))
timepoints = sort(timepoints)
plot_OT_heatmap(timepoints, 'spaDOT', file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_noremoval_results'), tp_color_palette=tp_color_palette)


# load SVGP results
# SVGP_pred_labels = read.csv(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_noremoval_results', 'SVGPVAE_predictions.csv'), header=T, row.names=1)
SVGP_pred_labels = read.csv(file.path(result_dir, 'explore_SpaDOT', 'ATAC-Visium_noremoval_results', 'SVGPVAE_predictions.csv'), header=T, row.names=1)
latent_adata$obs$SVGP_pred_labels = SVGP_pred_labels[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$SVGP_pred_labels = as.factor(latent_adata$obs$SVGP_pred_labels)
plot_pred_label(latent_adata$obs, color_col='SVGP_pred_labels', color_legend_title='spaDOT SVGP Predicted Labels', file_name=file.path(spatial_layout_dir, 'spaDOT_SVGP_prediction.png'))
# load GAT results
# GAT_pred_labels = read.csv(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_noremoval_results', 'GATVAE_predictions.csv'), header=T, row.names=1)
GAT_pred_labels = read.csv(file.path(result_dir, 'explore_SpaDOT', 'ATAC-Visium_noremoval_results', 'GATVAE_predictions.csv'), header=T, row.names=1)
latent_adata$obs$GAT_pred_labels = GAT_pred_labels[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$GAT_pred_labels = as.factor(latent_adata$obs$GAT_pred_labels)
plot_pred_label(latent_adata$obs, color_col='GAT_pred_labels', color_legend_title='spaDOT GAT Predicted Labels', file_name=file.path(spatial_layout_dir, 'spaDOT_GAT_prediction.png'))

# calculate LISI for SVGP / GAT
library(lisi)
for (tp in unique(latent_adata$obs$timepoint)) {
    tp_latent_adata = latent_adata[latent_adata$obs$timepoint == tp, ]
    locations = data.frame(x=tp_latent_adata$obs$pixel_x, y=tp_latent_adata$obs$pixel_y)
    lisi = compute_lisi(locations, tp_latent_adata$obs, c('SVGP_pred_labels', 'GAT_pred_labels'))
    median_lisi = median(lisi[, 1])
    cat(tp, median(lisi[, 1]), median(lisi[, 2]), '\n')
}


# --- load all methods' results
# spatialPCA
spatialPCA_latent = read_h5ad(file.path(result_dir, 'explore_spatialPCA', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
spatialPCA_latent$obs_names = sub("^t[0-9]_", "", spatialPCA_latent$obs_names)
latent_adata$obs$spatialPCA_pred_labels = spatialPCA_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='spatialPCA_pred_labels', color_legend_title='spatialPCA Predicted Labels', file_name=file.path(spatial_layout_dir, 'spatialPCA_prediction.png'))

# GraphST
GraphST_latent = read_h5ad(file.path(result_dir, 'explore_GraphST', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$GraphST_pred_labels = GraphST_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='GraphST_pred_labels', color_legend_title='GraphST Predicted Labels', file_name=file.path(spatial_layout_dir, 'GraphST_prediction.png'))

# SEDR
SEDR_latent = read_h5ad(file.path(result_dir, 'explore_SEDR', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$SEDR_pred_labels = SEDR_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='SEDR_pred_labels', color_legend_title='SEDR Predicted Labels', file_name=file.path(spatial_layout_dir, 'SEDR_prediction.png'))

# spaGCN
spaGCN_latent = read_h5ad(file.path(result_dir, 'explore_spaGCN', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$spaGCN_pred_labels = spaGCN_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='spaGCN_pred_labels', color_legend_title='spaGCN Predicted Labels', file_name=file.path(spatial_layout_dir, 'spaGCN_prediction.png'))

# spaVAE
spaVAE_latent = read_h5ad(file.path(result_dir, 'explore_spaVAE', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$spaVAE_pred_labels = spaVAE_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='spaVAE_pred_labels', color_legend_title='spaVAE Predicted Labels', file_name=file.path(spatial_layout_dir, 'spaVAE_prediction.png'))

# STAligner
STAligner_latent = read_h5ad(file.path(result_dir, 'explore_STAligner', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$STAligner_pred_labels = STAligner_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='STAligner_pred_labels', color_legend_title='STAligner Predicted Labels', file_name=file.path(spatial_layout_dir, 'STAligner_prediction.png'))

# PRECAST
PRECAST_latent = read_h5ad(file.path(result_dir, 'explore_PRECAST', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$PRECAST_pred_labels = PRECAST_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='PRECAST_pred_labels', color_legend_title='PRECAST Predicted Labels', file_name=file.path(spatial_layout_dir, 'PRECAST_prediction.png'))

# spatiAlign
spatiAlign_latent = read_h5ad(file.path(result_dir, 'explore_spatiAlign', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$spatiAlign_pred_labels = spatiAlign_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='spatiAlign_pred_labels', color_legend_title='spatiAlign Predicted Labels', file_name=file.path(spatial_layout_dir, 'spatiAlign_prediction.png'))

# SPIRAL
SPIRAL_latent = read_h5ad(file.path(result_dir, 'explore_SPIRAL', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$SPIRAL_pred_labels = SPIRAL_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='SPIRAL_pred_labels', color_legend_title='SPIRAL Predicted Labels', file_name=file.path(spatial_layout_dir, 'SPIRAL_prediction.png'))

# STADIA
STADIA_latent = read_h5ad(file.path(result_dir, 'explore_STADIA', 'ATAC-Visium_noremoval_results', 'prediction_adata.h5ad'))
latent_adata$obs$STADIA_pred_labels = STADIA_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata$obs, color_col='STADIA_pred_labels', color_legend_title='STADIA Predicted Labels', file_name=file.path(spatial_layout_dir, 'STADIA_prediction.png'))


# --- plot metrics and compute LISI
domain_metrics_df = read.csv(file.path(dirname(spatial_layout_dir), "final_ATAC-Visium_noremoval_domain_metrics.csv"), header=T)
domain_metrics_df$X = NULL
methods = unique(domain_metrics_df$method)
# method_color_palette = c("#F21A00", "#F3DF6C", "#00A08A", "#1b1811ff", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#446455")
method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae", "#ABDDDE")
customized_method_colors = method_color_palette[1:length(unique(domain_metrics_df$method))]
names(customized_method_colors) = unique(domain_metrics_df$method)
library(lisi)
for (method in methods) {
    pred_adata = read_h5ad(file.path(result_dir, paste0('explore_', method), 'ATAC-Visium_noremoval_results','prediction_adata.h5ad'))
    for (tp in unique(pred_adata$obs$timepoint)) {
        tp_pred_adata = pred_adata[pred_adata$obs$timepoint == tp, ]
        locations = data.frame(x=tp_pred_adata$obs$pixel_x, y=tp_pred_adata$obs$pixel_y)
        lisi = compute_lisi(locations, tp_pred_adata$obs, c('pred_labels'))
        median_lisi = median(lisi[, 1])
        domain_metrics_df[nrow(domain_metrics_df)+1,] = c('median-LISI', tp, median_lisi, method)
    }
}
# domain_metrics_df$method = factor(domain_metrics_df$method, levels=unique(domain_metrics_df$method))
domain_metrics_df$method = factor(domain_metrics_df$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST', 'STAMP'))
domain_metrics_df$value = as.numeric(domain_metrics_df$value)
for (metric in unique(domain_metrics_df$index)) {
    metric_df = domain_metrics_df[domain_metrics_df$index == metric, ]
    if (metric == 'median-LISI') {
        avg_metric_df = metric_df %>% group_by(method) %>% summarize(value=mean(as.numeric(value)), .groups="drop") %>% as.data.frame()
        avg_metric_df$index = metric
        avg_metric_df$timepoint = 'avg'
        avg_metric_df = avg_metric_df[, colnames(domain_metrics_df)]
        domain_metrics_df = rbind(domain_metrics_df, avg_metric_df)
    } else {
        avg_metric_df = metric_df[metric_df$timepoint == 'avg', ]
    }
    timepoint_metric_df = metric_df[metric_df$timepoint != 'avg', ]
    p = ggplot() +
        geom_hline(yintercept = avg_metric_df$value[avg_metric_df$method == "SpaDOT"],
           linetype = "dashed",
           color = customized_method_colors["SpaDOT"],
           size = 1) + 
        geom_point(data=timepoint_metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.3) +
        geom_point(data=avg_metric_df, aes(x=method, y=value, color=method, group=method), size=6, alpha=0.8) +
        scale_color_manual(values=customized_method_colors) + 
        theme_pubr() +
        theme(axis.text.x = element_blank()) + xlab(NULL) + ylab(metric)
    p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
    legend = get_legend(p)
    ggsave(file.path(dirname(spatial_layout_dir), 'legend.pdf'), plot=legend, width=6, height=2)
    p = p + theme(legend.position="none")
    ggsave(file.path(dirname(spatial_layout_dir), paste0(metric, ".pdf")), plot=p, width=4, height=3.5)
}
write.csv(domain_metrics_df, file.path(dirname(spatial_layout_dir), "SupplementaryTable4_metrics.csv"), row.names=F, quote=F)

q()

## derive cluster specific genes
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

# --- obtain marker gene overlap for each cluster called in each method
# cluster_marker_genes(Seurat_obj, 'spaDOT', method_result_dir=file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'spatialPCA', method_result_dir=file.path(result_dir, 'explore_spatialPCA', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'GraphST', method_result_dir=file.path(result_dir, 'explore_GraphST', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'SEDR', method_result_dir=file.path(result_dir, 'explore_SEDR', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'spaGCN', method_result_dir=file.path(result_dir, 'explore_spaGCN', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'spaVAE', method_result_dir=file.path(result_dir, 'explore_spaVAE', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'STAligner', method_result_dir=file.path(result_dir, 'explore_STAligner', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'PRECAST', method_result_dir=file.path(result_dir, 'explore_PRECAST', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'spatiAlign', method_result_dir=file.path(result_dir, 'explore_spatiAlign', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'SPIRAL', method_result_dir=file.path(result_dir, 'explore_SPIRAL', 'ATAC-Visium_noremoval_results'))
# cluster_marker_genes(Seurat_obj, 'STADIA', method_result_dir=file.path(result_dir, 'explore_STADIA', 'ATAC-Visium_noremoval_results'))


timepoints = as.vector(unique(Seurat_obj@meta.data$timepoint))
timepoints = timepoints[order(timepoints)]
# --- calculate overlap and OT -> does not have a good performance
spaDOT_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'spaDOT', method_result_dir=file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_noremoval_results'), perform_GO=TRUE, tp_color_palette=tp_color_palette)


# spatialPCA_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'spatialPCA', method_result_dir=file.path(result_dir, 'explore_spatialPCA', 'ATAC-Visium_noremoval_results'))
# GraphST_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'GraphST', method_result_dir=file.path(result_dir, 'explore_GraphST', 'ATAC-Visium_noremoval_results'))
# SEDR_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'SEDR', method_result_dir=file.path(result_dir, 'explore_SEDR', 'ATAC-Visium_noremoval_results'))
# spaGCN_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'spaGCN', method_result_dir=file.path(result_dir, 'explore_spaGCN', 'ATAC-Visium_noremoval_results'))
# spaVAE_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'spaVAE', method_result_dir=file.path(result_dir, 'explore_spaVAE', 'ATAC-Visium_noremoval_results'))
# STAligner_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'STAligner', method_result_dir=file.path(result_dir, 'explore_STAligner', 'ATAC-Visium_noremoval_results'))
# PRECAST_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'PRECAST', method_result_dir=file.path(result_dir, 'explore_PRECAST', 'ATAC-Visium_noremoval_results'))
# spatiAlign_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'spatiAlign', method_result_dir=file.path(result_dir, 'explore_spatiAlign', 'ATAC-Visium_noremoval_results'))
# SPIRAL_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'SPIRAL', method_result_dir=file.path(result_dir, 'explore_SPIRAL', 'ATAC-Visium_noremoval_results'))
# STADIA_spearman_coor = spearman_corr_markeroverlap_vs_ot(timepoints, 'STADIA', method_result_dir=file.path(result_dir, 'explore_STADIA', 'ATAC-Visium_noremoval_results'))


# --- plot gene expression
# paper marker genes for calculating enrichment score
marker_genes = c('Sox1', '3110039M20Rik', 'Foxg1', 'Pou3f2', 'Insm1', # 5: Forebrain
            'Crmp1', 'Myt1l1', 'Celf31', 'Nacad', 'Lemd11', # 7: Midbrain
            'Myt1l2', 'Rbfox32', 'Ank31', 'Ncam11', 'Xkr61', # 2: Hindbrain
            'Gm10800', 'Gm10801', 'Gm11168', 'Gm10722', 'Chl1', # 9: Periventricular
            'Gm108001', 'Gm108011', 'Chl11', 'Gm107221', 'Gm217381', # 3: Meningeal/PNS
            'Ino80', 'Adam19', 'Rasgrf22', 'Slc4a42', 'Kcnh82', # 4: Meningeal/PNS
            'Skint63', 'Vps13d', # 8: Internal
            'Foxp2', 'Nkx2-5', 'Tbx20', 'Actn2', # 1: Heart -> by guess, because original entry was entered as 33482
            'Rarg', 'Tcf7l1', 'Wwp2', 'Anxa2', 'Plec', # 0: Facial/bone
            'Hoxc5', 'Il17rd', 'Poc1a', 'Col26a11', 'Stra61', # 10: limb
            'Hbq1b', 'Zfpm1', 'Slc4a1', 'Nfe2', 'Eral1', 'Afp', 'Alb', 'Hnf4a', # 6: Liver,
            'Twist1', 'Snail', 'Vimentin' # mesenchyme
)

marker_genes = c('Vim', 'Ccnd2')
for (gene in marker_genes) {
    gene_idx = which(rownames(Seurat_obj) == gene)
    if (length(gene_idx) == 0) next
    gene_exprs = Seurat_obj@assays$RNA@layers$data[gene_idx, ]
    gene_exprs_df = Seurat_obj@meta.data
    gene_exprs_df[, gene] = gene_exprs
    ggplot(gene_exprs_df, aes(x=pixel_x, y=pixel_y, color=.data[[gene]])) + geom_point(size=0.5) +
        facet_wrap(~timepoint_original, nrow=1, scale='free') +
        theme_pubr() + 
        theme(axis.text =element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) + 
        scale_color_gradient(low = "grey", high = "red")
    ggsave(file.path(downstream_dir, paste0(gene, '_geneexprs.png')),height=4, width=10)    
}

# paper marker genes for calculating enrichment score
# marker_genes = c('Sox1'='CNS', 'Foxg1'='CNS', 'Notch1'='CNS', 'Fyn'='CNS',
#                  'Sox2'='CNS', 'Lhx2'='CNS', 'Emx1'='CNS', 'Neurog1'='CNS',
#                  'Cux2'='CNS', 'Neurod1'='CNS', 'Ascl1'='CNS', # CNS
#                 'Myh9'='Mesenchyme', 'Col1a1'='Mesenchyme', 'Smad3'='Mesenchyme', # Mesenchyme
#                 'Hba-a1'='Liver', 'Tal1'='Liver', 'Sptb'='Liver',
#                 'Gata4'='Liver', 'Gata6'='Liver', # Liver
#                 'Smad3'='BoneMuscle', 'Twist1'='BoneMuscle', 'Myog'='BoneMuscle') # bone and muscle

df = data.frame(matrix(nrow=0, ncol=5))
colnames(df) = c('timepoint', 'cluster', 'gene', 'region', 'enrichment')
for (tp in c('t0', 't1', 't2')) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    Idents(tp_Seurat_obj) = tp_Seurat_obj$pred_labels
    # Idents(tp_Seurat_obj) = tp_Seurat_obj$true_labels
    # --- paper marker gene enrichment
    for (cluster in levels(Idents(tp_Seurat_obj))) {
        cluster_cells_idx = Idents(tp_Seurat_obj) == cluster
        for (gene in names(marker_genes)) {
            avg_exprs = mean(tp_Seurat_obj[gene, cluster_cells_idx]$RNA@layers$counts)
            other_exprs = mean(tp_Seurat_obj[gene, !cluster_cells_idx]$RNA@layers$counts)
            fc_enrichment = (avg_exprs-other_exprs)/other_exprs
            df[nrow(df) + 1,] = c(tp, cluster, gene, unname(marker_genes[gene]), fc_enrichment)
        }
    }
    # --- DEGs among clusters
    # tp_markers = FindAllMarkers(tp_Seurat_obj)
    # write.csv(tp_markers, file.path(downstream_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')), quote=F)
    tp_markers = read.csv(file.path(downstream_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')),
                          header=TRUE, row.names=1)
    # sort markers by log2FC per group
    tp_top_markers = tp_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% top_n(5, avg_log2FC)

    subset_tp_markers = tp_markers[tp_markers$avg_log2FC > 1 & tp_markers$p_val_adj < 0.05,]
    # write.csv(subset_tp_markers, file.path(downstream_dir, paste0(tp, '_Seurat_FindAllMarkers_log2FC_0_adjp_0.01_original_cluster.csv')), quote=F)
    # gene enrichment analysis
    DEG_list = list()
    for (cluster in levels(Idents(tp_Seurat_obj))) {
        cluster_tp_markers = subset_tp_markers[subset_tp_markers$cluster == cluster, ]
        # Perform enrichment
        # enrich_kegg <- enrichKEGG(gene = cluster_tp_markers$gene, 
        #                         organism = 'mmu', 
        #                         pvalueCutoff = 0.05)
        enrich_go <- enrichGO(gene = cluster_tp_markers$gene, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
        enrich_go
        g = dotplot(enrich_go, showCategory=10)

        ggsave(file.path(downstream_dir, paste0(tp, '_cluster', cluster, '_enrichment_clusterprofiler.png')), width=8, height=6)
        # DEG_list[[as.character(cluster)]] = cluster_tp_markers$gene
        # DEG_ranks = sign(cluster_tp_markers$avg_log2FC) * -log10(cluster_tp_markers$p_val_adj)
        # names(DEG_ranks) = cluster_tp_markers$gene
        # inf_idx = which(is.infinite(DEG_ranks))
        # if (length(inf_idx) > 0) {
        #     DEG_ranks[inf_idx] = max(DEG_ranks[!is.infinite(DEG_ranks)])
        # }
        # DEG_ranks = DEG_ranks[order(-DEG_ranks)]
        # fgsea_res = fgsea(pathways = msigdbr_list, stats = DEG_ranks, minSize=15, maxSize=500)
        # # collapse pathway
        # #collapsedPathways = collapsePathways(fgsea_res[order(pval)][padj < 0.01],
        # #        msigdbr_list, DEG_ranks)
        # #mainPathways_res = fgsea_res[pathway %in% collapsedPathways$mainPathways][order(padj), ]
        # fgsea_res = fgsea_res[order(padj), ]
        # g_data = fgsea_res[1:10]
        # #g_data = mainPathways_res[1:10]
        # g_data$pathway = gsub("GOBP_", "", g_data$pathway)
        # g_data$pathway = gsub("_", " ", g_data$pathway)
        # g_data$pathway = factor(g_data$pathway, levels=rev(g_data$pathway))
        # ggplot(g_data, aes(y=pathway, x=-log10(padj), fill=-log10(padj))) +
        #     geom_bar(stat = "identity") +
        #     theme_pubr()
        # ggsave(file.path(downstream_dir, paste0(tp, '_cluster', cluster, '_fgsea_enrichment_original_cluster.png')), width=10, height=5)
    }
}
write.csv(df, file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster.csv'), quote=F)
# write.csv(df, file.path(downstream_dir, 'cluster_gene_enrichment_true_cluster.csv'), quote=F)
df$enrichment = as.numeric(df$enrichment)
df$gene_region = paste0(df$gene, '-', df$region)
df$gene_region = as.factor(df$gene_region)
df$gene_region = factor(df$gene_region,
                        levels=paste0(names(marker_genes),'-',unname(marker_genes)))
ggplot(data=df, aes(x=gene_region, y=enrichment, group=cluster)) +
    geom_line(aes(color=cluster)) + geom_point() +
    facet_wrap(~timepoint, ncol=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster.png'))
# ggsave(file.path(downstream_dir, 'cluster_gene_enrichment_true_cluster.png'))


# --- regenerate the ggplot based on enrichment and focus on marker genes provided by the study
# marker_genes = c('Sox1'='CNS', 'Foxg1'='CNS', 'Notch1'='CNS',
# 'Myh9'='Mesenchyme', 'Col1a1'='Mesenchyme', 'Smad3'='Mesenchyme',
# 'Hba-a1'='Liver', 'Tal1'='Liver', 'Sptb'='Liver')
# --- provided by chatGPT
# marker_genes = c('Tubb3'='CNS', 'Map2'='CNS', 'Sox2'='CNS', 'Sox1'='CNS', 'Nestin'='CNS', 'Olig2'='CNS', 'Gfap'='CNS', 'Mbp'='CNS', 'Snap25'='CNS', 'Syn1'='CNS', 'Syp'='CNS', 'S100b'='PNS', 'Mpz'='PNS', 'Plp1'='PNS', 'Ntrk1'='PNS', 'Ntrk2'='PNS', 'Ntrk3'='PNS', 'Isl1'='PNS', 'Phox2b'='PNS')

# --- provided by https://www.sciencedirect.com/science/article/pii/S009286741830789X?via%3Dihub#bib47
# PNS neurons, Figure S4, Usoskin et al.
PNS_neurons_marker_genes = c('Etv1', 'Gal', 'Cartpt', 'Nos1', 'Cox8b', 'Col12a1', 'Ly6e', 'Tac1', 'Trps1', 'Scgn', 'Ucn3', 'Nmu', 'Nog', 'Htr3a', 'Brip1', 'Ccne1', 'Hcrtr1', 'Cyp2j12', 'Pthlh', 'Ltk', 'Ctla2a', 'Vip', 'Slc18a3', 'Gda', 'Cgnl1', 'Sertm1', 'Itk', 'Dcdc2a', 'Prrx2', 'Angpt4', 'Pcdh19', 'Ntm', 'Cadps2', 'Pvalb', 'Gprc5b', 'Ptgfr', 'Cd34', 'Pkib', 'Lpar3', 'Otoa', 'Car8', 'Nppb')
names(PNS_neurons_marker_genes) = rep('PNS', length(PNS_neurons_marker_genes))
# Figure 6 and supplementary table 
CNS_neurons_marker_genes = c(
'Cplx3', 'Nxph3', 'Hs3st4', 'Ipcef1', 'Mpped1',
'Nr4a2', 'Col24a1', 'Oprk1', 'B3gat2', 'Kcnv1',
'Dkkl1', 'Galntl6', 'Pkib', 'Satb2', 'Nrep',
'Scube1', 'Tekt5', 'Fam19a2', 'Hrk', 'Plcxd2',
'Lypd1', 'Prss23', 'Kcng1', 'C1ql3', 'Myh7',
'Six3', 'Gm2115', 'Esrrg', 'Gad1', 'Ubash3b',
'Slc6a3', 'Chrna6', 'Chrnb3', 'Ddc', 'Aldh1a1',
'Pax2', 'Nyap2', 'Hoxb6', 'Tmem255a', 'Unc13c',
'Npnt', 'Gpc3', 'Tmem132c', 'Kcnd3', 'Tpd52l1',
'Syt2', 'Oscp1', 'Vamp1', 'Cd83', 'Scn1a'
)
CNS_neurons_marker_genes = unique(CNS_neurons_marker_genes)
names(CNS_neurons_marker_genes) = rep('CNS', length(CNS_neurons_marker_genes))
true_df = generate_marker_enrichment(Seurat_obj, label_col='true_labels', marker_genes=c(PNS_neurons_marker_genes, CNS_neurons_marker_genes))
write.csv(true_df, file.path(downstream_dir, 'cluster_gene_enrichment_true_cluster_CNSPNS.csv'), quote=F)
pred_df = generate_marker_enrichment(Seurat_obj, label_col='pred_labels', marker_genes=c(PNS_neurons_marker_genes, CNS_neurons_marker_genes))
write.csv(pred_df, file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster_CNSPNS.csv'), quote=F)

true_df = read.csv(file.path(downstream_dir, 'cluster_gene_enrichment_true_cluster_CNSPNS.csv'), header=T, row.names=1)
true_df = na.omit(true_df)
true_df$enrichment = as.numeric(true_df$enrichment)
true_df = true_df[abs(true_df$enrichment) > 1,]
true_df = true_df[!is.infinite(true_df$enrichment),]
# --- remove genes highly enriched in Heart, Cox8b, Nppb, Myh7, Itk
true_df = true_df[!true_df$gene %in% c('Cox8b', 'Nppb', 'Myh7', 'Itk'),]
true_kept_genes = c()
for (gene in unique(true_df$gene)) {
    true_df_gene = true_df[true_df$gene == gene,]   
    if (length(unique(true_df_gene$timepoint)) == 3) {
        true_kept_genes = c(true_kept_genes, gene)
    }
}

pred_df = read.csv(file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster_CNSPNS.csv'), header=T, row.names=1)
pred_df = pred_df %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "E12.5",
        timepoint == 't1' ~ "E13.5",
        timepoint == 't2' ~ "E15.5"
    ))
pred_df$pred_cluster = paste0(pred_df$timepoint_original, '_', pred_df$cluster)
pred_df = na.omit(pred_df)
pred_df$enrichment = as.numeric(pred_df$enrichment)
pred_df = pred_df[abs(pred_df$enrichment) > 1,]
pred_df = pred_df[!is.infinite(pred_df$enrichment),]
pred_kept_genes = c()
for (gene in unique(pred_df$gene)) {
    pred_df_gene = pred_df[pred_df$gene == gene,]   
    if (length(unique(pred_df_gene$timepoint)) == 3) {
        pred_kept_genes = c(pred_kept_genes, gene)
    }
}

common_marker_genes = intersect(true_kept_genes, pred_kept_genes)
true_df = true_df[true_df$gene %in% common_marker_genes,]
pred_df = pred_df[pred_df$gene %in% common_marker_genes,]

for (tp in unique(true_df$timepoint)) {
    tp_true_df = true_df[true_df$timepoint == tp, ]
    ggplot(data=tp_true_df, aes(x=gene, y=enrichment, group=cluster)) +
    geom_point(aes(color=cluster)) +
    facet_wrap(~region, ncol=3, scales='free') +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust=1))
    ggsave(file.path(downstream_dir, paste0(tp,'_cluster_gene_enrichment_true_cluster.png')), height=3, width=10)
}

for (tp in unique(pred_df$timepoint)) {
    tp_pred_df = pred_df[pred_df$timepoint == tp, ]
    tp_pred_df$cluster = as.factor(tp_pred_df$cluster)
    ggplot(data=tp_pred_df, aes(x=gene, y=enrichment, group=pred_cluster)) +
    geom_point(aes(color=pred_cluster)) +
    scale_color_manual(values=tp_color_palette) +  
    facet_wrap(~region, ncol=3, scales='free') +
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1))
    ggsave(file.path(downstream_dir, paste0(tp,'_cluster_gene_enrichment_original_cluster.png')), height=3, width=10)
}


# --- liver, heart, mesenchyme
liver_marker_genes = c('Alb', 'Cdh5', 'Clec4f', 'Spp1', 'Reln', 'Lyve1') # https://www.nature.com/articles/s41467-021-27354-w
names(liver_marker_genes) = rep('liver', length(liver_marker_genes))
heart_marker_genes = c('Ttn', 'Sln', 'Myl2', 'Pecam1', 'Npr3', 'Fabp4', 'Wt1', 'Postn', 'Pdgfrb', 'C1qa', 'Hba-a1') # https://www.nature.com/articles/s41467-022-35691-7
names(heart_marker_genes) = rep('heart', length(heart_marker_genes))
mesenchyme_marker_genes = c('Col3a1',  'Pcolce',  'Cnn2',  'Serpinh1', 'Rhoc',  'Spint2', 'Epcam',  'Krt8',  'Cldn6', 'Krt18')
#'Maged2', 'Fbn2', 'Cdh11', 'Col1a2', 'Sparc') # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1416-2
names(mesenchyme_marker_genes) = rep('mesenchyme', length(mesenchyme_marker_genes))
true_df = generate_marker_enrichment(Seurat_obj, label_col='true_labels', marker_genes=c(liver_marker_genes, heart_marker_genes, mesenchyme_marker_genes))
write.csv(true_df, file.path(downstream_dir, 'cluster_gene_enrichment_true_cluster_liver_heart_mesenchyme.csv'), quote=F)
pred_df = generate_marker_enrichment(Seurat_obj, label_col='pred_labels', marker_genes=c(liver_marker_genes, heart_marker_genes, mesenchyme_marker_genes))
write.csv(pred_df, file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster_liver_heart_mesenchyme.csv'), quote=F)

true_df = read.csv(file.path(downstream_dir, 'cluster_gene_enrichment_true_cluster_liver_heart_mesenchyme.csv'), header=T, row.names=1)
true_df = na.omit(true_df)
true_df$enrichment = as.numeric(true_df$enrichment)
true_df = true_df[!is.infinite(true_df$enrichment),]
true_kept_genes = c()
for (gene in unique(true_df$gene)) {
    true_df_gene = true_df[true_df$gene == gene,]   
    if (length(unique(true_df_gene$timepoint)) == 3) {
        true_kept_genes = c(true_kept_genes, gene)
    }
}

pred_df = read.csv(file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster_liver_heart_mesenchyme.csv'), header=T, row.names=1)
pred_df = pred_df %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "E12.5",
        timepoint == 't1' ~ "E13.5",
        timepoint == 't2' ~ "E15.5"
    ))
pred_df$pred_cluster = paste0(pred_df$timepoint_original, '_', pred_df$cluster)
pred_df = na.omit(pred_df)
pred_df$enrichment = as.numeric(pred_df$enrichment)
pred_df = pred_df[!is.infinite(pred_df$enrichment),]
pred_kept_genes = c()
for (gene in unique(pred_df$gene)) {
    pred_df_gene = pred_df[pred_df$gene == gene,]   
    if (length(unique(pred_df_gene$timepoint)) == 3) {
        pred_kept_genes = c(pred_kept_genes, gene)
    }
}

common_marker_genes = intersect(true_kept_genes, pred_kept_genes)
true_df = true_df[true_df$gene %in% common_marker_genes,]
pred_df = pred_df[pred_df$gene %in% common_marker_genes,]

for (tp in unique(true_df$timepoint)) {
    tp_true_df = true_df[true_df$timepoint == tp, ]
    ggplot(data=tp_true_df, aes(x=gene, y=enrichment, group=cluster)) +
    geom_point(aes(color=cluster)) +
    facet_wrap(~region, ncol=3, scales='free') +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 90, hjust=1))
    ggsave(file.path(downstream_dir, paste0(tp,'_cluster_gene_enrichment_true_cluster_heart_liver_mesenchyme.png')), height=3, width=10)
}

for (tp in unique(pred_df$timepoint)) {
    tp_pred_df = pred_df[pred_df$timepoint == tp, ]
    tp_pred_df$cluster = as.factor(tp_pred_df$cluster)
    ggplot(data=tp_pred_df, aes(x=gene, y=enrichment, group=pred_cluster)) +
    geom_point(aes(color=pred_cluster)) +
    scale_color_manual(values=tp_color_palette) +  
    facet_wrap(~region, ncol=3, scales='free') +
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1))
    ggsave(file.path(downstream_dir, paste0(tp,'_cluster_gene_enrichment_original_cluster_heart_liver_mesenchyme.png')), height=3, width=10)
}
