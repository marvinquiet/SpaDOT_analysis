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
library(fgsea)
library(msigdbr)
C5_gene_sets = msigdbr(species = "human", category = "C5", subcategory = 'GO:BP') # subcategory=
msigdbr_list = split(x = C5_gene_sets$gene_symbol, f = C5_gene_sets$gs_name)
source('Axolotl_utils.R')

base_dir = "/media/sda/wenjinma"
result_dir = file.path(base_dir, "projects/spatialATAC-explore")
spatial_layour_dir = file.path(result_dir, "summarize_benchmark/Axolotl_figures/spatial_layout")
dir.create(spatial_layour_dir, showWarnings = FALSE, recursive = TRUE)
downstream_dir = file.path(result_dir, "summarize_benchmark/Axolotl_figures/Axolotl_downstream")
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)
latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'Axolotl_results', 'prediction_adata.h5ad'))


color_palette = c(
    '#292444',
    '#ECCBAE',
    '#EE7E81',
    '#971E24',
    '#F3DF6C',
    '#F98400',
    '#55597A',
    '#A6C5D0',
    '#567B87',
    '#74AAB8')
annotations = unique(latent_adata$obs$true_labels)
customized_color_palette = color_palette[1:length(annotations)]
names(customized_color_palette) = annotations

# --- draw ground truth
latent_adata$obs = latent_adata$obs %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "St.44",
        timepoint == 't1' ~ "St.54",
        timepoint == 't2' ~ "St.57"
    ))
plot_pred_label(latent_adata, color_col='true_labels', file_name=file.path(spatial_layour_dir, 'Ground_Truth.pdf'), tp_color_palette=customized_color_palette)

# --- customize color palette for spaDOT prediction
# set color palette for different timepoint and plot spaDOT prediction
cluster_mapping = list(
    'St.44'=c("0"="6", "1"="5", "2"="3", "3"="2", "4"="4", "5"="1"),
    'St.54'=c("0"="1", "1"="5", "2"="3", "3"="4-1", "4"="4-2", "5"="6", "6"="2"),
    'St.57'=c("0"="4-1", "1"="2", "2"="4-2-1", "3"="6", "4"="5", "5"="4-2-2", "6"="3", "7"="1")
)

latent_adata$obs = latent_adata$obs %>% mutate(
    spaDOT_pred_labels = case_when(
        timepoint_original == 'St.44' ~ cluster_mapping[['St.44']][as.character(pred_labels)],
        timepoint_original == 'St.54' ~ cluster_mapping[['St.54']][as.character(pred_labels)],
        timepoint_original == 'St.57' ~ cluster_mapping[['St.57']][as.character(pred_labels)]
    )
)

tp_color_palette=c(
    "1"="#EE7E81",
    "2"="#F3DF6C",
    "3"="#971E24",
    "4"="#292444",
    "4-1"="#A6C5D0",
    "4-2"="#55597A",
    "4-2-1"="#567B87",
    "4-2-2"="#74AAB8",
    "5"="#F98400",
    "6"="#ECCBAE"
)

# latent_adata$obs$spaDOT_pred_labels = paste(latent_adata$obs$timepoint_original, latent_adata$obs$pred_labels, sep='_')
# tp_color_palette = list(
#     "St.44"=c("0"="#ECCBAE", "1"="#F98400", "2"="#971E24", "3"="#F3DF6C", "4"="#292444", "5"="#EE7E81"),
#     "St.54"=c("0"="#EE7E81", "1"="#F98400", "2"="#971E24", "3"="#A6C5D0", "4"="#55597A", "5"="#ECCBAE", "6"="#F3DF6C"),
#     "St.57"=c("0"="#A6C5D0", "1"="#F3DF6C", "2"="#567B87", "3"="#ECCBAE", "4"="#F98400", "5"="#74AAB8", "6"="#971E24", "7"="#EE7E81")
# )
plot_pred_label(latent_adata, color_col='spaDOT_pred_labels', file_name=file.path(spatial_layour_dir, 'spaDOT_prediction.pdf'), tp_color_palette=tp_color_palette)

# --- load all methods's results
# spatialPCA
spatialPCA_latent = read_h5ad(file.path(result_dir, 'explore_spatialPCA', 'Axolotl_results', 'prediction_adata.h5ad'))
spatialPCA_latent$obs_names = paste0(spatialPCA_latent$obs_names, '_', spatialPCA_latent$obs$timepoint)
latent_adata$obs$spatialPCA_pred_labels = spatialPCA_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spatialPCA_pred_labels', file_name=file.path(spatial_layour_dir, 'spatialPCA_prediction.png'))

# GraphST
GraphST_latent = read_h5ad(file.path(result_dir, 'explore_GraphST', 'Axolotl_results', 'prediction_adata.h5ad'))
latent_adata$obs$GraphST_pred_labels = GraphST_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='GraphST_pred_labels', file_name=file.path(spatial_layour_dir, 'GraphST_prediction.png'))

# SEDR
SEDR_latent = read_h5ad(file.path(result_dir, 'explore_SEDR', 'Axolotl_results', 'prediction_adata.h5ad'))
latent_adata$obs$SEDR_pred_labels = SEDR_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='SEDR_pred_labels', file_name=file.path(spatial_layour_dir, 'SEDR_prediction.png'))

# spaGCN
spaGCN_latent = read_h5ad(file.path(result_dir, 'explore_spaGCN', 'Axolotl_results', 'prediction_adata.h5ad'))
latent_adata$obs$spaGCN_pred_labels = spaGCN_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spaGCN_pred_labels', file_name=file.path(spatial_layour_dir, 'spaGCN_prediction.png'))

# spaVAE
spaVAE_latent = read_h5ad(file.path(result_dir, 'explore_spaVAE', 'Axolotl_results', 'prediction_adata.h5ad'))
latent_adata$obs$spaVAE_pred_labels = spaVAE_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spaVAE_pred_labels', file_name=file.path(spatial_layour_dir, 'spaVAE_prediction.png'))

# STAligner
STAligner_latent = read_h5ad(file.path(result_dir, 'explore_STAligner', 'Axolotl_results', 'prediction_adata.h5ad'))
latent_adata$obs$STAligner_pred_labels = STAligner_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='STAligner_pred_labels', file_name=file.path(spatial_layour_dir, 'STAligner_prediction.png'))

# PRECAST
PRECAST_latent = read_h5ad(file.path(result_dir, 'explore_PRECAST', 'Axolotl_results', 'prediction_adata.h5ad'))
PRECAST_latent$obs_names = paste0(PRECAST_latent$obs_names, '_', PRECAST_latent$obs$timepoint)
latent_adata$obs$PRECAST_pred_labels = PRECAST_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='PRECAST_pred_labels', file_name=file.path(spatial_layour_dir, 'PRECAST_prediction.png'))

# spatiAlign
spatiAlign_latent = read_h5ad(file.path(result_dir, 'explore_spatiAlign', 'Axolotl_results', 'prediction_adata.h5ad'))
latent_adata$obs$spatiAlign_pred_labels = spatiAlign_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='spatiAlign_pred_labels', file_name=file.path(spatial_layour_dir, 'spatiAlign_prediction.png'))

# SPIRAL
SPIRAL_latent = read_h5ad(file.path(result_dir, 'explore_SPIRAL', 'Axolotl_results', 'prediction_adata.h5ad'))
latent_adata$obs$SPIRAL_pred_labels = SPIRAL_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='SPIRAL_pred_labels', file_name=file.path(spatial_layour_dir, 'SPIRAL_prediction.png'))

# STADIA
STADIA_latent = read_h5ad(file.path(result_dir, 'explore_STADIA', 'Axolotl_results', 'prediction_adata.h5ad'))
STADIA_latent$obs_names = sub("^[^.]*\\.", "", STADIA_latent$obs_names)
STADIA_latent$obs_names = paste0(STADIA_latent$obs_names, '_', STADIA_latent$obs$timepoint)
latent_adata$obs$STADIA_pred_labels = STADIA_latent$obs[latent_adata$obs_names, 'pred_labels']
plot_pred_label(latent_adata, color_col='STADIA_pred_labels', file_name=file.path(spatial_layour_dir, 'STADIA_prediction.png'))

# load SVGP results
SVGP_latent = read.csv(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'Axolotl_results', 'SVGPVAE_predictions.csv'), header=T, row.names=1)
latent_adata$obs$SVGP_pred_labels = SVGP_latent[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$SVGP_pred_labels = as.factor(latent_adata$obs$SVGP_pred_labels)
plot_pred_label(latent_adata, color_col='SVGP_pred_labels', file_name=file.path(spatial_layour_dir, 'spaDOT_SVGP_prediction.png'))

GAT_pred_labels = read.csv(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'Axolotl_results', 'GATVAE_predictions.csv'), header=T, row.names=1)
latent_adata$obs$GAT_pred_labels = GAT_pred_labels[latent_adata$obs_names, 'pred_labels']
latent_adata$obs$GAT_pred_labels = as.factor(latent_adata$obs$GAT_pred_labels)
plot_pred_label(latent_adata, color_col='GAT_pred_labels', file_name=file.path(spatial_layour_dir, 'spaDOT_GAT_prediction.png'))

# --- compute IoU for MP&DP&LP, MP&DP/LP, MP/DP/LP
t0_latent_adata = latent_adata[latent_adata$obs$timepoint == 't0', ]
t1_latent_adata = latent_adata[latent_adata$obs$timepoint == 't1', ]
t2_latent_adata = latent_adata[latent_adata$obs$timepoint == 't2', ]
calculate_iou = function(y_true, pred_labels) {
    y_true <- as.logical(y_true)
    # obtain the one with highest overlap
    unique_pred_labels = unique(pred_labels)
    pred_labels_overlap_df = data.frame(matrix(nrow=0, ncol=3))
    colnames(pred_labels_overlap_df) = c('pred_labels', 'overlap', 'IoU')
    for (pred_label in unique_pred_labels) {
        y_pred = as.logical(pred_labels == pred_label)
        intersection = sum(y_true & y_pred)
        union = sum(y_true | y_pred)
        if (union == 0) {
            IoU = NA
        } else {
            IoU = intersection / union
        }
        pred_labels_overlap_df[nrow(pred_labels_overlap_df)+1,] = c(pred_label, intersection/sum(y_pred), IoU)
    }
    return(pred_labels_overlap_df)
    # if (nrow(pred_labels_overlap_df) == 0) {
    #     return(list(pred_labels=NA, overlap=NA, IoU=NA))
    # }
    # else {
    #     pred_labels_overlap_df = pred_labels_overlap_df[order(pred_labels_overlap_df$overlap, decreasing=TRUE), ]
    #     return(as.vector(pred_labels_overlap_df[1, ]))
    # }
}
# --- for all annotations in true labels
methods = c('SpaDOT', 'SEDR', 'STAligner', 'SPIRAL', 'GraphST', 'spaGCN', 'spaVAE', 'PRECAST', 'STADIA', 'spatiAlign', 'spatialPCA')
method_IoU_df = data.frame(matrix(nrow=0, ncol=6))
colnames(method_IoU_df) = c('method', 'pred_labels', 'overlap', 'IoU', 'timepoint', 'annotation')
for (tp in unique(latent_adata$obs$timepoint)) {
    tp_latent_adata = latent_adata[latent_adata$obs$timepoint == tp, ]
    for (true_label in unique(tp_latent_adata$obs$true_labels)) {
        tp_latent_adata_true_onehot = ifelse(tp_latent_adata$obs$true_labels == true_label, 1, 0)
        for (method in methods) {
            if (method == 'SpaDOT') {
                method_IoU = calculate_iou(tp_latent_adata_true_onehot, tp_latent_adata$obs[['pred_labels']])
            } else {
                method_IoU = calculate_iou(tp_latent_adata_true_onehot, tp_latent_adata$obs[[paste0(method, '_pred_labels')]])
            }
            method_IoU$timepoint = tp
            method_IoU$annotation = true_label
            method_IoU$method = method
            method_IoU = method_IoU[, colnames(method_IoU_df)]
            method_IoU_df = rbind(method_IoU_df, method_IoU)
        }
    }
}
write.csv(method_IoU_df, file=file.path(downstream_dir, 'IoU_summary.csv'), row.names=F)

# --- try Hungarian algorithm
library(dplyr)
library(tidyr)
library(clue)
method_IoU_df = method_IoU_df %>% mutate(pred_labels=as.character(pred_labels))
method_IoU_df$overlap = as.numeric(method_IoU_df$overlap)
assign_labels = function(sub_df) {
    # Create overlap matrix
    overlap_mat <- sub_df %>%
        select(pred_labels, annotation, overlap) %>%
        pivot_wider(names_from = annotation, values_from = overlap, values_fill = 0) %>%
        tibble::column_to_rownames("pred_labels") %>%
        as.matrix()
    
    # Solve the assignment problem (maximize overlap → minimize -overlap)
    assignment <- solve_LSAP(1-overlap_mat)
    matched <- tibble(
        pred_labels = rownames(overlap_mat),
        assigned_annotation = colnames(overlap_mat)[assignment]
    )
    
    # Join the matched annotation back
    sub_df %>%
        left_join(matched, by = "pred_labels")
}
result = method_IoU_df %>%
    group_by(method, timepoint) %>%
    group_modify(~ assign_labels(.x)) %>%
    ungroup()

filtered_result = as.data.frame(result[result$annotation == result$assigned_annotation, ])
# filtered_result = filtered_result %>% ## check duplicates
#   group_by(timepoint, method, pred_labels) %>%
#   mutate(
#     pred_labels = ifelse(n() > 1, NA, pred_labels),
#     overlap = ifelse(n() > 1, NA, overlap),
#     IoU = ifelse(n() > 1, NA, IoU)
#   ) %>%
#   ungroup() %>% as.data.frame()

filtered_result$overlap = as.numeric(filtered_result$overlap)
filtered_result$IoU = as.numeric(filtered_result$IoU)
filtered_result$annotation = as.factor(filtered_result$annotation)
# method_color_palette = c("#F21A00", "#F3DF6C", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#446455")
method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae")
names(method_color_palette) = c('SpaDOT', 'spatialPCA', 'GraphST', 'SEDR', 'spaGCN', 'spaVAE', 'STAligner', 'PRECAST', 'spatiAlign', 'SPIRAL', 'STADIA')
filtered_result$timepoint_annotation = paste(filtered_result$timepoint, filtered_result$annotation, sep=' - ')
# subset to only MP, DP, LP related
selected_annotations = c('t0 - MP&DP&LP', 't1 - MP&DP', 't1 - LP', 't2 - MP', 't2 - DP', 't2 - LP')
sub_filtered_result = filtered_result[filtered_result$timepoint_annotation %in% selected_annotations, ]
sub_filtered_result$timepoint_annotation = factor(sub_filtered_result$timepoint_annotation, levels=selected_annotations)
sub_filtered_result$method = factor(sub_filtered_result$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST'))
ggplot(sub_filtered_result, aes(x=overlap, y=IoU, color=method)) +
    geom_point(size=2, aes(shape=timepoint)) +
    facet_wrap(~timepoint_annotation, nrow=2, scales='free') +
    scale_color_manual(values=method_color_palette) +
    theme_pubr() + theme(legend.position='none') + xlab('Overlap Precentage') + ylab('IoU')
ggsave(file.path(downstream_dir, 'IoU_vs_Percentage.pdf'), width=8, height=4)




t0_latent_adata$obs$MPDPLP = ifelse(t0_latent_adata$obs$true_labels == 'MP&DP&LP', 1, 0)
for (method in methods) {
    if (method == 'SpaDOT') {
        method_IoU = calculate_iou(t0_latent_adata$obs$MPDPLP, t0_latent_adata$obs[['pred_labels']])
    } else {
        method_IoU = calculate_iou(t0_latent_adata$obs$MPDPLP, t0_latent_adata$obs[[paste0(method, '_pred_labels')]])
    }
    method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't0', 'MP&DP&LP')
}

# --- obtain for t1
t1_latent_adata$obs$MPDP = ifelse(t1_latent_adata$obs$true_labels == 'MP&DP', 1, 0)
for (method in methods) {
    if (method == 'SpaDOT') {
        method_IoU = calculate_iou(t1_latent_adata$obs$MPDP, t1_latent_adata$obs[['pred_labels']])
    } else {
        method_IoU = calculate_iou(t1_latent_adata$obs$MPDP, t1_latent_adata$obs[[paste0(method, '_pred_labels')]])
    }
    method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't1', 'MP&DP')
}
t1_latent_adata$obs$LP = ifelse(t1_latent_adata$obs$true_labels == 'LP', 1, 0)
for (method in methods) {
    if (method == 'SpaDOT') {
        method_IoU = calculate_iou(t1_latent_adata$obs$LP, t1_latent_adata$obs[['pred_labels']])
    } else {
        method_IoU = calculate_iou(t1_latent_adata$obs$LP, t1_latent_adata$obs[[paste0(method, '_pred_labels')]])
    }
    method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't1', 'LP')
}

t2_latent_adata$obs$MP = ifelse(t2_latent_adata$obs$true_labels == 'MP', 1, 0)
for (method in methods) {
    if (method == 'SpaDOT') {
        method_IoU = calculate_iou(t2_latent_adata$obs$MP, t2_latent_adata$obs[['pred_labels']])
    } else {
        method_IoU = calculate_iou(t2_latent_adata$obs$MP, t2_latent_adata$obs[[paste0(method, '_pred_labels')]])
    }
    method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't2', 'MP')
}
t2_latent_adata$obs$DP = ifelse(t2_latent_adata$obs$true_labels == 'DP', 1, 0)
for (method in methods) {
    if (method == 'SpaDOT') {
        method_IoU = calculate_iou(t2_latent_adata$obs$DP, t2_latent_adata$obs[['pred_labels']])
    } else {
        method_IoU = calculate_iou(t2_latent_adata$obs$DP, t2_latent_adata$obs[[paste0(method, '_pred_labels')]])
    }
    method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't2', 'DP')
}
t2_latent_adata$obs$LP = ifelse(t2_latent_adata$obs$true_labels == 'LP', 1, 0)
for (method in methods) {
    if (method == 'SpaDOT') {
        method_IoU = calculate_iou(t2_latent_adata$obs$LP, t2_latent_adata$obs[['pred_labels']])
    } else {
        method_IoU = calculate_iou(t2_latent_adata$obs$LP, t2_latent_adata$obs[[paste0(method, '_pred_labels')]])
    }
    method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't2', 'LP')
}
# write.csv(method_IoU_df, file=file.path(downstream_dir, 'MPDPLP_IoU_summary.csv'), row.names=F)

# --- perform the same analysis for VZ
# method_IoU_df = data.frame(matrix(nrow=0, ncol=6))
# colnames(method_IoU_df) = c('method', 'pred_labels', 'overlap', 'IoU', 'timepoint', 'annotation')
# t0_latent_adata$obs$VZ = ifelse(t0_latent_adata$obs$true_labels == 'VZ', 1, 0)
# for (method in methods) {
#     if (method == 'SpaDOT') {
#         method_IoU = calculate_iou(t0_latent_adata$obs$VZ, t0_latent_adata$obs[['pred_labels']])
#     } else {
#         method_IoU = calculate_iou(t0_latent_adata$obs$VZ, t0_latent_adata$obs[[paste0(method, '_pred_labels')]])
#     }
#     method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't0', 'VZ')
# }
# t1_latent_adata$obs$VZ = ifelse(t1_latent_adata$obs$true_labels == 'VZ', 1, 0)
# for (method in methods) {
#     if (method == 'SpaDOT') {
#         method_IoU = calculate_iou(t1_latent_adata$obs$VZ, t1_latent_adata$obs[['pred_labels']])
#     } else {
#         method_IoU = calculate_iou(t1_latent_adata$obs$VZ, t1_latent_adata$obs[[paste0(method, '_pred_labels')]])
#     }
#     method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't1', 'VZ')
# }
# t2_latent_adata$obs$VZ = ifelse(t2_latent_adata$obs$true_labels == 'VZ', 1, 0)
# for (method in methods) {
#     if (method == 'SpaDOT') {
#         method_IoU = calculate_iou(t2_latent_adata$obs$VZ, t2_latent_adata$obs[['pred_labels']])
#     } else {
#         method_IoU = calculate_iou(t2_latent_adata$obs$VZ, t2_latent_adata$obs[[paste0(method, '_pred_labels')]])
#     }
#     method_IoU_df[nrow(method_IoU_df)+1,] = c(method, method_IoU$pred_labels, method_IoU$overlap, method_IoU$IoU, 't2', 'VZ')
# }
# write.csv(method_IoU_df, file=file.path(downstream_dir, 'VZ_IoU_summary.csv'), row.names=F)

# --- visualize MP DP
MPDPLP_method_IoU_df = read.csv(file.path(downstream_dir, 'MPDPLP_IoU_summary.csv'), header=T)
# --- remove those methods that identify the same pred_label across annotations within time point
MPDPLP_method_IoU_dedup_df = MPDPLP_method_IoU_df %>%
  group_by(timepoint, method, pred_labels) %>%
  mutate(
    pred_labels = ifelse(n() > 1, NA, pred_labels),
    overlap = ifelse(n() > 1, NA, overlap),
    IoU = ifelse(n() > 1, NA, IoU)
  ) %>%
  ungroup() %>% as.data.frame()

MPDPLP_method_IoU_dedup_df$overlap = as.numeric(MPDPLP_method_IoU_dedup_df$overlap)
MPDPLP_method_IoU_dedup_df$IoU = as.numeric(MPDPLP_method_IoU_dedup_df$IoU)
MPDPLP_method_IoU_dedup_df$annotation = as.factor(MPDPLP_method_IoU_dedup_df$annotation)
# method_color_palette = c("#F21A00", "#F3DF6C", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#446455")
method_color_palette = c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0", "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae")
names(method_color_palette) = c('SpaDOT', 'spatialPCA', 'GraphST', 'SEDR', 'spaGCN', 'spaVAE', 'STAligner', 'PRECAST', 'spatiAlign', 'SPIRAL', 'STADIA')
MPDPLP_method_IoU_dedup_df$timepoint_annotation = paste(MPDPLP_method_IoU_dedup_df$timepoint, MPDPLP_method_IoU_dedup_df$annotation, sep=' - ')
MPDPLP_method_IoU_dedup_df$timepoint_annotation = factor(MPDPLP_method_IoU_dedup_df$timepoint_annotation, levels=c('t0 - MP&DP&LP', 't1 - MP&DP', 't1 - LP', 't2 - MP', 't2 - DP', 't2 - LP'))
MPDPLP_method_IoU_dedup_df$method = factor(MPDPLP_method_IoU_dedup_df$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST'))
ggplot(MPDPLP_method_IoU_dedup_df, aes(x=overlap, y=IoU, color=method)) +
    geom_point(size=2, aes(shape=timepoint)) +
    facet_wrap(~timepoint_annotation, nrow=2, scales='free') +
    scale_color_manual(values=method_color_palette) +
    theme_pubr() + theme(legend.position='none') + xlab('Overlap Precentage') + ylab('IoU')
ggsave(file.path(downstream_dir, 'MPDPLP_IoU.pdf'), width=8, height=4)

# --- visualize VZ_IOU in two dimensions in a dotplot
# VZ_method_IoU_df = read.csv(file.path(downstream_dir, 'VZ_IoU_summary.csv'), header=T)
# VZ_method_IoU_df$overlap = as.numeric(VZ_method_IoU_df$overlap)
# VZ_method_IoU_df$IoU = as.numeric(VZ_method_IoU_df$IoU)
# method_color_palette = c("#F21A00", "#F3DF6C", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6", "#ECCBAE", "#046C9A", "#D69C4E", "#ABDDDE", "#446455")
# names(method_color_palette) = methods
# VZ_method_IoU_df$method = factor(VZ_method_IoU_df$method, levels=methods)
# ggplot(VZ_method_IoU_df, aes(x=overlap, y=IoU, color=method)) +
#     geom_point(size=2) +
#     facet_wrap(~timepoint, nrow=1, scale='free') +
#     scale_color_manual(values=method_color_palette) +
#     theme_pubr() + xlab('Overlap Precentage') + ylab('IoU') + 
#     ylim(0, 1)
# ggsave(file.path(downstream_dir, 'VZ_IoU.png'), width=8, height=4)



# --- compute lisi
domain_metrics_df = read.csv(file.path(dirname(spatial_layour_dir), "final_Axolotl_domain_metrics.csv"), header=T)
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
domain_metrics_df$method = factor(domain_metrics_df$method, levels=c('SpaDOT', 'spatiAlign', 'spatialPCA', 'spaVAE', 'SPIRAL', 'SEDR', 'STADIA', 'PRECAST', 'spaGCN', 'STAligner', 'GraphST'))
metric_df = domain_metrics_df[domain_metrics_df$index == 'median-LISI', ]
metric_df$value = as.numeric(metric_df$value)
# compute average
avg_metric_df = metric_df %>% group_by(method) %>% summarize(value=mean(value), .groups="drop") %>% as.data.frame()
avg_metric_df$timepoint = 'avg'
avg_metric_df$index = 'median-LISI'
avg_metric_df = avg_metric_df[, colnames(domain_metrics_df)]
domain_metrics_df = rbind(domain_metrics_df, avg_metric_df)
write.csv(domain_metrics_df, file=file.path(dirname(spatial_layour_dir), 'SupplementaryTable6_metric.csv'), row.names=F)


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
    p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
    legend = get_plot_component(p, "guide-box", return_all=TRUE)[[4]]
    ggsave(file.path(dirname(spatial_layour_dir), 'legend.pdf'), plot=legend, width=6, height=2)
    p = p + theme(legend.position="none")
    ggsave(file.path(dirname(spatial_layour_dir), paste0(metric, ".pdf")), plot=p, width=4, height=3.5)
}


# p = ggplot() +
#     geom_point(data=avg_metric_df, aes(x=method, y=value, color=method, group=method), size=6) +
#     scale_color_manual(values=customized_method_colors) + 
#     geom_point(data=metric_df, aes(x=method, y=value, group=timepoint, shape=timepoint),size=3, alpha=0.5) +
#     theme_pubr() +
#     theme(axis.text.x = element_blank()) + xlab(NULL) + ylab('median-LISI')
# p = p + guides(color = guide_legend(nrow = 3, title=NULL), shape=guide_legend(nrow=3, title=NULL))
# p = p + theme(legend.position="none")
# ggsave(file.path(dirname(spatial_layour_dir), "median-LISI.png"), plot=p, width=4, height=3.5)

# --- plot heatmap
timepoints = as.vector(unique(latent_adata$obs$timepoint))
timepoints = sort(timepoints)
plot_OT_heatmap(timepoints, 'spaDOT', file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'Axolotl_results'), tp_color_palette=tp_color_palette)

# --- explore region proportion changes
true_prop_summary = latent_adata$obs %>% group_by(timepoint_original, true_labels) %>% summarize(Count=n(), .groups="drop") %>% group_by(timepoint_original) %>% mutate(Percentage=Count/sum(Count)*100) %>% as.data.frame()
true_prop_summary = merge(true_prop_summary, as.data.frame(customized_color_palette), by.x='true_labels', by.y='row.names')
true_prop_summary$color_palette = true_prop_summary$customized_color_palette

pred_prop_summary = latent_adata$obs %>% group_by(timepoint_original, pred_labels) %>% summarize(Count=n(), .groups="drop") %>% group_by(timepoint_original) %>% mutate(Percentage=Count/sum(Count)*100) %>% as.data.frame()
pred_colors = c()
for (i in 1:nrow(pred_prop_summary)) {
    pred_colors = c(pred_colors, unname(tp_color_palette[[pred_prop_summary$timepoint_original[i]]][as.character(pred_prop_summary$pred_labels[i])]))
}
pred_prop_summary$color_palette = pred_colors

# merge true and pred proportion summary
merged_prop_summary = merge(true_prop_summary, pred_prop_summary, by=c('timepoint_original', 'color_palette'))
cor_data = merged_prop_summary %>% group_by(timepoint_original) %>% summarize(cor=cor(Percentage.x, Percentage.y), .groups="drop") %>% as.data.frame()

g = ggplot(merged_prop_summary, aes(x=Percentage.x, y=Percentage.y, color=true_labels)) + geom_point(size=3) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~timepoint_original, nrow=1) +
    geom_text(data=cor_data, aes(label=paste0("r = ", round(cor, 2)), x = 30, y = 10), size=6, inherit.aes = FALSE) + 
    theme_pubr() + 
    theme(
    legend.title = element_text(size=15),
    legend.text = element_text(size=12)) + 
    scale_color_manual(name='Region Annotation', values=customized_color_palette) + 
    labs(x='True Proportion (%)', y='Predicted Proportion (%)')
ggsave(file.path(downstream_dir, 'true_vs_pred_proportion.png'), height=4, width=9)


## derive cluster specific genes
data_dir = file.path('/media/sda/wenjinma/data/spatial_Artista')
adata = read_h5ad(file.path(data_dir, 'Development_addCP.h5ad'))
adata$obs_names = paste0(adata$obs_names, '_', adata$obs$timepoint)
## extract gene name
extract_hs_genes <- function(entries) {
    hs_genes = c()
    for (entry in entries) {
        if (grepl(('[hs]'), entry)) {
            hs_gene = str_extract(entry, "(?<=\\|)[A-Z0-9]+(?=\\[)")
            hs_genes = c(hs_genes, hs_gene)
        } else if (!grepl('[nr]', entry) && grepl('\\|', entry)) {
            hs_gene = str_extract(entry, "\\b[A-Z0-9]+\\b(?=\\s*\\|)")
            hs_genes = c(hs_genes, hs_gene)
        } else {
            hs_genes = c(hs_genes, NA)
        }
    }
    return(hs_genes)
}
hs_genes = extract_hs_genes(adata$var_names)
adata$var$hs_genes = hs_genes
common_cells = intersect(rownames(adata), rownames(latent_adata))
adata = adata[common_cells]
latent_adata = latent_adata[common_cells, ]
# --- create Seurat object for analysis
Seurat_obj = CreateSeuratObject(counts=t(adata$X), meta.data=latent_adata$obs)
gene_metadata = data.frame(GeneName=rownames(Seurat_obj), hs_genes=hs_genes, row.names=rownames(Seurat_obj))
Seurat_obj@assays$RNA[[]] = gene_metadata
Seurat_obj[['spaDOT_latent']] = CreateAssayObject(data = t(latent_adata$X))
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
# Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
# Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
# Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)

# --- derive cluster specific genes
# cluster_marker_genes(Seurat_obj, 'spaDOT', method_result_dir=file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'Axolotl_results'))

# --- perform GO
GO_analysis(timepoints, gene_metadata, 'spaDOT', file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'Axolotl_results'), tp_color_palette=tp_color_palette)

St57_markers = read.csv(file.path(method_result_dir, paste0('t2', '_Seurat_FindAllMarkers_original_cluster.csv')), row.names=1)

# --- perform DEGs between t2 2 and 5
Idents(Seurat_obj) = "spaDOT_pred_labels"
# St57_2_vs_5_DE_genes = FindMarkers(Seurat_obj, ident.1 = "St.57_2", ident.2 = "St.57_5") # here, we want to consider 5 as baseline
St57_2_vs_5_DE_genes = FindMarkers(Seurat_obj, ident.1 = "4-2-1", ident.2 = "4-2-2" )

# volcano plot
St57_2_vs_5_DE_genes$negLog10adjP = -log10(St57_2_vs_5_DE_genes$p_val_adj)
St57_2_vs_5_DE_genes$significant = ifelse(St57_2_vs_5_DE_genes$p_val_adj < 0.05 & abs(St57_2_vs_5_DE_genes$avg_log2FC) > 0.5, 'Significant', 'Not Significant')
St57_2_vs_5_DE_genes = St57_2_vs_5_DE_genes[order(St57_2_vs_5_DE_genes$p_val_adj), ]
labeled_genes = rownames(St57_2_vs_5_DE_genes)[which(abs(St57_2_vs_5_DE_genes$avg_log2FC) > 0.5)[1:20]]
St57_2_vs_5_DE_genes$label = ifelse(St57_2_vs_5_DE_genes$significant == "Significant", rownames(St57_2_vs_5_DE_genes), NA)
ggplot(St57_2_vs_5_DE_genes, aes(x = avg_log2FC, y = negLog10adjP)) +
    geom_point(aes(color = significant), size = 2) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(
        title = "Volcano Plot",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
    geom_text_repel(aes(label = label), 
                  box.padding = 0.35, 
                  point.padding = 0.3, 
                  max.overlaps = 10)  # Adjust `max.overlaps` as needed
ggsave(file.path(downstream_dir, 'St57_2_vs_5_DEGs_volcano.pdf'), height=5, width=6)

# obtain upregulated genes
St57_2_vs_5_DE_genes_merged = merge(St57_2_vs_5_DE_genes, gene_metadata, by='row.names', all.x=TRUE)
St57_2_vs_5_DE_genes_merged = St57_2_vs_5_DE_genes_merged[!is.na(St57_2_vs_5_DE_genes_merged$hs_genes), ]
# up-regulated genes
St57_2_vs_5_DE_up_genes = St57_2_vs_5_DE_genes_merged[St57_2_vs_5_DE_genes_merged$avg_log2FC > 0 &St57_2_vs_5_DE_genes_merged$p_val_adj < 0.05, 'hs_genes']
St57_2_vs_5_DE_up_genes_GO = enrichGO(St57_2_vs_5_DE_up_genes, keyType = 'SYMBOL',
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP", # Biological process
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)
g = dotplot(St57_2_vs_5_DE_up_genes_GO, showCategory=10)
ggsave(file.path(downstream_dir,'St57_2_vs_5_UP_enrichment_clusterprofiler.png'), width=8, height=6)
write(St57_2_vs_5_DE_up_genes, file=file.path(downstream_dir, 'St57_2_vs_5_upregulated_genes.txt'))

St_57_2_vs_5_DE_down_genes = St57_2_vs_5_DE_genes_merged[St57_2_vs_5_DE_genes_merged$avg_log2FC < 0 & St57_2_vs_5_DE_genes_merged$p_val_adj < 0.05, 'hs_genes']
St_57_2_vs_5_DE_down_genes_GO = enrichGO(St_57_2_vs_5_DE_down_genes, keyType = 'SYMBOL',
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP", # Biological process
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)
g = dotplot(St_57_2_vs_5_DE_down_genes_GO, showCategory=10)
ggsave(file.path(downstream_dir,'St57_2_vs_5_DOWN_enrichment_clusterprofiler.png'), width=8, height=6)
write(St_57_2_vs_5_DE_down_genes, file=file.path(downstream_dir, 'St57_2_vs_5_downregulated_genes.txt'))

common_GOs = intersect(St57_2_vs_5_DE_up_genes_GO$ID, St_57_2_vs_5_DE_down_genes_GO$ID)

# --- St57_2_vs_5 gene expression
St57_2_vs_5_dir = file.path(downstream_dir, 'St57_2_vs_5_DEGs')
dir.create(St57_2_vs_5_dir, showWarnings = FALSE, recursive = TRUE)
St57_2_vs_5_significant_genes = rownames(St57_2_vs_5_DE_genes[St57_2_vs_5_DE_genes$significant == "Significant", ])
St57_Seurat_obj = subset(Seurat_obj, ident=c("St.57_2", "St.57_5"))
for (gene in St57_2_vs_5_significant_genes) {
    gene_idx = which(rownames(St57_Seurat_obj) == gene)
    if (length(gene_idx) == 0) next
    gene_exprs = St57_Seurat_obj@assays$RNA@layers$data[gene_idx, ]
    gene_exprs_df = St57_Seurat_obj@meta.data
    gene_exprs_df[, gene] = gene_exprs
    ggplot(gene_exprs_df, aes(x=pixel_x, y=pixel_y, color=.data[[gene]])) + geom_point(size=0.5) +
        theme_pubr() + 
        theme(axis.text =element_blank(), axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12)) + 
        scale_color_gradient(low = "grey", high = "red")
    ggsave(file.path(St57_2_vs_5_dir, paste0(gene, '_geneexprs.png')),height=4, width=10)    
}

# # --- St57_MP_vs_DP gene expression
# St57_MP_vsDP_dir = file.path(downstream_dir, 'St57_MP_vs_DP_DEGs')
# dir.create(St57_MP_vsDP_dir, showWarnings = FALSE, recursive = TRUE)
# Seurat_obj$tp_true_labels = paste(Seurat_obj$timepoint_original, Seurat_obj$true_labels, sep='_')
# Idents(Seurat_obj) = "tp_true_labels"
# St57_MP_vs_DP_DE_genes = FindMarkers(Seurat_obj, ident.1 = "St.57_MP", ident.2 = "St.57_DP")
# # volcano plot
# St57_MP_vs_DP_DE_genes$negLog10adjP = -log10(St57_MP_vs_DP_DE_genes$p_val_adj)
# St57_MP_vs_DP_DE_genes$significant = ifelse(St57_MP_vs_DP_DE_genes$p_val_adj < 0.05 & abs(St57_MP_vs_DP_DE_genes$avg_log2FC) > 1, 'Significant', 'Not Significant')
# St57_MP_vs_DP_DE_genes$label = ifelse(St57_MP_vs_DP_DE_genes$significant == "Significant", rownames(St57_MP_vs_DP_DE_genes), NA)
# ggplot(St57_MP_vs_DP_DE_genes, aes(x = avg_log2FC, y = negLog10adjP)) +
#     geom_point(aes(color = significant), size = 2) +
#     scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
#     theme_minimal() +
#     labs(
#         title = "Volcano Plot",
#         x = "Log2 Fold Change",
#         y = "-Log10 Adjusted P-value"
#     ) +
#     theme(
#         plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank()
#     ) +
#     geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
#     geom_text_repel(aes(label = label), 
#                   box.padding = 0.35, 
#                   point.padding = 0.3, 
#                   max.overlaps = 10)  # Adjust `max.overlaps` as needed
# ggsave(file.path(downstream_dir, 'St57_MP_vs_DP_DEGs_volcano.png'), height=20, width=20)

# # obtain upregulated genes
# St57_MP_vs_DP_DE_genes_merged = merge(St57_MP_vs_DP_DE_genes, gene_metadata, by='row.names', all.x=TRUE)
# St57_MP_vs_DP_DE_genes_merged = St57_MP_vs_DP_DE_genes_merged[!is.na(St57_MP_vs_DP_DE_genes_merged$hs_genes), ]
# # up-regulated genes
# St57_MP_vs_DP_DE_up_genes = St57_MP_vs_DP_DE_genes_merged[St57_MP_vs_DP_DE_genes_merged$avg_log2FC > 0.1 & St57_MP_vs_DP_DE_genes_merged$p_val_adj < 0.05, 'hs_genes']
# St57_MP_vs_DP_DE_up_genes_GO = enrichGO(St57_MP_vs_DP_DE_up_genes, keyType = 'SYMBOL',
#                                 OrgDb = org.Hs.eg.db, 
#                                 ont = "BP", # Biological process
#                                 pAdjustMethod = "BH", 
#                                 qvalueCutoff = 0.05)
# g = dotplot(St57_MP_vs_DP_DE_up_genes_GO, showCategory=10)
# ggsave(file.path(downstream_dir,'St57_MP_vs_DP_UP_enrichment_clusterprofiler.png'), width=8, height=6)
# write(St57_MP_vs_DP_DE_up_genes, file=file.path(downstream_dir, 'St57_MP_vs_DP_upregulated_genes.txt'))

# St57_MP_vs_DP_DE_down_genes = St57_MP_vs_DP_DE_genes_merged[St57_MP_vs_DP_DE_genes_merged$avg_log2FC < -0.1 & St57_MP_vs_DP_DE_genes_merged$p_val_adj < 0.05, 'hs_genes']
# St57_MP_vs_DP_DE_down_genes_GO = enrichGO(St57_MP_vs_DP_DE_down_genes, keyType = 'SYMBOL',
#                                 OrgDb = org.Hs.eg.db, 
#                                 ont = "BP", # Biological process
#                                 pAdjustMethod = "BH", 
#                                 qvalueCutoff = 0.05)
# g = dotplot(St57_MP_vs_DP_DE_down_genes_GO, showCategory=10)
# ggsave(file.path(downstream_dir,'St57_MP_vs_DP_DOWN_enrichment_clusterprofiler.png'), width=8, height=6)
# write(St57_MP_vs_DP_DE_down_genes, file=file.path(downstream_dir, 'St57_MP_vs_DP_downregulated_genes.txt'))



# # --- St57_MP_vs_DP gene expression
# St57_MP_vs_DP_significant_genes = rownames(St57_MP_vs_DP_DE_genes[St57_MP_vs_DP_DE_genes$significant == "Significant", ])
# St57_Seurat_obj = subset(Seurat_obj, ident=c("St.57_MP", "St.57_DP"))
# for (gene in St57_MP_vs_DP_significant_genes) {
#     gene_idx = which(rownames(St57_Seurat_obj) == gene)
#     if (length(gene_idx) == 0) next
#     gene_exprs = St57_Seurat_obj@assays$RNA@layers$data[gene_idx, ]
#     gene_exprs_df = St57_Seurat_obj@meta.data
#     gene_exprs_df[, gene] = gene_exprs
#     ggplot(gene_exprs_df, aes(x=pixel_x, y=pixel_y, color=.data[[gene]])) + geom_point(size=0.5) +
#         theme_pubr() + 
#         theme(axis.text =element_blank(), axis.title=element_blank(),
#         axis.ticks = element_blank(),
#         axis.line=element_blank(),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=12)) + 
#         scale_color_gradient(low = "grey", high = "red")
#     ggsave(file.path(St57_MP_vsDP_dir, paste0(gene, '_geneexprs.png')),height=4, width=10)    
# }


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


# --- plot marker genes
# https://www.science.org/doi/10.1126/science.abp9262, Figure2B (Elmo1, Satb1, Rbfox3)
# https://elifesciences.org/articles/13998, Figure 1D (Er81 (ETV1 in human), Fezf2)
marker_genes = c("ELMO1 - AMEX60DD038068", "SATB1 - AMEX60DD022782", "RBFOX3 - AMEX60DD031407", "FEZF2 - AMEX60DD023808", "ETV1 - AMEX60DD022427", "SST - AMEX60DD003175", "TMIE - AMEX60DD020531", "C1QL1 - AMEX60DD009856", "CBLN2 - AMEX60DD039277", "SLC1A3 - AMEX60DD041912", "PTN[nr] - AMEX60DD007022", "GFAP - AMEX60DD009861", 'CALB2 - AMEX60DD015523', "VIM - AMEX60DD022108", "RHO - AMEX60DD023501", "BCL6 - AMEX60DD003177", "RBP1 - AMEX60DD002752")
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

marker_genes = c('SLC1A3 - AMEX60DD041912', 'GFAP - AMEX60DD009861', 'VIM - AMEX60DD022108', 'FABP7[nr] - AMEX60DD034100', 'SOX3[nr]-SOX2[hs] - AMEX60DD037395', 'WNT8B - AMEX60DD051636', 'SFRP1 - AMEX60DD002950')
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

# plot marker gene enrichment
VZ_markers = c('SLC1A3 - AMEX60DD041912', 'GFAP - AMEX60DD009861', 'VIM - AMEX60DD022108', 'FABP7[nr] - AMEX60DD034100', 'SOX3[nr]-SOX2[hs] - AMEX60DD037395', 'WNT8B - AMEX60DD051636', 'SFRP1 - AMEX60DD002950', 'FZD5 - AMEX60DD055073', 'SOX1 - AMEX60DD047894', 'LOC101342876[nr]-STMN2[hs] - AMEX60DD039970', 'DCX - AMEX60DD037111', 'NES - AMEX60DD016851')
names(VZ_markers) = rep('VZ', length(VZ_markers))
pred_df = generate_marker_enrichment(Seurat_obj, label_col='pred_labels', marker_genes=VZ_markers)
write.csv(pred_df, file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster_VZ.csv'), quote=F)

pred_df = read.csv(file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster_VZ.csv'), header=T, row.names=1)
pred_df = pred_df %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "St.44",
        timepoint == 't1' ~ "St.54",
        timepoint == 't2' ~ "St.57"
    ))
# pred_df$pred_cluster = paste0(pred_df$timepoint_original, '_', pred_df$cluster)
# --- mutate labels
pred_df = pred_df %>% mutate(
    pred_cluster = case_when(
        timepoint_original == 'St.44' ~ cluster_mapping[['St.44']][as.character(cluster)],
        timepoint_original == 'St.54' ~ cluster_mapping[['St.54']][as.character(cluster)],
        timepoint_original == 'St.57' ~ cluster_mapping[['St.57']][as.character(cluster)]
    )
)
pred_df = na.omit(pred_df)
pred_df$enrichment = as.numeric(pred_df$enrichment)
pred_df = pred_df[!is.infinite(pred_df$enrichment),]
pred_df$gene = sapply(strsplit(pred_df$gene, ' - '), `[`, 1)

pred_kept_genes = c()
for (gene in unique(pred_df$gene)) {
    pred_df_gene = pred_df[pred_df$gene == gene,]   
    if (length(unique(pred_df_gene$timepoint)) == 3) {
        pred_kept_genes = c(pred_kept_genes, gene)
    }
}
pred_df = pred_df[pred_df$gene %in% pred_kept_genes,]
pred_df = pred_df[!pred_df$gene %in% c('LOC101342876[nr]-STMN2[hs]', 'VIM', 'WNT8B', 'DCX'),]
for (tp in unique(pred_df$timepoint_original)) {
    tp_pred_df = pred_df[pred_df$timepoint_original == tp, ]
    tp_pred_df$pred_cluster = as.factor(tp_pred_df$pred_cluster)
    ggplot(data=tp_pred_df, aes(x=gene, y=enrichment, group=pred_cluster)) +
    geom_point(aes(color=pred_cluster)) +
    scale_color_manual(values=tp_color_palette) +  
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))
    ggsave(file.path(downstream_dir, paste0(tp,'_cluster_gene_enrichment_VZ.pdf')), height=4, width=4)
}


marker_genes = c("LOC115074257[nr]-NPY[hs] - AMEX60DD009927", "KCNS3 - AMEX60DD033135", "BHLHE22 - AMEX60DD039815", "C1QL1 - AMEX60DD009856", "TLL1 - AMEX60DD045228", "NEUROD2 - AMEX60DD010297", )
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

# --- plot gene similarity
require(dplyr)
require(ggplot2)
require(reshape2)
library(pheatmap)
timepoints = as.vector(unique(Seurat_obj@meta.data$timepoint))
timepoints = timepoints[order(timepoints)]
# --- compute only for MD&PD&LD -> true label
selected_idx = (Seurat_obj$timepoint == 't0' & Seurat_obj$true_labels %in% c('MP&DP&LP')) | (Seurat_obj$timepoint == 't1' & Seurat_obj$true_labels %in% c('MP&DP', 'LP')) | (Seurat_obj$timepoint == 't2' & Seurat_obj$true_labels %in% c('MP', 'DP', 'LP'))
sub_Seurat_obj = Seurat_obj[, selected_idx]
sub_Seurat_obj = SCTransform(sub_Seurat_obj, verbose = FALSE)
# average genes
avg_genes = rowMeans(sub_Seurat_obj@assays$SCT@data)
# remove those < 0.5
avg_genes = avg_genes[avg_genes > 0.5]
sub_count_t = as.data.frame(t(sub_Seurat_obj@assays$SCT@data[names(avg_genes), ]))
sub_count_t$true_labels = paste0(sub_Seurat_obj@meta.data[rownames(sub_count_t), 'timepoint_original'], '_', sub_Seurat_obj@meta.data[rownames(sub_count_t), 'true_labels'])
group_avg_exprs = sub_count_t %>%
    group_by(true_labels) %>% summarise(across(everything(), mean))
group_avg_exprs_df = as.data.frame(t(group_avg_exprs[-1]))
colnames(group_avg_exprs_df) = group_avg_exprs$true_labels
cor_mat = cor(group_avg_exprs_df, method='spearman')
pheatmap(cor_mat, color= colorRampPalette(c("blue", "white", "red"))(100), main="Gene expression correlation heatmap", cluster_rows = T, cluster_cols = T, filename=(file.path(method_result_dir, 'MDPDLP_truelabels_gene_correlation_heatmap.png')))


# -- compare for predicted labels
selected_idx = (Seurat_obj$timepoint == 't0' & Seurat_obj$pred_labels %in% c(4)) | (Seurat_obj$timepoint == 't1' & Seurat_obj$pred_labels %in% c(3,4)) | (Seurat_obj$timepoint == 't2' & Seurat_obj$pred_labels %in% c(0,2,5))
sub_Seurat_obj = Seurat_obj[, selected_idx]
sub_Seurat_obj = SCTransform(sub_Seurat_obj, verbose = FALSE)
# --- plot gene count
for (tp in unique(sub_Seurat_obj$timepoint_original)) {
    sub_tp_Seurat_obj = subset(sub_Seurat_obj, subset=timepoint_original==tp)
    ggplot(sub_tp_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=nCount_RNA)) + geom_point(size=0.5) +
            theme_pubr() + 
            theme(axis.text =element_blank(), axis.title=element_blank(),
            axis.ticks = element_blank(),
            axis.line=element_blank(),
            legend.title = element_text(size=15),
            legend.text = element_text(size=12)) + 
            scale_color_gradient(low = "grey", high = "red")
    ggsave(file.path(downstream_dir, paste0(tp, '_subset_nCount_SCT_geneexprs.png')),height=4, width=10)
    # --- histogram of each cluster within each time point
    ggplot(sub_tp_Seurat_obj@meta.data, aes(x=log(nCount_SCT), fill=spaDOT_pred_labels)) + geom_histogram(bins=30, alpha=0.5) + theme_pubr()
    ggsave(file.path(downstream_dir, paste0(tp, '_subset_nCount_SCT_histogram.png')),height=4, width=4)
}

# average genes
avg_genes = rowMeans(sub_Seurat_obj@assays$SCT@data)
# remove those < 0.5
avg_genes = avg_genes[avg_genes > 0.5]
sub_count_t = as.data.frame(t(sub_Seurat_obj@assays$SCT@data[names(avg_genes), ]))
# sub_count_t$pred_labels = paste0(sub_Seurat_obj@meta.data[rownames(sub_count_t), 'timepoint_original'], '_', sub_Seurat_obj@meta.data[rownames(sub_count_t), 'pred_labels'])
sub_count_t$pred_labels = sub_Seurat_obj$spaDOT_pred_labels
group_avg_exprs = sub_count_t %>%
    group_by(pred_labels) %>% summarise(across(everything(), mean))
group_avg_exprs_df = as.data.frame(t(group_avg_exprs[-1]))
colnames(group_avg_exprs_df) = group_avg_exprs$pred_labels
cor_mat = cor(group_avg_exprs_df, method='spearman')
pheatmap(cor_mat, color= colorRampPalette(c("blue", "white", "red"))(100), main="Gene expression correlation heatmap", cluster_rows = T, cluster_cols = T, filename=(file.path(method_result_dir, 'MDPDLP_predlabels_gene_correlation_heatmap.pdf')))



tp_count_list = list()
for (tp in timepoints) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    tp_Seurat_obj = SCTransform(tp_Seurat_obj, verbose = FALSE)
    # average genes
    avg_genes = rowMeans(tp_Seurat_obj@assays$SCT@data)
    # remove those < 0.5
    avg_genes = avg_genes[avg_genes > 0.5]
    tp_count_list[[tp]] = tp_Seurat_obj@assays$SCT@data[names(avg_genes), ]
    tp_subset_count_t = as.data.frame(t(tp_count_list[[tp]]))
    tp_subset_count_t$spaDOT_pred_labels = tp_Seurat_obj@meta.data[rownames(tp_subset_count_t), 'spaDOT_pred_labels']
    tp_group_avg_exprs = tp_subset_count_t %>%
    group_by(spaDOT_pred_labels) %>% summarise(across(everything(), mean))
    tp_group_avg_exprs_df = as.data.frame(t(tp_group_avg_exprs[-1]))
    colnames(tp_group_avg_exprs_df) = tp_group_avg_exprs$spaDOT_pred_labels
    tp_cor_mat = cor(tp_group_avg_exprs_df, method='spearman')
    pheatmap(tp_cor_mat, color= colorRampPalette(c("blue", "white", "red"))(100), main="Gene expression correlation heatmap", cluster_rows = T, cluster_cols = T,
    #scale="row", 
    filename=file.path(method_result_dir, paste0(tp, '_gene_correlation_heatmap.png')))
}
