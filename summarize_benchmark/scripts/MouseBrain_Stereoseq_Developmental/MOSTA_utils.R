# compute_ASW <- function(location, cluster_label) {
#     require(cluster)
#     sil <- silhouette(cluster_label, dist(location))
#     mean(sil[, 3])  # Extract and return the average silhouette width
# }
plot_pred_label = function(latent_adata, color_col, file_name, tp_color_palette=NULL, legend_ncol=2, height=4, width=10) {
    require(cowplot)
    timepoints = unique(latent_adata$obs$timepoint_original)
    timepoints = sort(timepoints)
    plot_list = list()
    for (tp in timepoints) {
        g = ggplot(subset(latent_adata$obs, timepoint_original == tp), aes(x=pixel_x, y=pixel_y, color=!!sym(color_col))) + 
        geom_point(size=0.1) +
        guides(color = guide_legend(ncol=legend_ncol, title=NULL)) + 
        theme_void() + 
        ggtitle(tp) + 
        theme(legend.position = "bottom", plot.title=element_text(size=11, hjust=0.5))
        if (!is.null(tp_color_palette)) {
            g = g + scale_color_manual(values=tp_color_palette)
        }
        g = g + coord_fixed(ratio = 1)
        plot_list[[tp]] = g
    }
    plot_grid(plotlist=plot_list, rel_widths=c(0.5, 0.6, 0.7, 0.8, 0.9), ncol=length(timepoints), align = "hv", axis = "tblr")
    ggsave(file_name, height=height, width=width)

    # g = ggplot(df, aes(x=pixel_x, y=pixel_y, color=!!sym(color_col))) + geom_point(size=0.5) +
    #     facet_wrap(~timepoint_original, nrow=1, scale='free') +
    #     theme_pubr() + 
    #     theme(axis.text =element_blank(), axis.title=element_blank(),
    #     axis.ticks = element_blank(),
    #     axis.line=element_blank(),
    #     legend.title = element_text(size=15),
    #     legend.text = element_text(size=12)) + 
    #     scale_color_discrete(name=color_legend_title)
    # if (!is.null(tp_color_palette)) {
    #     g = g + scale_color_manual(values=tp_color_palette)
    # }
    # ggsave(file_name, height=height, width=width)
}

plot_OT_heatmap = function(timepoints, method, method_result_dir, tp_color_palette=NULL, prob_cutoff=0.2) {
    require(reshape2)
    if (method == 'spaDOT') {method = 'GATVAE-SVGP-VAE'}
    # calculate marker overlap
    timepoints_dict = c('t0'='E12.5', 't1'='E13.5', 't2'='E14.5', 't3'='E15.5', 't4'='E16.5')
    cluster_mapping = list(
    'E12.5'=c("1"="2", "2"="5", "3"="3", "4"="11", "5"="8", "6"="1", "7"="4"),
    'E13.5'=c("1"="3-2", "2"="11", "3"="1-1", "4"="3-1", "5"="1-2", "6"="2-2", "7"="2-1", "8"="4"),
    'E14.5'=c("1"="2-1", "2"="4", "3"="3-2", "4"="6", "5"="9", "6"="1-1", "7"="1-2", "8"="2-2", "9"="7", "10"="3-1"),
    'E15.5'=c("1"="1", "2"="11", "3"="2-2", "4"="3-1", "5"="2-1", "6"="10", "7"="6", "8"="4", "9"="3-2"),
    'E16.5'=c("1"="1-1", "2"="1-2", "3"="2-2", "4"="4-1", "5"="10", "6"="4-2", "7"="6", "8"="11", "9"="3-2", "10"="2-1")
    )
    order_vec = c("1", "1-1", "1-2", "2", "2-1", "2-2", "3", "3-1", "3-2", "4", "4-1", "4-2",
    "5", "6", "7", "8", "9", "10", "11")
    for (tp_i in 1:(length(timepoints)-1)) {
        cur_tp = timepoints[tp_i]
        next_tp = timepoints[tp_i+1]
        # obtain OT matrix
        OT_mat_row_norm = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '_row_norm.csv')), row.names=1)
        OT_mat_col_norm = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '_col_norm.csv')), row.names=1)
        # obtain minimum value for each entry
        OT_mat = data.frame(mapply(pmin, OT_mat_row_norm, OT_mat_col_norm))
        rownames(OT_mat) = rownames(OT_mat_row_norm)
        # OT_mat[OT_mat < prob_cutoff] = 0
        rownames(OT_mat) = as.character(sapply(strsplit(rownames(OT_mat), '_'), '[', 2))
        colnames(OT_mat) = as.character(sapply(strsplit(colnames(OT_mat), '_'), '[', 2))
        rownames(OT_mat) = unname(cluster_mapping[[timepoints_dict[cur_tp]]][rownames(OT_mat)])
        colnames(OT_mat) = unname(cluster_mapping[[timepoints_dict[next_tp]]][colnames(OT_mat)])
        sorted_rows_idx = match(order_vec, rownames(OT_mat))
        sorted_rows_idx = sorted_rows_idx[!is.na(sorted_rows_idx)]
        sorted_cols_idx = match(order_vec, colnames(OT_mat))
        sorted_cols_idx = sorted_cols_idx[!is.na(sorted_cols_idx)]
        OT_mat = OT_mat[sorted_rows_idx, sorted_cols_idx]
        x_axis_colors = unname(tp_color_palette[colnames(OT_mat)])
        y_axis_colors = unname(tp_color_palette[rownames(OT_mat)])

        OT_mat$ID = rownames(OT_mat)
        OT_mat_melt = melt(OT_mat, id.vars = "ID")
        OT_mat_melt$ID = factor( OT_mat_melt$ID, levels=order_vec)
        OT_mat_melt$variable = factor(OT_mat_melt$variable, levels=order_vec)
        colnames(OT_mat_melt) = c(cur_tp, next_tp, 'value')
        OT_mat_melt$value_edited = ifelse(OT_mat_melt$value < prob_cutoff, NA, OT_mat_melt$value)
        g = ggplot(OT_mat_melt, aes(!!sym(next_tp), !!sym(cur_tp), size=value, fill=value_edited)) +
            geom_point(shape=21) + 
            scale_size_continuous(range = c(1, 8), breaks=c(0.25, 0.5, 0.75), limits=c(0, 1), name = "Transition Probability") +
            scale_fill_gradientn(
                colors = c("white", "red"),
                values = c(0, 1),
                limits = c(0, 1),
                na.value = "grey",
                name = "Color scale"
            ) +
            labs(x = paste0(unname(timepoints_dict[next_tp]), " Domains"), 
                y = paste0(unname(timepoints_dict[cur_tp]), " Domains")) +
            theme_pubr() + 
            theme(legend.position="right",
                axis.text.x = element_text(color=x_axis_colors, angle=45, hjust=1),
                axis.text.y = element_text(color=y_axis_colors)
            )
        legend = get_legend(g)
        g = g + theme(legend.position="none")
        ggsave(file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_OT_heatmap_minimum.pdf')), width=4, height=4)
        ggsave(file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_OT_heatmap_legend_minimum.pdf')), cowplot::plot_grid(legend), width = 2, height = 4)
    }
}


cluster_marker_genes = function(Seurat_obj, method, method_result_dir) {
    require(Seurat)
    timepoints = as.vector(unique(Seurat_obj@meta.data$timepoint))
    timepoints = timepoints[order(timepoints)]

    if (method == 'spaDOT') {method_col='pred_labels'}
    else {method_col=paste0(method, '_pred_labels')}

    # remove ribosomal genes if any
    ribosomal_genes = grep("^RP[SL]", rownames(Seurat_obj), value=T)
    Seurat_obj =  Seurat_obj[!rownames(Seurat_obj) %in% ribosomal_genes, ]
    mt_genes = grep("^MT-", rownames(Seurat_obj), value=T)
    Seurat_obj =  Seurat_obj[!rownames(Seurat_obj) %in% mt_genes, ]
    for (tp in timepoints) {
        tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
        Idents(tp_Seurat_obj) = tp_Seurat_obj@meta.data[[method_col]]
        tp_markers = FindAllMarkers(tp_Seurat_obj)
        write.csv(tp_markers, file.path(method_result_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')), quote=F)    
    }
}

GO_analysis = function(timepoints, method, method_result_dir, tp_color_palette=NULL) {
    require(clusterProfiler)
    require(org.Gg.eg.db)
    require(reshape2)
    for (tp in timepoints) {
        tp_marker = read.csv(file.path(method_result_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')), row.names=1)
        # sub_tp_marker = tp_marker[tp_marker$avg_log2FC > 1 & tp_marker$p_val_adj < 0.05,]
        sub_tp_marker = tp_marker[tp_marker$p_val_adj < 0.05,] # relax on log2FC because of low FC
        for (cluster in unique(sub_tp_marker$cluster)) {
            sub_tp_cluster_marker = sub_tp_marker[sub_tp_marker$cluster==cluster, 'gene']
            write(sub_tp_cluster_marker, file=file.path(method_result_dir, paste0(tp, '_', cluster, '_marker_genes.txt')))
            tp_cluster_GO = enrichGO(gene = sub_tp_cluster_marker, 
                            keyType = 'SYMBOL',
                            OrgDb = org.Gg.eg.db, 
                            ont = "BP", # Biological process
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05)
            write.csv(as.data.frame(tp_cluster_GO), file.path(method_result_dir, paste0(tp, '_', cluster, '_enrichment_clusterprofiler.csv')), quote=F)
            g = dotplot(tp_cluster_GO, showCategory=10)
            ggsave(file.path(method_result_dir, paste0(tp, '_', cluster, '_enrichment_clusterprofiler.png')), width=8, height=6)
            # tp_cluster_GO = enrichGO(gene = sub_tp_cluster_marker[1:min(100, length(sub_tp_cluster_marker))], 
            #                 keyType = 'SYMBOL',
            #                 OrgDb = org.Gg.eg.db, 
            #                 ont = "BP", # Biological process
            #                 pAdjustMethod = "BH", 
            #                 qvalueCutoff = 0.05)
            # write.csv(as.data.frame(tp_cluster_GO), file.path(method_result_dir, paste0(tp, '_', cluster, '_top100_enrichment_clusterprofiler.csv')), quote=F)
            # g = dotplot(tp_cluster_GO, showCategory=10)
            # ggsave(file.path(method_result_dir, paste0(tp, '_', cluster, '_top100_enrichment_clusterprofiler.png')), width=8, height=6)
        }
    }

}


generate_marker_enrichment = function(Seurat_obj, label_col, marker_genes) {
    timepoints = as.vector(unique(Seurat_obj$timepoint))
    timepoints = sort(timepoints)
    df = data.frame(matrix(nrow=0, ncol=5))
    colnames(df) = c('timepoint', 'cluster', 'gene', 'region', 'enrichment')
    for (tp in timepoints) {
        tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
        Idents(tp_Seurat_obj) = tp_Seurat_obj@meta.data[[label_col]]
        for (cluster in levels(Idents(tp_Seurat_obj))) {
            cluster_cells_idx = Idents(tp_Seurat_obj) == cluster
            for (gene_idx in 1:length(marker_genes)) {
                gene = marker_genes[gene_idx]
                if (!unname(gene) %in% rownames(tp_Seurat_obj)) {
                    print(paste0(gene, ' not in Seurat object'))
                    next
                }
                avg_exprs = mean(tp_Seurat_obj[unname(gene), cluster_cells_idx]$RNA@layers$counts)
                other_exprs = mean(tp_Seurat_obj[unname(gene), !cluster_cells_idx]$RNA@layers$counts)
                fc_enrichment = (avg_exprs-other_exprs)/other_exprs
                df[nrow(df) + 1,] = c(tp, cluster, unname(gene), names(gene), fc_enrichment)
            }
        }
    }
    return(df)
}

gene_similarity = function(Seurat_obj, method_result_dir, tp_color_palette=NULL) {
    options(future.globals.maxSize = 8000 * 1024^2)
    require(dplyr)
    require(ggplot2)
    require(reshape2)
    timepoints = as.vector(unique(Seurat_obj@meta.data$timepoint))
    timepoints = timepoints[order(timepoints)]
    tp_count_list = list()
    for (tp in timepoints) {
        tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
        # tp_Seurat_obj = SCTransform(tp_Seurat_obj, verbose = FALSE) ## according to axolotl brain Fig S6D
        # average genes
        avg_genes = rowMeans(tp_Seurat_obj@assays$RNA@layers$data)
        avg_genes_idx = which(avg_genes > 0.5)
        tp_count_list[[tp]] = tp_Seurat_obj@assays$RNA@layers$data[avg_genes_idx, ]

        tp_subset_count_t = as.data.frame(t(tp_count_list[[tp]]))
        tp_subset_count_t$spaDOT_pred_labels = tp_Seurat_obj@meta.data$'spaDOT_pred_labels'
        tp_group_avg_exprs = tp_subset_count_t %>%
        group_by(spaDOT_pred_labels) %>% summarise(across(everything(), mean))

        tp_group_avg_exprs_df = as.data.frame(t(tp_group_avg_exprs[-1]))
        colnames(tp_group_avg_exprs_df) = tp_group_avg_exprs$spaDOT_pred_labels
        tp_cor_mat = cor(tp_group_avg_exprs_df, method='spearman')
        library(pheatmap)
        pheatmap(tp_cor_mat, color= colorRampPalette(c("blue", "white", "red"))(100), main="Gene expression correlation heatmap", cluster_rows = T, cluster_cols = T,
        #scale="row", 
        filename=file.path(method_result_dir, paste0(tp, '_gene_correlation_heatmap.png')))
    }
}