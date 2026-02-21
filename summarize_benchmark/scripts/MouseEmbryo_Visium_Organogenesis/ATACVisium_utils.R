
plot_pred_label = function(df, color_col, color_legend_title, file_name, tp_color_palette=NULL, height=4, width=10) {
    g = ggplot(df, aes(x=pixel_x, y=pixel_y, color=!!sym(color_col))) + geom_point(size=0.5) +
        facet_wrap(~timepoint_original, nrow=1, scale='free') +
        theme_pubr(base_size = 14) + 
        #theme_bw()+
        theme(
            axis.text =element_blank(), axis.title=element_blank(),
            axis.ticks = element_blank(),
            axis.line=element_blank(),
            legend.title = element_text(size = 18, face = "bold"),
            legend.text = element_text(size=16),
            strip.background = element_blank(),
            strip.text = element_text(size = 14, face = "bold"),
            #axis.title = element_blank(),
            #axis.text = element_blank(),
            #axis.ticks = element_blank(),
            #panel.grid = element_blank(),
            #strip.background = element_blank(),
            #strip.text = element_text(size=16, face = "bold"),
            #legend.title = element_text(size=18),
            #legend.text = element_text(size=16),
            #plot.title = element_text(hjust = 0.5)
        ) + 
        scale_color_discrete(name=color_legend_title) + 
        guides(color = guide_legend(
                                    override.aes = list(size = 3)
        ))
    if (!is.null(tp_color_palette)) {
        g = g + scale_color_manual(values=tp_color_palette, name=color_legend_title)
    }
    ggsave(file_name, height=height, width=width)
}

plot_OT_heatmap = function(timepoints, method, method_result_dir, tp_color_palette=NULL, prob_cutoff=0.2) {
    require(reshape2)
    if (method == 'spaDOT') {method = 'GATVAE-SVGP-VAE'}
    # calculate marker overlap
    timepoints_dict = c('t0'='E12.5', 't1'='E13.5', 't2'='E15.5')
    cluster_mapping = list('t0'=c('0'='5', '1'='1', '2'='4', '3'='3', '4'='2'),
                            't1'=c('0'='4', '1'='1', '2'='3', '3'='5', '4'='2'),
                            't2'=c('0'='5-1', '1'='1', '2'='5-2', '3'='4', '4'='2'))
    for (tp_i in 1:(length(timepoints)-1)) {
        cur_tp = timepoints[tp_i]
        next_tp = timepoints[tp_i+1]
        OT_mat_row_norm = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '_row_norm.csv')), row.names=1)
        OT_mat_col_norm = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '_col_norm.csv')), row.names=1)
        # obtain minimum value for each entry
        OT_mat = data.frame(mapply(pmin, OT_mat_row_norm, OT_mat_col_norm))
        write.csv(OT_mat, file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_OT_mat_minimum.csv')), quote=F)
        # OT_mat[OT_mat < prob_cutoff] = 0
        rownames(OT_mat) = as.character(1:nrow(OT_mat)-1)
        colnames(OT_mat) = as.character(1:ncol(OT_mat)-1)
        # rename 
        rownames(OT_mat) = unname(cluster_mapping[[cur_tp]][rownames(OT_mat)])
        colnames(OT_mat) = unname(cluster_mapping[[next_tp]][colnames(OT_mat)])
        # sort rownames and colnames
        OT_mat = OT_mat[order(rownames(OT_mat)), order(colnames(OT_mat))]
        x_axis_colors = unname(tp_color_palette[colnames(OT_mat)])
        y_axis_colors = unname(tp_color_palette[rownames(OT_mat)])
        OT_mat$ID = rownames(OT_mat)
        OT_mat_melt = melt(OT_mat, id.vars = "ID")
        colnames(OT_mat_melt) = c(cur_tp, next_tp, 'value')
        OT_mat_melt$significant = ifelse(OT_mat_melt$value < prob_cutoff, 0, 1)
        # OT_mat_melt$color <- ifelse(OT_mat_melt$value < 0.2, "< 0.2", ">= 0.2")
        OT_mat_melt$value_edited = ifelse(OT_mat_melt$value < prob_cutoff, NA, OT_mat_melt$value)
        g = ggplot(OT_mat_melt, aes(!!sym(next_tp), !!sym(cur_tp), size=value, fill=value_edited)) +
            # geom_tile() + 
            geom_point(shape=21) + 
            # geom_text(aes(label = ifelse(value < prob_cutoff, "", round(value, 3))), color = "black") +
            scale_size_continuous(range = c(1, 8), name = "Transition Probability") +
            # scale_color_manual(values = c(">= 0.2" = "red", "< 0.2" = "grey")) +
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
        ggsave(file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_OT_heatmap.pdf')), width=4, height=4)
        ggsave(file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_OT_heatmap_legend.pdf')), cowplot::plot_grid(legend), width = 2, height = 4)
    }
}


cluster_marker_genes = function(Seurat_obj, method, method_result_dir) {
    require(Seurat)
    timepoints = as.vector(unique(Seurat_obj@meta.data$timepoint))
    timepoints = timepoints[order(timepoints)]

    if (method == 'spaDOT') {method_col='pred_labels'}
    else {method_col=paste0(method, '_pred_labels')}

    # remove ribosomal genes if any
    ribosomal_genes = grep("^Rp[sl]", rownames(Seurat_obj), value=T)
    Seurat_obj =  Seurat_obj[!rownames(Seurat_obj) %in% ribosomal_genes, ]
    mt_genes = grep("^mt-", rownames(Seurat_obj), value=T)
    Seurat_obj =  Seurat_obj[!rownames(Seurat_obj) %in% mt_genes, ]
    for (tp in timepoints) {
        tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
        Idents(tp_Seurat_obj) = tp_Seurat_obj@meta.data[[method_col]]
        tp_markers = FindAllMarkers(tp_Seurat_obj)
        write.csv(tp_markers, file.path(method_result_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')), quote=F)    
    }
}

spearman_corr_markeroverlap_vs_ot = function(timepoints, method, method_result_dir, perform_GO=FALSE, tp_color_palette=NULL) {
    require(clusterProfiler)
    require(org.Mm.eg.db)
    require(reshape2)
    # calculate overlap between marker genes
    tp_marker_list = list()
    cluster_mapping = list('t0'=c('0'='5', '1'='1', '2'='4', '3'='3', '4'='2'),
                            't1'=c('0'='4', '1'='1', '2'='3', '3'='5', '4'='2'),
                            't2'=c('0'='5-1', '1'='1', '2'='5-2', '3'='4', '4'='2'))
    for (tp in timepoints) {
        tp_marker_list[[tp]] = list()
        tp_marker = read.csv(file.path(method_result_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')), row.names=1)
        sub_tp_marker = tp_marker[tp_marker$avg_log2FC > 1 & tp_marker$p_val_adj < 0.05,]
        for (cluster in unique(sub_tp_marker$cluster)) {
            tp_marker_list[[tp]][[as.character(cluster)]] = sub_tp_marker[sub_tp_marker$cluster==cluster, 'gene']
            if (perform_GO) {
                # tp_cluster_GO = enrichGO(gene = sub_tp_marker[sub_tp_marker$cluster==cluster, 'gene'],
                #                 keyType = 'SYMBOL',
                #                 OrgDb = org.Mm.eg.db, 
                #                 ont = "BP", # Biological process
                #                 pAdjustMethod = "BH", 
                #                 qvalueCutoff = 0.05)
                # g = dotplot(tp_cluster_GO, showCategory=10)
                # ggsave(file.path(method_result_dir, paste0(tp, '_', cluster, '_enrichment_clusterprofiler.png')), width=8, height=6)

                tp_cluster_GO = enrichGO(gene = sub_tp_marker[sub_tp_marker$cluster==cluster, 'gene'][1:100], 
                                keyType = 'SYMBOL',
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP", # Biological process
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05)
                # --- write results into file
                write.csv(as.data.frame(tp_cluster_GO), file.path(method_result_dir, paste0(tp, '_', cluster, '_top100_enrichment_clusterprofiler.csv')), quote=F)
                g = dotplot(tp_cluster_GO, showCategory=10)
                ggsave(file.path(method_result_dir, paste0(tp, '_', cluster, '_top100_enrichment_clusterprofiler.png')), width=8, height=6)
            }
        }
    }

    # --- calculating number of overlapping GO terms
    for (tp_i in 1:(length(timepoints)-1)) {
        cur_tp = timepoints[tp_i]
        cur_tp_GO_terms_files = list.files(method_result_dir, pattern=paste0(cur_tp, '_.*_top100_enrichment_clusterprofiler.csv'), full.names=T)
        next_tp = timepoints[tp_i+1]
        next_tp_GO_terms_files = list.files(method_result_dir, pattern=paste0(next_tp, '_.*_top100_enrichment_clusterprofiler.csv'), full.names=T)
        
        cur_tp_next_tp_mat = matrix(nrow=length(cur_tp_GO_terms_files), ncol=length(next_tp_GO_terms_files))
        rownames(cur_tp_next_tp_mat) = unname(sapply(cur_tp_GO_terms_files, function(x) {gsub('_top100_enrichment_clusterprofiler.csv', '', basename(x))}))
        colnames(cur_tp_next_tp_mat) = unname(sapply(next_tp_GO_terms_files, function(x) {gsub('_top100_enrichment_clusterprofiler.csv', '', basename(x))}))

        for (cur_tp_GO_term_file in cur_tp_GO_terms_files) {
            cur_tp_GO = read.csv(cur_tp_GO_term_file, header=T)[1:100, ]
            cur_tp_cluster = gsub('_top100_enrichment_clusterprofiler.csv', '', basename(cur_tp_GO_term_file))
            for (next_tp_GO_term_file in next_tp_GO_terms_files) {
                next_tp_GO = read.csv(next_tp_GO_term_file, header=T)[1:100, ]
                next_tp_cluster = gsub('_top100_enrichment_clusterprofiler.csv', '', basename(next_tp_GO_term_file))
                overlap_GO = unique(intersect(cur_tp_GO$ID, next_tp_GO$ID))
                cur_tp_next_tp_mat[cur_tp_cluster, next_tp_cluster] = length(overlap_GO)
            }
        }
        timepoints_dict = c('t0'='E12.5', 't1'='E13.5', 't2'='E15.5')
        cur_tp_next_tp_mat = as.data.frame(cur_tp_next_tp_mat)
        # --- rename row and column names
        rownames(cur_tp_next_tp_mat) = sapply(strsplit(rownames(cur_tp_next_tp_mat), '_'), '[', 2)
        colnames(cur_tp_next_tp_mat) = sapply(strsplit(colnames(cur_tp_next_tp_mat), '_'), '[', 2)
        rownames(cur_tp_next_tp_mat) = unname(cluster_mapping[[cur_tp]][rownames(cur_tp_next_tp_mat)])
        colnames(cur_tp_next_tp_mat) = unname(cluster_mapping[[next_tp]][colnames(cur_tp_next_tp_mat)])
        # sort rownames and colnames
        cur_tp_next_tp_mat = cur_tp_next_tp_mat[order(rownames(cur_tp_next_tp_mat)), order(colnames(cur_tp_next_tp_mat))]
        cur_tp_next_tp_mat$ID = rownames(cur_tp_next_tp_mat)
        cur_tp_next_tp_mat_melt = melt(cur_tp_next_tp_mat, id.vars = "ID")
        colnames(cur_tp_next_tp_mat_melt) = c(cur_tp, next_tp, 'value')
        cur_tp_next_tp_mat_melt[[cur_tp]] = gsub(cur_tp, unname(timepoints_dict[cur_tp]), cur_tp_next_tp_mat_melt[[cur_tp]])
        cur_tp_next_tp_mat_melt[[next_tp]] = gsub(next_tp, unname(timepoints_dict[next_tp]), cur_tp_next_tp_mat_melt[[next_tp]])
        
        x_axis_colors = unname(tp_color_palette[colnames(cur_tp_next_tp_mat)])
        y_axis_colors = unname(tp_color_palette[rownames(cur_tp_next_tp_mat)])
        g = ggplot(cur_tp_next_tp_mat_melt, aes(!!sym(next_tp), !!sym(cur_tp), size=value, fill=value)) +
            geom_point(shape=21) + 
            scale_size_continuous(range = c(1, 8), breaks=c(25, 50, 75), name = "No. GOs") +
            scale_fill_gradientn(
                colors = c("white", "red"),
                limits = c(0, 100),
                na.value = "grey",
                name = "Color scale"
            ) +
            labs(x = paste0(unname(timepoints_dict[next_tp]), " Domains"), 
                y = paste0(unname(timepoints_dict[cur_tp]), " Domains")) +
            theme_pubr() + 
            theme(legend.position="bottom",
                axis.text.x = element_text(color=x_axis_colors, angle=45, hjust=1),
                axis.text.y = element_text(color=y_axis_colors)
            )
        legend = get_legend(g)
        g = g + theme(legend.position="none")
        ggsave(file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_GO_overlap_heatmap.pdf')), width=4, height=4)
        ggsave(file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_GO_overlap_heatmap_legend.pdf')), cowplot::plot_grid(legend), width = 6, height = 2)
    }

    # --- Venn diagram on certain clusters at certain timepoints..
    t0_2_go = read.csv(file.path(method_result_dir, paste0('t0', '_', 2, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t0_2_go = t0_2_go[t0_2_go$p.adjust < 0.05,]
    t1_0_go = read.csv(file.path(method_result_dir, paste0('t1', '_', 0, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t1_0_go = t1_0_go[t1_0_go$p.adjust < 0.05,]
    t2_3_go = read.csv(file.path(method_result_dir, paste0('t2', '_', 3, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t2_3_go = t2_3_go[t2_3_go$p.adjust < 0.05,]
    # remove NA
    go_list = list(t0_2_go$ID, t1_0_go$ID, t2_3_go$ID)
    go_list = lapply(go_list, function(x) {
        x = x[!grepl("NA", x)]
        x = x[!is.na(x)]
        x = unique(x)
        x
    })
    library(ggvenn)
    names(go_list) = c('E12.5_2', 'E13.5_0', 'E15.5_3')
    g = ggvenn(go_list)
    ggsave(file.path(method_result_dir, 'E12.5_2-E13.5_0-E15.5_3_venn.png'), width=4, height=4)

    # --- permutation test
    list1 = go_list[[1]]
    list2 = go_list[[2]]
    list3 = go_list[[3]]
    perm_N = 1000
    size1 = length(list1)
    size2 = length(list2)
    size3 = length(list3)
    total_GOs = Reduce(union, list(list1, list2, list3))
    overlap_perm = rep(NA, perm_N)
    for (perm_i in 1:perm_N) {
        perm_list1 = sample(total_GOs, size1)
        perm_list2 = sample(total_GOs, size2)
        perm_list3 = sample(total_GOs, size3)
        overlap_perm[perm_i] = length(intersect(perm_list1, intersect(perm_list2, perm_list3)))
    }
    # Z statistics
    z = (length(Reduce(intersect, go_list)) - mean(overlap_perm)) / sd(overlap_perm)
    p = 2 * (1 - pnorm(abs(z)))


    contingency_table = matrix(c(
        length(go_list[['E12.5_2']]), length(intersect(go_list[['E12.5_2']], go_list[['E13.5_0']])), length(intersect(go_list[['E12.5_2']], go_list[['E15.5_3']])),
        length(go_list[['E13.5_0']]), length(intersect(go_list[['E13.5_0']], go_list[['E15.5_3']])), length(unique(Reduce(intersect, go_list))),
        length(go_list[['E15.5_3']]), length(intersect(go_list[['E15.5_3']], go_list[['E12.5_2']])), length(intersect(go_list[['E15.5_3']], go_list[['E13.5_0']]))
    ), ncol=3)
    chi_square_result <- chisq.test(contingency_table)
    print(chi_square_result)

    t0_1_go = read.csv(file.path(method_result_dir, paste0('t0', '_', 1, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t0_1_go = t0_1_go[t0_1_go$p.adjust < 0.05,]
    t1_1_go = read.csv(file.path(method_result_dir, paste0('t1', '_', 1, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t1_1_go = t1_1_go[t1_1_go$p.adjust < 0.05,]
    t2_1_go = read.csv(file.path(method_result_dir, paste0('t2', '_', 1, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t2_1_go = t2_1_go[t2_1_go$p.adjust < 0.05,]
    # remove NA
    go_list = list(t0_1_go$ID, t1_1_go$ID, t2_1_go$ID)
    go_list = lapply(go_list, function(x) {
        x = x[!grepl("NA", x)]
        x = x[!is.na(x)]
        x = unique(x)
        x
    })
    names(go_list) = c('E12.5_1', 'E13.5_1', 'E15.5_1')
    g = ggvenn(go_list)
    ggsave(file.path(method_result_dir, 'E12.5_1-E13.5_1-E15.5_1_venn.png'), width=4, height=4)

    # --- permutation test
    list1 = go_list[[1]]
    list2 = go_list[[2]]
    list3 = go_list[[3]]
    perm_N = 1000
    size1 = length(list1)
    size2 = length(list2)
    size3 = length(list3)
    total_GOs = Reduce(union, list(list1, list2, list3))
    overlap_perm = rep(NA, perm_N)
    for (perm_i in 1:perm_N) {
        perm_list1 = sample(total_GOs, size1)
        perm_list2 = sample(total_GOs, size2)
        perm_list3 = sample(total_GOs, size3)
        overlap_perm[perm_i] = length(intersect(perm_list1, intersect(perm_list2, perm_list3)))
    }
    # Z statistics
    z = (length(Reduce(intersect, go_list)) - mean(overlap_perm)) / sd(overlap_perm)
    p = 2 * (1 - pnorm(abs(z)))

    # # generate contingency table
    # contingency_table = matrix(c(
    #     length(go_list[['E12.5_1']]), length(intersect(go_list[['E12.5_1']], go_list[['E13.5_1']])), length(intersect(go_list[['E12.5_1']], go_list[['E15.5_1']])),
    #     length(go_list[['E13.5_1']]), length(intersect(go_list[['E13.5_1']], go_list[['E15.5_1']])), length(Reduce(intersect, go_list)),
    #     length(go_list[['E15.5_1']]), length(intersect(go_list[['E15.5_1']], go_list[['E12.5_1']])), length(intersect(go_list[['E15.5_1']], go_list[['E13.5_1']]))
    # ), ncol=3)
    # chi_square_result <- chisq.test(contingency_table)
    # print(chi_square_result)



    t2_0_go = read.csv(file.path(method_result_dir, paste0('t2', '_', 0, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t2_0_go = t2_0_go[t2_0_go$p.adjust < 0.05,]
    t2_2_go = read.csv(file.path(method_result_dir, paste0('t2', '_', 2, '_top100_enrichment_clusterprofiler.csv')), header=T)
    t2_2_go = t2_2_go[t2_2_go$p.adjust < 0.05,]
    t2_4_go = read.csv(file.path(method_result_dir, paste0('t2', '_', 4, '_top100_enrichment_clusterprofiler.csv')), header=T) 
    t2_4_go = t2_4_go[t2_4_go$p.adjust < 0.05,]
    # remove NA
    go_list = list(t2_0_go$ID, t2_2_go$ID, t2_4_go$ID)
    go_list = lapply(go_list, function(x) {
        x = x[!grepl("NA", x)]
        x = unique(x)
        x
    })
    names(go_list) = c('E15.5_0', 'E15.5_2', 'E15.5_4')
    g = ggvenn(go_list)
    ggsave(file.path(method_result_dir, 'E15.5_0-E15.5_2-E15.5_4_venn.png'), width=4, height=4)

    # find unique GO terms in E15.5_4
    unique_go = setdiff(go_list[['E15.5_4']], Reduce(union, go_list[-3]))
    t2_4_go_sub = t2_4_go[t2_4_go$ID %in% unique_go,]
    t2_4_go_sub = t2_4_go_sub[t2_4_go_sub$p.adjust < 0.05,]

    if (method == 'spaDOT') {method = 'GATVAE-SVGP-VAE'}
    # calculate marker overlap
    overlap_spear = list()
    for (tp_i in 1:(length(timepoints)-1)) {
        cur_tp = timepoints[tp_i]
        next_tp = timepoints[tp_i+1]
        cur_tp_marker_list = tp_marker_list[[cur_tp]]
        next_tp_marker_list = tp_marker_list[[next_tp]]
        overlap_markers = sapply(next_tp_marker_list, function(cluster2) {
            sapply(cur_tp_marker_list, function(cluster1) {
                length(intersect(cluster1, cluster2)) # Find overlapping genes
            })
        })
        # obtain OT matrix
        OT_mat = read.csv(file.path(method_result_dir, paste0(method, '_gamma_', 'tp', cur_tp, '_', 'tp', next_tp, '.csv')), row.names=1)
        colnames(OT_mat) = 1:ncol(OT_mat)-1
        # calculate spearman correlation per column between overlap_markers and OT_mat
        # spearman_corr = apply(overlap_markers, 2, function(x) {
        #     apply(OT_mat, 2, function(y) {
        #         cor(x, y, method='spearman')
        #     })
        # })
        overlap_spear[[paste0(cur_tp, '_', next_tp)]] = cor(as.vector(as.matrix(overlap_markers)), as.vector(as.matrix(OT_mat)), method='pearson') 
    }
    return(overlap_spear)
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
    require(dplyr)
    require(ggplot2)
    require(reshape2)
    timepoints = as.vector(unique(Seurat_obj@meta.data$timepoint))
    timepoints = timepoints[order(timepoints)]
    tp_count_list = list()
    for (tp in timepoints) {
        tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
        tp_Seurat_obj = SCTransform(tp_Seurat_obj, verbose = FALSE) ## according to axolotl brain Fig S6D
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
        library(pheatmap)
        pheatmap(tp_cor_mat, color= colorRampPalette(c("blue", "white", "red"))(100), main="Gene expression correlation heatmap", cluster_rows = T, cluster_cols = T,
        #scale="row", 
        filename=file.path(method_result_dir, paste0(tp, '_gene_correlation_heatmap.png')))

    }


    # calculate intersected genes
    common_genes = Reduce(intersect, lapply(tp_count_list, rownames))
    subset_count_list = lapply(tp_count_list, function(mat) mat[common_genes, , drop = FALSE])
    # calculate per group
    tp_group_avg_exprs_res = NULL
    for (tp in timepoints) {
        tp_subset_count = subset_count_list[[tp]]
        tp_subset_count_t = as.data.frame(t(tp_subset_count))
        tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
        tp_subset_count_t$spaDOT_pred_labels = tp_Seurat_obj@meta.data[rownames(tp_subset_count_t), 'spaDOT_pred_labels']
        tp_group_avg_exprs = tp_subset_count_t %>%
        group_by(spaDOT_pred_labels) %>% summarise(across(everything(), mean))
        tp_group_avg_exprs_df = as.data.frame(t(tp_group_avg_exprs[-1]))
        colnames(tp_group_avg_exprs_df) = tp_group_avg_exprs$spaDOT_pred_labels
        if (is.null(tp_group_avg_exprs_res)) {
            tp_group_avg_exprs_res = tp_group_avg_exprs_df
        } else {
            tp_group_avg_exprs_res = cbind(tp_group_avg_exprs_res, tp_group_avg_exprs_df)
        }
        # calculate correlation and plot for each timepoint
        cor_mat = cor(tp_group_avg_exprs_df)
        x_axis_colors = unname(tp_color_palette[rownames(cor_mat)])
        y_axis_colors = unname(tp_color_palette[colnames(cor_mat)])
        cor_mat_melt = melt(cor_mat)
        g = ggplot(cor_mat_melt, aes(x=Var1, y=Var2, fill = value)) +
            geom_tile() +
            geom_text(aes(label = round(value, 3)), color = "white") +
            scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
            labs(x = "",  y = "") +
            theme_pubr() + 
            theme(legend.position="right",
                axis.text.x = element_text(color=x_axis_colors, angle=45, hjust=1),
                axis.text.y = element_text(color=y_axis_colors)
            )
        ggsave(file.path(method_result_dir, paste0(tp, '_gene_correlation_heatmap.png')), width=5, height=4)
    }
    # # calculate correlation and plot
    cor_mat = cor(tp_group_avg_exprs_res, method='spearman')
    library(pheatmap)
    pheatmap(cor_mat, color= colorRampPalette(c("blue", "white", "red"))(100), main="Gene expression correlation heatmap", cluster_rows = T, cluster_cols = T,
    #scale="row", 
    filename=file.path(method_result_dir, 'gene_correlation_heatmap.png'))




    x_axis_colors = unname(tp_color_palette[rownames(cor_mat)])
    y_axis_colors = unname(tp_color_palette[colnames(cor_mat)])
    cor_mat_melt = melt(cor_mat)
    g = ggplot(cor_mat_melt, aes(x=Var1, y=Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
    labs(x = "",  y = "") +
    theme_pubr() + 
    theme(legend.position="right",
        axis.text.x = element_text(color=x_axis_colors, angle=45, hjust=1),
        axis.text.y = element_text(color=y_axis_colors)
    )
    ggsave(file.path(method_result_dir, 'gene_correlation_heatmap.png'), width=5, height=4)






}
