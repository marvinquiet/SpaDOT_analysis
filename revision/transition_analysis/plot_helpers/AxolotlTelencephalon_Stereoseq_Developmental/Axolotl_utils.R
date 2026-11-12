# Dataset-specific plot helper used for the manuscript revision.
plot_OT_heatmap = function(timepoints, method, method_result_dir, tp_color_palette=NULL, prob_cutoff=0.2) {
    require(reshape2)
    if (method == 'spaDOT') {method = 'GATVAE-SVGP-VAE'}
    # calculate marker overlap
    timepoints_dict = c('t0'='St.44', 't1'='St.54', 't2'='St.57')
    cluster_mapping = list(
    'St.44'=c("0"="6", "1"="5", "2"="3", "3"="2", "4"="4", "5"="1"),
    'St.54'=c("0"="1", "1"="5", "2"="3", "3"="4-1", "4"="4-2", "5"="6", "6"="2"),
    'St.57'=c("0"="4-1", "1"="2", "2"="4-2-1", "3"="6", "4"="5", "5"="4-2-2", "6"="3", "7"="1")
    )
    for (tp_i in 1:(length(timepoints)-1)) {
        cur_tp = timepoints[tp_i]
        next_tp = timepoints[tp_i+1]
        # obtain OT matrix
        OT_mat_row_norm = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '_row_norm.csv')), row.names=1)
        OT_mat_col_norm = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '_col_norm.csv')), row.names=1)
        # obtain minimum value for each entry
        OT_mat = data.frame(mapply(pmin, OT_mat_row_norm, OT_mat_col_norm))
        # write.csv(OT_mat, file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_OT_mat_minimum.csv')), quote=F)
        # OT_mat[OT_mat < prob_cutoff] = 0
        rownames(OT_mat) = 1:nrow(OT_mat)-1
        colnames(OT_mat) = 1:ncol(OT_mat)-1
        rownames(OT_mat) = unname(cluster_mapping[[timepoints_dict[cur_tp]]][rownames(OT_mat)])
        colnames(OT_mat) = unname(cluster_mapping[[timepoints_dict[next_tp]]][colnames(OT_mat)])
        OT_mat = OT_mat[order(rownames(OT_mat)), order(colnames(OT_mat))]
        x_axis_colors = unname(tp_color_palette[colnames(OT_mat)])
        y_axis_colors = unname(tp_color_palette[rownames(OT_mat)])

        OT_mat$ID = rownames(OT_mat)
        OT_mat_melt = melt(OT_mat, id.vars = "ID")
        colnames(OT_mat_melt) = c(cur_tp, next_tp, 'value')
        OT_mat_melt$value_edited = ifelse(OT_mat_melt$value < prob_cutoff, NA, OT_mat_melt$value)
        g = ggplot(OT_mat_melt, aes(!!sym(next_tp), !!sym(cur_tp), size=value, fill=value_edited)) +
            # geom_tile() + 
            geom_point(shape=21) + 
            # geom_text(aes(label = ifelse(value < prob_cutoff, "", round(value, 3))), color = "black") +
            scale_size_continuous(range = c(1, 8), breaks=c(0.25, 0.50, 0.75), name = "Transition Probability") +
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
