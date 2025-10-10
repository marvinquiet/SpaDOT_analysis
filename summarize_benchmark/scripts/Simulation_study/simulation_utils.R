rename_GT_labels = function(latent_adata) {
    latent_adata$obs = latent_adata$obs %>%
    mutate(renamed_GT_labels = case_when(
        timepoint == 't1' & true_labels == 1 ~ "2",
        timepoint == 't1' & true_labels == 2 ~ "1-2",
        timepoint == 't1' & true_labels == 3 ~ "1-1",
        timepoint == 't1' & true_labels == 4 ~ "3",
        timepoint == 't2' & true_labels == 1 ~ "2",
        timepoint == 't2' & true_labels == 2 ~ "3-1",
        timepoint == 't2' & true_labels == 3 ~ "1-2",
        timepoint == 't2' & true_labels == 4 ~ "1-1",
        timepoint == 't2' & true_labels == 5 ~ "3-2",
        TRUE ~ as.character(true_labels)  # Default case
    ))
    return(latent_adata$obs$renamed_GT_labels)
}

GT_custom_palette = function() {
    return(
        c('1'= '#EA3423',
        '2'= '#76FC4C',
        '3'= '#FEFC53',
        '1-1'= '#891AF5',
        '1-2'='#EA33F4',
        '3-1'= '#716A1E',
        '3-2'= '#F5C242')
    )
}

custom_palette_2030 = c(
    't0_0'='#EA3423',
    't0_1'='#76FC4C',
    't0_2'='#FEFC53',
    't1_0'='#EA33F4',
    't1_1'='#FEFC53',
    't1_2'='#891AF5',
    't1_3'='#76FC4C',
    't2_0'='#76FC4C',
    't2_4'='#F5C242',
    't2_2'='#891AF5',
    't2_1'='#716A1E',
    't2_3'='#EA33F4'
)

custome_palette_list = list(
    '1993/20_30'=c(
        't0_0'='#EA3423',
        't0_1'='#76FC4C',
        't0_2'='#FEFC53',
        't1_0'='#EA33F4',
        't1_1'='#FEFC53',
        't1_2'='#891AF5',
        't1_3'='#76FC4C',
        't2_0'='#76FC4C',
        't2_4'='#F5C242',
        't2_2'='#891AF5',
        't2_1'='#716A1E',
        't2_3'='#EA33F4'
        ),
    '1993/20_50'=c(
        't0_0'='#EA3423',
        't0_1'='#FEFC53',
        't0_2'='#76FC4C',
        't1_0'='#EA33F4',
        't1_1'='#FEFC53',
        't1_2'='#76FC4C',
        't1_3'='#891AF5',
        't2_0'='#716A1E',
        't2_4'='#F5C242',
        't2_2'='#76FC4C',
        't2_1'='#891AF5',
        't2_3'='#EA33F4'
        ),
    '1993/30_30'=c(
        't0_0'='#EA3423',
        't0_1'='#FEFC53',
        't0_2'='#76FC4C',
        't1_0'='#EA33F4',
        't1_1'='#FEFC53',
        't1_2'='#76FC4C',
        't1_3'='#891AF5',
        't2_0'='#76FC4C',
        't2_4'='#F5C242',
        't2_2'='#716A1E',
        't2_1'='#891AF5',
        't2_3'='#EA33F4'
        ),
    '1993/30_50'=c(
        't0_0'='#EA3423',
        't0_1'='#FEFC53',
        't0_2'='#76FC4C',
        't1_0'='#EA33F4',
        't1_1'='#FEFC53',
        't1_2'='#76FC4C',
        't1_3'='#891AF5',
        't2_0'='#716A1E',
        't2_4'='#F5C242',
        't2_2'='#EA33F4',
        't2_1'='#891AF5',
        't2_3'='#76FC4C'
    )
)

plot_OT_heatmap = function(timepoints, method_result_dir, tp_color_palette=NULL, prob_cutoff=0.2, perm_N = 100) {
    require(reshape2)
    require(ggsankey) ## networkD3 visualization is not good
    # calculate marker overlap
    sankey_df = NULL
    for (tp_i in 1:(length(timepoints)-1)) {
        # obtain OT matrix
        # OT_mat = read.csv(file.path(method_result_dir, paste0(method, '_gamma_', 'tp', cur_tp, '_', 'tp', next_tp, '.csv')), row.names=1)
        cur_tp = timepoints[tp_i]
        next_tp = timepoints[tp_i+1]
        OT_mat = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '.csv')), row.names=1)
        OT_mat[OT_mat < prob_cutoff] = 0
        OT_mat$ID = rownames(OT_mat)
        OT_mat_melt = melt(OT_mat, id.vars = "ID")
        colnames(OT_mat_melt) = c(cur_tp, next_tp, 'value')
        x_axis_colors = unname(tp_color_palette[as.vector(sort(unique(OT_mat_melt[,next_tp])))])
        y_axis_colors = unname(tp_color_palette[as.vector(sort(unique(OT_mat_melt[,cur_tp])))])
        ggplot(OT_mat_melt, aes(!!sym(next_tp), !!sym(cur_tp), fill = value)) +
            geom_tile() +
            geom_text(aes(label = ifelse(value < prob_cutoff, "", round(value, 3))), color = "black") +
            scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
            labs(x = paste0(next_tp, " Clusters"), 
                y = paste0(cur_tp, " Clusters")) +
            theme_pubr() + 
            theme(legend.position="right",
                axis.text.x = element_text(color=x_axis_colors, angle=45, hjust=1),
                axis.text.y = element_text(color=y_axis_colors)
            )
        ggsave(file.path(method_result_dir, paste0(cur_tp, '_', next_tp, '_OT_heatmap.png')), width=5, height=4)

        # --- plot sankey diagram
        tp_i_sankey_df = OT_mat_melt
        tp_i_sankey_df = tp_i_sankey_df[tp_i_sankey_df$value > prob_cutoff, ]
        tp_i_sankey_df$times = round(tp_i_sankey_df$value * perm_N)
        tp_i_sankey_df = as.data.frame(lapply(tp_i_sankey_df, rep, tp_i_sankey_df$times))
        tp_i_sankey_df$value = NULL
        tp_i_sankey_df$times = NULL
        tp_i_sankey_df = tp_i_sankey_df %>% make_long(!!sym(cur_tp), !!sym(next_tp))
        if (is.null(sankey_df)) {
            sankey_df = tp_i_sankey_df
        }
        else {
            sankey_df = rbind(sankey_df, tp_i_sankey_df)
        }      
    }
    ggplot(sankey_df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
        geom_sankey() +
        theme_sankey(base_size = 16)
    ggsave(file.path(method_result_dir, 'sankey.png'), width=5, height=4)
}

plot_Sankey = function(timepoints, method_result_dir, tp_color_palette=NULL, prob_cutoff=0.2, perm_N = 100) {
    require(reshape2)
    require(ggsankey) ## networkD3 visualization is not good
    # calculate marker overlap
    sankey_df = NULL
    for (tp_i in 1:(length(timepoints)-1)) {
        # obtain OT matrix
        cur_tp = timepoints[tp_i]
        next_tp = timepoints[tp_i+1]
        OT_mat = read.csv(file.path(method_result_dir, paste0('transition_prob_', gsub('t', '', cur_tp), '_', gsub('t', '', next_tp), '.csv')), row.names=1)
        OT_mat[OT_mat < prob_cutoff] = 0
        OT_mat$ID = rownames(OT_mat)
        OT_mat_melt = melt(OT_mat, id.vars = "ID")
        colnames(OT_mat_melt) = c(cur_tp, next_tp, 'value')
        # --- plot sankey diagram
        tp_i_sankey_df = OT_mat_melt
        tp_i_sankey_df = tp_i_sankey_df[tp_i_sankey_df$value > prob_cutoff, ]
        tp_i_sankey_df$times = round(tp_i_sankey_df$value * perm_N)
        tp_i_sankey_df = as.data.frame(lapply(tp_i_sankey_df, rep, tp_i_sankey_df$times))
        tp_i_sankey_df$value = NULL
        tp_i_sankey_df$times = NULL
        tp_i_sankey_df = tp_i_sankey_df %>% make_long(!!sym(cur_tp), !!sym(next_tp))
        if (is.null(sankey_df)) {
            sankey_df = tp_i_sankey_df
        }
        else {
            sankey_df = rbind(sankey_df, tp_i_sankey_df)
        }      
    }
    ggplot(sankey_df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
        geom_sankey() +
        theme_sankey(base_size = 16)
    ggsave(file.path(method_result_dir, 'sankey.png'), width=5, height=4)
}