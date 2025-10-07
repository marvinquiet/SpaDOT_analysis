## conda activate /media/sda/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(ggplot2)
library(ggpubr)

prefix_dir = "/net/mulan/home/wenjinma"
data_dir = file.path(prefix_dir, 'data/spatial_scRNAseq_GSE214989')
spatialPCA_result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_GSE214989_multiplesample_combinedlabels_removedcells')
spaDOT_result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_GSE214989_multiplesample_combinedlabels/fgsea_for_spaDOT_labels_afterremoval')
spatialPCA_res_df = read.csv(file.path(spatialPCA_result_dir, 'spatialPCA_results.csv'), header=T, row.names=1)
rownames(spatialPCA_res_df) = spatialPCA_res_df$cell_names
spatialPCA_res_df$section = sapply(strsplit(spatialPCA_res_df$cell_names, split='_'), '[', 2)
section_to_tp = c('1'='t2', '2'='t1', '3'='t0')
spatialPCA_res_df$timepoint = unname(section_to_tp[spatialPCA_res_df$section])

# match for each timepont
for (timepoint in c('t0', 't1', 't2')) {
    tp_spaDOT = read.csv(file.path(spaDOT_result_dir, paste0(timepoint, '_Louvain_labels.csv')), header=T, row.names=1)
    tp_spatialPCA = spatialPCA_res_df[spatialPCA_res_df$timepoint == timepoint, ]
    common_cells = intersect(rownames(tp_spaDOT), rownames(tp_spatialPCA))
    if (length(common_cells) == 0) {
        cat('There must be some matching errors between spatialPCA and spaDOT, please check!')
        next
    }
    tp_spaDOT = tp_spaDOT[common_cells, ]
    tp_spatialPCA = tp_spatialPCA[common_cells, ]
    cluster_tables = table(tp_spaDOT$louvain, tp_spatialPCA$spatialPCA_labels)
    print(timepoint)
    print(cluster_tables)
    cluster_df = as.data.frame(cluster_tables)
    colnames(cluster_df) = c('spaDOT', 'spatialPCA', 'value')
    ggplot(cluster_df, aes(x=spaDOT, y=spatialPCA, fill=value)) + geom_tile() + 
        geom_text(aes(label = value), color = "white", size = 4)
    ggsave(file.path(spatialPCA_result_dir, paste0(timepoint, '_spatialPCA_spaDOT_heatmap.png')))
}
