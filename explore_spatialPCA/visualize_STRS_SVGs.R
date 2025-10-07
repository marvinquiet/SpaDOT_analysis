.libPaths("/net/mulan/home/wenjinma/Rlib")
library(SpatialPCA)
library(SeuratDisk)
library(Seurat)
#library(Signac)
library(aricode)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(tidyverse)
library(VennDiagram)

prefix_dir = "/net/mulan/home/wenjinma"

data_dir = file.path(prefix_dir, 'data/spatial_STRS/')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_STRS_multiplesample_combinedlabels')
dir.create(result_dir, showWarnings = F)

Seurat_obj = readRDS(file.path(data_dir, 'STRS.RDS'))
Seurat_obj = JoinLayers(Seurat_obj)
ggplot(Seurat_obj@meta.data, aes(factor(timepoint), log10(nFeature_RNA+1))) + geom_violin()
ggsave(file.path(result_dir, 'nFeature_RNA_per_timepoint_log10.png'))
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)

# Visualize ground truth
ggplot(Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=injury.zones)) + geom_point() + 
    facet_wrap(~timepoint, scales="free") + theme_pubr()
ggsave(file.path(result_dir, 'Groundtruth_annotation.png'), height=5, width=12)



# --- explore marker genes in 2dpi and 5dpi
tp_region_markers_list = list()
for (tp in c('t0', 't1')) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    Idents(tp_Seurat_obj) = tp_Seurat_obj$injury.zones
    tp_markers = FindAllMarkers(tp_Seurat_obj)
    write.csv(tp_markers, file.path(result_dir, paste0(tp, '_Seurat_FindAllMarkers.csv')), quote=F)
    subset_tp_markers = tp_markers[abs(tp_markers$avg_log2FC) > 0 & tp_markers$p_val_adj < 0.01, ]
    write.csv(subset_tp_markers, file.path(result_dir, paste0(tp, '_Seurat_FindAllMarkers_log2FC_0_adjp_0.01.csv')), quote=F)
    # region specific marker genes
    region_markers_list = list()
    for (region in unique(subset_tp_markers$cluster)) {
        region_markers = subset_tp_markers[subset_tp_markers$cluster == region, ]$gene[1:100] ## all DE genes -> top 100 genes
        region_markers_list[[region]] = region_markers
    }
    tp_region_markers_list[[tp]] = region_markers_list
    venn.plot = venn.diagram(
      x = region_markers_list,
      category.names = names(region_markers_list), 
      filename = file.path(result_dir, paste0(tp, '_marker_Venn_diagram_spaDOT.png'))
    )
}

t0_region_markers_list = tp_region_markers_list[['t0']]
t1_region_markers_list = tp_region_markers_list[['t1']]
overlap_mat = matrix(0, nrow=3, ncol=3)
rownames(overlap_mat) = colnames(overlap_mat) = unique(subset_tp_markers$cluster)
for (region_i in unique(subset_tp_markers$cluster)) {
    for (region_j in unique(subset_tp_markers$cluster)) {
        common_genes = intersect(t0_region_markers_list[[region_i]], t1_region_markers_list[[region_j]])
        overlap_mat[region_i, region_j] = length(common_genes)
    }
}

# visualze common genes for injury_zone in t2
t2_Seurat_obj = subset(Seurat_obj, subset=timepoint=='t2')
t2_Seurat_obj = SCTransform(t2_Seurat_obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
common_genes = intersect(t0_region_markers_list[['injury_zone']], t1_region_markers_list[['injury_zone']])
cluster_counts = t2_Seurat_obj@assays$SCT@data[common_genes, ]
avg_counts = colMeans(cluster_counts)
t2_Seurat_obj$avg_counts = avg_counts
ggplot(t2_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red') + 
    theme_pubr()
ggsave(file.path(result_dir, paste0('t2', '_injury_zones_common_genes.png')))
 
common_genes = intersect(t0_region_markers_list[['myofiber']], t1_region_markers_list[['myofiber']])
cluster_counts = t2_Seurat_obj@assays$SCT@data[common_genes, ]
avg_counts = colMeans(cluster_counts)
t2_Seurat_obj$avg_counts = avg_counts
ggplot(t2_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red') + 
    theme_pubr()
ggsave(file.path(result_dir, paste0('t2', '_myofiber_common_genes.png')))
 
common_genes = intersect(t0_region_markers_list[['injury_border']], t1_region_markers_list[['injury_border']])
cluster_counts = t2_Seurat_obj@assays$SCT@data[common_genes, ]
avg_counts = colMeans(cluster_counts)
t2_Seurat_obj$avg_counts = avg_counts
ggplot(t2_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red') + 
    theme_pubr()
ggsave(file.path(result_dir, paste0('t2', '_injury_border_common_genes.png')))
 
# visualize in t1
t1_Seurat_obj = subset(Seurat_obj, subset=timepoint=='t1')
t1_Seurat_obj = SCTransform(t1_Seurat_obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
common_genes = intersect(t0_region_markers_list[['injury_zone']], t1_region_markers_list[['injury_zone']])
cluster_counts = t1_Seurat_obj@assays$SCT@data[common_genes, ]
avg_counts = colMeans(cluster_counts)
t1_Seurat_obj$avg_counts = avg_counts
ggplot(t1_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red') + 
    theme_pubr()
ggsave(file.path(result_dir, paste0('t1', '_injury_zones_common_genes.png')))
 
common_genes = intersect(t0_region_markers_list[['myofiber']], t1_region_markers_list[['myofiber']])
cluster_counts = t1_Seurat_obj@assays$SCT@data[common_genes, ]
avg_counts = colMeans(cluster_counts)
t1_Seurat_obj$avg_counts = avg_counts
ggplot(t1_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red') + 
    theme_pubr()
ggsave(file.path(result_dir, paste0('t1', '_myofiber_common_genes.png')))
 
common_genes = intersect(t0_region_markers_list[['injury_border']], t1_region_markers_list[['injury_border']])
cluster_counts = t1_Seurat_obj@assays$SCT@data[common_genes, ]
avg_counts = colMeans(cluster_counts)
t1_Seurat_obj$avg_counts = avg_counts
ggplot(t1_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red') + 
    theme_pubr()
ggsave(file.path(result_dir, paste0('t1', '_injury_border_common_genes.png')))
 



# --- top 100 SVGs
for (tp in c('t0', 't1', 't2')) {
    tp_SVGs = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx.csv')), header=T, row.names=1)
    tp_SVGs_cluster = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx_clustered.csv')), header=T, row.names=1)
    tp_Seurat_obj = subset(Seurat_obj, subset=section == unname(section_dict[tp]))
    tp_df = data.frame(matrix(nrow=0, ncol=4))
    colnames(tp_df) = c('timepoint', 'annotation', 'SVG_gene', 'enrichment')
    Idents(tp_Seurat_obj) = tp_Seurat_obj$annotation
    for (cluster in unique(tp_Seurat_obj$annotation)) {
        cluster_cells_idx = tp_Seurat_obj$annotation == cluster
        for (gene in rownames(tp_SVGs)[1:100]) {
            avg_exprs = mean(tp_Seurat_obj[gene, cluster_cells_idx]$RNA@layers$counts)
            other_exprs = mean(tp_Seurat_obj[gene, !cluster_cells_idx]$RNA@layers$counts)
            fc_enrichment = (avg_exprs-other_exprs)/other_exprs
            tp_df[nrow(tp_df) + 1,] = c(tp, cluster, gene, fc_enrichment)
        }
    }
    write.csv(tp_df, file.path(result_dir, paste0(tp, '_SVG_gene_enrichment.csv')), quote=F)
    tp_df$enrichment = as.numeric(tp_df$enrichment)
    tp_df$SVG_gene = factor(tp_df$SVG_gene, levels=rownames(tp_SVGs)[1:100])
    ggplot(data=tp_df, aes(x=SVG_gene, y=enrichment, group=annotation)) +
        geom_line(aes(color=annotation)) + geom_point() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(file.path(result_dir, paste0(tp, '_SVG_gene_enrichment.png')))
}

# --- all SVGs
#for (tp in c('t0', 't1', 't2')) {
#    tp_SVGs = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx.csv')), header=T, row.names=1)
#    tp_SVGs_cluster = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx_clustered.csv')), header=T, row.names=1)
#    tp_Seurat_obj = subset(Seurat_obj, subset=section == unname(section_dict[tp]))
#    tp_df = data.frame(matrix(nrow=0, ncol=4))
#    colnames(tp_df) = c('timepoint', 'annotation', 'SVG_gene', 'enrichment')
#    Idents(tp_Seurat_obj) = tp_Seurat_obj$annotation
#    for (cluster in unique(tp_Seurat_obj$annotation)) {
#        cluster_cells_idx = tp_Seurat_obj$annotation == cluster
#        for (gene in rownames(tp_SVGs)) {
#            avg_exprs = mean(tp_Seurat_obj[gene, cluster_cells_idx]$RNA@layers$counts)
#            other_exprs = mean(tp_Seurat_obj[gene, !cluster_cells_idx]$RNA@layers$counts)
#            fc_enrichment = (avg_exprs-other_exprs)/other_exprs
#            tp_df[nrow(tp_df) + 1,] = c(tp, cluster, gene, fc_enrichment)
#        }
#    }
#    write.csv(tp_df, file.path(result_dir, paste0(tp, '_SVG_gene_enrichment_all.csv')), quote=F)
#}

# check clustering results alignment with region information
for (tp in c('t0', 't1', 't2')) {
    tp_df = read.csv(file.path(result_dir, paste0(tp, '_SVG_gene_enrichment_all.csv')), header=T, row.names=1)
    tp_SVGs_cluster = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx_clustered_louvain.csv')), header=T, row.names=1)
    print(tp)
    print(table(tp_SVGs_cluster$cluster))
    tp_df$SVG_cluster = rep(tp_SVGs_cluster$cluster, length(unique(tp_df$annotation)))
    tp_df$SVG_cluster = as.factor(tp_df$SVG_cluster)
    ggplot(tp_df, aes(x=SVG_cluster, y=enrichment)) + geom_violin() + facet_wrap(~annotation)
    ggsave(file.path(result_dir, paste(tp, 'SVG_gene_enrichment_louvain_violin.png', sep='_')))
}


# --- top 10 metagene per cluster
for (tp in c('t0', 't1', 't2')) {
    tp_SVGs = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx.csv')), header=T, row.names=1)
    tp_SVGs_cluster = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx_clustered.csv')), header=T, row.names=1)
    tp_Seurat_obj = subset(Seurat_obj, subset=section == unname(section_dict[tp]))
    tp_Seurat_obj = SCTransform(tp_Seurat_obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    tp_SVGs$cluster = tp_SVGs_cluster$cluster
    # select top 10 genes per cluster
    tp_SVGs$gene = rownames(tp_SVGs)
    top_genes <- tp_SVGs %>%
          group_by(cluster) %>%
          top_n(10, wt = -adjustedPval)
    for (cluster in unique(top_genes$cluster)) {
        cluster_genes = top_genes[top_genes$cluster == cluster, ]$gene
        #cluster_counts = tp_Seurat_obj[cluster_genes]$RNA@layers$data
        ## standardize and average
        #if (is.null(dim(cluster_counts))) { # only one gene
        #    standardized_counts = t(scale(cluster_counts))
        #} else {
        #standardized_counts = t(scale(t(as.matrix(cluster_counts))))
        #}
        #avg_counts = colMeans(standardized_counts)
        cluster_counts = tp_Seurat_obj@assays$SCT@counts[cluster_genes, ]
        if (is.null(dim(cluster_counts))) {
            avg_counts = cluster_counts
        } else {
            avg_counts = colMeans(cluster_counts)
        }
        tp_Seurat_obj$avg_counts = avg_counts
        ggplot(tp_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red')
        ggsave(file.path(result_dir, paste0(cluster, '_', tp, '_aggregated_top10SVGs_per_cluster.png')))
    }
}
# --- all megagene per cluster
for (tp in c('t0', 't1', 't2')) {
    tp_SVGs = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx.csv')), header=T, row.names=1)
    tp_SVGs_cluster = read.csv(file.path(result_dir, paste0(unname(SVG_dict[tp]), '_SVG_sparkx_clustered_louvain.csv')), header=T, row.names=1)
    tp_Seurat_obj = subset(Seurat_obj, subset=section == unname(section_dict[tp]))
    tp_Seurat_obj = SCTransform(tp_Seurat_obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    tp_SVGs$cluster = tp_SVGs_cluster$cluster
    tp_SVGs$gene = rownames(tp_SVGs)
    for (cluster in unique(tp_SVGs$cluster)) {
        cluster_genes = tp_SVGs[tp_SVGs$cluster == cluster, ]$gene
        #cluster_counts = tp_Seurat_obj[cluster_genes]$RNA@layers$data
        # standardize and average
        #if (is.null(dim(cluster_counts))) { # only one gene
        #    standardized_counts = t(scale(cluster_counts))
        #} else {
        #    standardized_counts = t(scale(t(as.matrix(cluster_counts))))
        #}
        #avg_counts = colMeans(standardized_counts)
        cluster_counts = tp_Seurat_obj@assays$SCT@counts[cluster_genes, ]
        avg_counts = colMeans(cluster_counts)
        tp_Seurat_obj$avg_counts = avg_counts
        ggplot(tp_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=avg_counts)) + geom_point() + scale_colour_gradient(low='grey', high='red')
        ggsave(file.path(result_dir, paste0(cluster, '_', tp, '_aggregated_allSVGs_per_cluster_louvain.png')))
    }
}



