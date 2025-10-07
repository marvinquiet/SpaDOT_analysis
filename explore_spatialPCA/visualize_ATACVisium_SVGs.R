.libPaths("/net/mulan/home/wenjinma/Rlib")
library(SpatialPCA)
library(SeuratDisk)
library(Seurat)
#library(Signac)
library(aricode)
library(ggplot2)
library(rhdf5)
library(tidyverse)

prefix_dir = "/net/mulan/home/wenjinma"

# --- visualize ATAC-visium data
data_dir = file.path(prefix_dir, 'data/spatial_scRNAseq_GSE214989')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_GSE214989_multiplesample_combinedlabels')
dir.create(result_dir, showWarnings = F)

Seurat_obj = readRDS(file.path(data_dir, 'GSE214989_processed_SeuratObj.RDS'))
SVG_dict = c('t0'='1', 't1'='2', 't2'='3')
section_dict = c('t0'='3', 't1'='2', 't2'='1')

# --- plot muscle related genes
# https://pubmed.ncbi.nlm.nih.gov/26748348/ Expression of FGFs during early mouse tongue development
tongue_marker_dir = file.path(result_dir, 'tongue_markers')
dir.create(tongue_marker_dir, showWarnings = F)
for (tp in c('t0', 't1', 't2')) {
    tp_Seurat_obj = subset(Seurat_obj, subset=section == unname(section_dict[tp]))
    tp_Seurat_obj = SCTransform(tp_Seurat_obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    marker_genes = rownames(tp_Seurat_obj)[grep('Fgf', rownames(tp_Seurat_obj))]
    cluster_counts = tp_Seurat_obj@assays$SCT@counts[marker_genes, ]
    for (marker in marker_genes) {
        df = data.frame('pixel_x'=tp_Seurat_obj$pixel_x, 'pixel_y'=tp_Seurat_obj$pixel_y, 'norm_count'=cluster_counts[marker, ])
        ggplot(df, aes(x=pixel_x, y=pixel_y, color=norm_count)) + geom_point()
        ggsave(file.path(tongue_marker_dir, paste(tp, marker, '_gene_exprs.png', sep='_')))
    }
}

# Nuclear factor 1 X-type-associated regulation of myogenesis in developing mouse tongue
markers = c('Nfix', 'Mck', 'Myh8', 'Myh3', 'Myh7')
for (tp in c('t0', 't1', 't2')) {
    tp_Seurat_obj = subset(Seurat_obj, subset=section == unname(section_dict[tp]))
    tp_Seurat_obj = SCTransform(tp_Seurat_obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    for (marker in markers) {
        if (marker %in% rownames(tp_Seurat_obj@assays$SCT@counts)) {
            df = data.frame('pixel_x'=tp_Seurat_obj$pixel_x, 'pixel_y'=tp_Seurat_obj$pixel_y, 'norm_count'=tp_Seurat_obj@assays$SCT@counts[marker, ])
            ggplot(df, aes(x=pixel_x, y=pixel_y, color=norm_count)) + geom_point()
            ggsave(file.path(tongue_marker_dir, paste(tp, marker, '_gene_exprs.png', sep='_')))
        }
    }
}

# --- plot forebrain development related marker genes
forebrain_marker_dir = file.path(result_dir, 'forebrain')
dir.create(forebrain_marker_dir, showWarnings = F)
marker_genes = c('Pax6')
for (tp in c('t0', 't1', 't2')) {
    tp_Seurat_obj = subset(Seurat_obj, subset=section == unname(section_dict[tp]))
    tp_Seurat_obj = SCTransform(tp_Seurat_obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    cluster_counts = tp_Seurat_obj@assays$SCT@counts[marker_genes, ]
    for (marker in marker_genes) {
        df = data.frame('pixel_x'=tp_Seurat_obj$pixel_x, 'pixel_y'=tp_Seurat_obj$pixel_y, 'norm_count'=cluster_counts[marker, ])
        ggplot(df, aes(x=pixel_x, y=pixel_y, color=norm_count)) + geom_point()
        ggsave(file.path(forebrain_marker_dir, paste(tp, marker, '_gene_exprs.png', sep='_')))
    }
}






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



