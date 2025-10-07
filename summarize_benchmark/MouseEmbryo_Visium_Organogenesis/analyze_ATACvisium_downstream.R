Sys.setenv(RETICULATE_PYTHON = "/media/sda/wenjinma/envs/spatialATAC/bin/python")
library(reticulate)

library(VennDiagram)
library(anndata)
library(aricode)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(Seurat)
library(fgsea)
library(msigdbr)
C5_gene_sets = msigdbr(species = "mouse", category = "C5", subcategory = 'GO:BP') # subcategory=
msigdbr_list = split(x = C5_gene_sets$gene_symbol, f = C5_gene_sets$gs_name)

base_dir = "/media/sda/wenjinma"
result_dir = file.path(base_dir, "projects/spatialATAC-explore")
downstream_dir = file.path(result_dir, "summarize_benchmark/ATAC-Visium_figures/ATAC-Visium_downstream")
dir.create(downstream_dir, showWarnings = FALSE, recursive = TRUE)
latent_adata = read_h5ad(file.path(result_dir, 'explore_GATVAE-SVGP-VAE', 'ATAC-Visium_results', 'prediction_adata.h5ad'))

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

# paper marker genes for calculating enrichment score
marker_genes = c('Sox1'='CNS', 'Foxg1'='CNS', 'Notch1'='CNS', 'Fyn'='CNS',
                 'Sox2'='CNS', 'Lhx2'='CNS', 'Emx1'='CNS', 'Neurog1'='CNS',
                 'Cux2'='CNS', 'Neurod1'='CNS', 'Ascl1'='CNS', # CNS
                'Myh9'='Mesenchyme', 'Col1a1'='Mesenchyme', 'Smad3'='Mesenchyme', # Mesenchyme
                'Hba-a1'='Liver', 'Tal1'='Liver', 'Sptb'='Liver',
                'Gata4'='Liver', 'Gata6'='Liver', # Liver
                'Smad3'='BoneMuscle', 'Twist1'='BoneMuscle', 'Myog'='BoneMuscle') # bone and muscle
df = data.frame(matrix(nrow=0, ncol=5))
colnames(df) = c('timepoint', 'cluster', 'gene', 'region', 'enrichment')
for (tp in c('t0', 't1', 't2')) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    Idents(tp_Seurat_obj) = tp_Seurat_obj$pred_labels
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
    tp_markers = FindAllMarkers(tp_Seurat_obj)
    write.csv(tp_markers, file.path(downstream_dir, paste0(tp, '_Seurat_FindAllMarkers_original_cluster.csv')), quote=F)
    subset_tp_markers = tp_markers[tp_markers$avg_log2FC > 1 & tp_markers$p_val_adj < 0.05, ]
    write.csv(subset_tp_markers, file.path(downstream_dir, paste0(tp, '_Seurat_FindAllMarkers_log2FC_0_adjp_0.01_original_cluster.csv')), quote=F)
    # gene enrichment analysis
    DEG_list = list()
    for (cluster in levels(Idents(tp_Seurat_obj))) {
        cluster_tp_markers = subset_tp_markers[subset_tp_markers$cluster == cluster, ]
        DEG_list[[as.character(cluster)]] = cluster_tp_markers$gene
        DEG_ranks = sign(cluster_tp_markers$avg_log2FC) * -log10(cluster_tp_markers$p_val_adj)
        names(DEG_ranks) = cluster_tp_markers$gene
        inf_idx = which(is.infinite(DEG_ranks))
        if (length(inf_idx) > 0) {
            DEG_ranks[inf_idx] = max(DEG_ranks[!is.infinite(DEG_ranks)])
        }
        DEG_ranks = DEG_ranks[order(-DEG_ranks)]
        fgsea_res = fgsea(pathways = msigdbr_list, stats = DEG_ranks, minSize=15, maxSize=500)
        # collapse pathway
        #collapsedPathways = collapsePathways(fgsea_res[order(pval)][padj < 0.01],
        #        msigdbr_list, DEG_ranks)
        #mainPathways_res = fgsea_res[pathway %in% collapsedPathways$mainPathways][order(padj), ]
        fgsea_res = fgsea_res[order(padj), ]
        g_data = fgsea_res[1:10]
        #g_data = mainPathways_res[1:10]
        g_data$pathway = gsub("GOBP_", "", g_data$pathway)
        g_data$pathway = gsub("_", " ", g_data$pathway)
        g_data$pathway = factor(g_data$pathway, levels=rev(g_data$pathway))
        ggplot(g_data, aes(y=pathway, x=-log10(padj), fill=-log10(padj))) +
            geom_bar(stat = "identity") +
            theme_pubr()
        ggsave(file.path(downstream_dir, paste0(tp, '_cluster', cluster, '_fgsea_enrichment_original_cluster.png')), width=10, height=5)
    }
}
write.csv(df, file.path(downstream_dir, 'cluster_gene_enrichment_original_cluster.csv'), quote=F)
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

# Adaptive Kmeans clustering

# for (tp in c('t0', 't1', 't2')) {
#     tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
#     tp_latent = t(tp_Seurat_obj@assays$spaDOT_latent$data)
#     # --- what about using MClust then
#     BIC = mclustBIC(tp_latent)
#     mod = Mclust(tp_latent, x = BIC)
#     tp_Seurat_obj$mclust = as.factor(mod$classification)
#     ggplot(tp_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=mclust)) + geom_point()
#     ggsave(file.path(downstream_dir, paste0(tp, '_mclust_clustering.png')))
#     # --- what about Silhouette score
#     pam_res = pamk(tp_latent, krange=2:20, usepam=T, scaling=F)
#     tp_Seurat_obj$pam = as.factor(pam_res$pamobject$clustering)
#     ggplot(tp_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=pam)) + geom_point()
#     ggsave(file.path(downstream_dir, paste0(tp, '_pam_clustering.png')))

#     asw = lapply(2:20, function(k) {
#         require(fpc)
#         # k_model = kmeans(tp_latent, centers=k, nstart=20)
#         pamk(tp_latent, k)$silinfo$avg.width
#     })
#     # Compute BIC for each model
#     bics <- sapply(kmeans_models, function(model) {
#         n <- nrow(tp_latent)
#         k <- length(model$centers)
#         withinss <- sum(model$withinss)
#         BIC <- n * log(withinss / n) + k * log(n)
#         return(BIC)
#     })
#     # Find the optimal number of clusters that minimizes BIC
#     optimal_k = k_values[which.min(bics)]
#     # kmeans clustering
#     kmeans_res = kmeans(tp_latent, centers=5, nstart=20)
#     tp_Seurat_obj$kmeans = as.factor(kmeans_res$cluster)
#     cat(tp, 'KMeans ARI:', ARI(tp_Seurat_obj$annotation, tp_Seurat_obj$kmeans), '\n')
#     ggplot(tp_Seurat_obj@meta.data, aes(x=pixel_x, y=pixel_y, color=kmeans)) + geom_point()
#     ggsave(file.path(downstream_dir, paste0(tp, '_kmeans_clustering.png')))
# }






