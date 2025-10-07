## conda activate /media/sda/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(SeuratDisk)
library(Seurat)
#library(Signac)
library(aricode)
library(ggplot2)
library(ggpubr)
library(rhdf5)
library(VennDiagram)

library(fgsea)
library(msigdbr)
C5_gene_sets = msigdbr(species = "mouse", category = "C5", subcategory = 'GO:BP') # subcategory=
msigdbr_list = split(x = C5_gene_sets$gene_symbol, f = C5_gene_sets$gs_name)

prefix_dir = "/net/mulan/home/wenjinma"
data_dir = file.path(prefix_dir, 'data/spatial_scRNAseq_GSE214989')
result_dir = file.path(prefix_dir, 'projects/spatialATAC-explore/results/SpatialPCA_GSE214989_multiplesample_combinedlabels')
dir.create(result_dir, showWarnings = F)

# load data
h5_file = file.path(data_dir, "GSE214989_counts_embryo_visium.h5")
data = Read10X_h5(h5_file)
metadata = read.csv(file.path(data_dir, 'GSE214989_metadata_embryo_visium.csv.gz'),
                    header=T, row.names=1, sep=';')
common_cells = intersect(colnames(data), rownames(metadata))
data = data[, common_cells]
metadata = metadata[common_cells, ] # 8255 cells

# perform clustering
Seurat_obj = CreateSeuratObject(counts=data, meta.data=metadata)
Seurat_obj = NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_obj = FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000)
Seurat_obj = ScaleData(Seurat_obj, features = rownames(Seurat_obj))
Seurat_obj = RunPCA(Seurat_obj, features = VariableFeatures(object = Seurat_obj))
Seurat_obj = RunUMAP(Seurat_obj, dims = 1:20)
# from data analysis code
Seurat_obj$annotation = ifelse(Seurat_obj$seurat_clusters_harmony %in% c("3", "2","4"), 
                            "CNS",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% c("1", "9", "11"),
                            "Men/PNS",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% c("0","5","6", "8"),
                            "Mesenchyme",
                            ifelse(Seurat_obj$seurat_clusters_harmony %in% "7",
                            "Liver", "Heart"))))
# load spaDOT results
spaDOT_res = NULL
for (timepoint in c('t0', 't1', 't2')) {
    # tp_res = read.csv(file.path(result_dir, 'fgsea_for_spaDOT_labels_afterremoval', paste0(timepoint, '_Louvain_labels.csv')), header=T, row.names=1)
    tp_res = read.csv(file.path(result_dir, 'fgsea_for_spaDOT_KMeans_labels', paste0(timepoint, '_KMeans_labels.csv')), header=T, row.names=1)
    if (is.null(spaDOT_res)) {
        spaDOT_res = tp_res
    } else {
        spaDOT_res = rbind(spaDOT_res, tp_res)
    }
}
common_cells = intersect(colnames(Seurat_obj), rownames(spaDOT_res))
Seurat_obj = Seurat_obj[,common_cells]
spaDOT_res = spaDOT_res[common_cells, ]
# Seurat_obj@meta.data['spaDOT_res'] = spaDOT_res[rownames(Seurat_obj@meta.data), 'louvain']
Seurat_obj@meta.data['spaDOT_res'] = spaDOT_res[rownames(Seurat_obj@meta.data), 'kmeans']

Seurat_obj@meta.data$timepoint = ''
t1_idx = which(Seurat_obj@meta.data$section == 3)
Seurat_obj@meta.data[t1_idx, 'timepoint'] = 't0'
t2_idx = which(Seurat_obj@meta.data$section == 2)
Seurat_obj@meta.data[t2_idx, 'timepoint'] = 't1'
t3_idx = which(Seurat_obj@meta.data$section == 1)
Seurat_obj@meta.data[t3_idx, 'timepoint'] = 't2'
## find markers
for (tp in c('t0', 't1', 't2')) {
    tp_Seurat_obj = subset(Seurat_obj, subset=timepoint==tp)
    Idents(tp_Seurat_obj) = tp_Seurat_obj$spaDOT_res
    tp_markers = FindAllMarkers(tp_Seurat_obj)
    write.csv(tp_markers, file.path(result_dir, paste0(tp, '_spaDOT_Seurat_FindAllMarkers.csv')), quote=F)
    subset_tp_markers = tp_markers[abs(tp_markers$avg_log2FC) > 0 & tp_markers$p_val_adj < 0.01, ]
    write.csv(subset_tp_markers, file.path(result_dir, paste0(tp, '_spaDOT_Seurat_FindAllMarkers_log2FC_0_adjp_0.01.csv')), quote=F)

    # subset_tp_markers = read.csv(file.path(result_dir, paste0(tp, '_spaDOT_Seurat_FindAllMarkers_log2FC_0_adjp_0.01.csv')), header=T, row.names=1)
    # subset_tp_markers = subset_tp_markers[subset_tp_markers$avg_log2FC > 0, ]
    DEG_list = list()
    for (cluster in unique(subset_tp_markers$cluster)) {
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
        ggsave(file.path(result_dir, paste0(tp, '_cluster', cluster, '_fgsea_enrichment.png')), width=10, height=5)
    }
    # venn diagram
    #venn.plot = venn.diagram(
    #  x = DEG_list,
    #  category.names = names(DEG_list), 
    #  filename = file.path(result_dir, paste0(tp, '_DEG_Venn_diagram_spaDOT.tiff'))
    #)
}

## calculate gene enrichment score
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
    Idents(tp_Seurat_obj) = tp_Seurat_obj$spaDOT_res
    for (cluster in unique(tp_Seurat_obj$spaDOT_res)) {
        cluster_cells_idx = tp_Seurat_obj$spaDOT_res == cluster
        for (gene in names(marker_genes)) {
            avg_exprs = mean(tp_Seurat_obj[gene, cluster_cells_idx]$RNA@layers$counts)
            other_exprs = mean(tp_Seurat_obj[gene, !cluster_cells_idx]$RNA@layers$counts)
            fc_enrichment = (avg_exprs-other_exprs)/other_exprs
            df[nrow(df) + 1,] = c(tp, cluster, gene, unname(marker_genes[gene]), fc_enrichment)
        }
    }
}
write.csv(df, file.path(result_dir, 'fgsea_for_spaDOT_labels_afterremoval', 'cluster_gene_enrichment.csv'), quote=F)

df$enrichment = as.numeric(df$enrichment)
df$gene_region = paste0(df$gene, '-', df$region)
df$gene_region = as.factor(df$gene_region)
df$gene_region = factor(df$gene_region, 
                        levels=paste0(names(marker_genes),'-',unname(marker_genes)))
ggplot(data=df, aes(x=gene_region, y=enrichment, group=cluster)) + 
    geom_line(aes(color=cluster)) + geom_point() + 
    facet_wrap(~timepoint, ncol=1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file.path(result_dir, 'fgsea_for_spaDOT_labels_afterremoval', 'cluster_gene_enrichment.png'))
 

