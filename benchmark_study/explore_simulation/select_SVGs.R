## conda activate /net/mulan/home/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(Seurat)
library(optparse)
library(SPARK)
library(tidyverse)
library(Matrix)
library(FNN)
library(igraph)
library(bluster)

# set optparse
option_list = list(
    make_option(c("--scenario_dir"), type='character', default=NULL, help="Simulation directory.")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)
scenario_dir = opt$scenario_dir
# load data
tps = c('t0', 't1', 't2')
for (tp in tps) {
    tp_data = readMM(file.path(scenario_dir, paste0(tp, '.mtx')))
    tp_metadata = read.csv(file.path(scenario_dir,paste0(tp, '_metadata.csv')), header=T, row.names=1)
    obj = CreateSeuratObject(tp_data, meta.data=tp_metadata)
    obj = SCTransform(obj, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.3)
    count_spark = obj@assays$RNA@layers$counts
    rownames(count_spark) = rownames(obj)
    colnames(count_spark) = colnames(obj)
    locations_spark = cbind(obj$x_coor, obj$y_coor)
    locations_spark = as.matrix(locations_spark)
    sparkX = sparkx(count_spark, locations_spark,numCores=5)
    significant_gene_number = sum(sparkX$res_mtest$adjustedPval<=0.05)
    SVGs = sparkX$res_mtest[order(sparkX$res_mtest$adjustedPval),][1:significant_gene_number, ]
    write.csv(SVGs, file.path(scenario_dir, paste0(tp, '_SVG_sparkx.csv')), quote=F)
    SVG_mat = obj@assays$SCT@scale.data[rownames(SVGs), ]
    SVG_pcs = prcomp(SVG_mat, center=T, scale=T)
    SVG_top_pcs = SVG_pcs$x[, 1:30]
    knn.norm = FNN::get.knn(as.matrix(SVG_top_pcs), k=30)
    knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
               k=30), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
    nw.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
    nw.norm = igraph::simplify(nw.norm)
    lc.norm = igraph::cluster_louvain(nw.norm)
    merged = bluster::mergeCommunities(nw.norm, lc.norm$membership, number=10)
    clusterlabel = as.character(as.integer(as.factor(paste0("cluster",merged))))
    SVGs$cluster = clusterlabel
    write.csv(SVGs, file.path(scenario_dir, paste0(tp, '_SVG_sparkx_clustered_louvain.csv')), quote=F)
}

