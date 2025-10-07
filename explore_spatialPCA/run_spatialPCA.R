## conda activate ~/gpu/envs/R_env

library(SpatialPCA)
library(SeuratDisk)
library(Seurat)
library(Signac)
library(aricode)
library(ggplot2)

data_dir = "/net/mulan/home/wenjinma/data/spatial_scATACseq_GSE214991/"
result_dir = "/net/mulan/home/wenjinma/projects/spatialATAC-explore/results/SpatialPCA/"
Seurat_files = list.files(data_dir, pattern="*.h5Seurat")
for (Seurat_file in Seurat_files) {
    ## process Seurat file
    sample = gsub('.h5Seurat', '', Seurat_file)
    cat("Processing", sample, '\n')
    data = LoadH5Seurat(file.path(data_dir, Seurat_file))

    ## create SpatialPCA object
    xy_coords = as.matrix(data@meta.data[, c('x', 'y')])
    ## run SPARK to keep more spatially variable genes
    spark_obj = SPARK::CreateSPARKObject(counts=GetAssayData(data, 'counts'), location=as.data.frame(xy_coords), 
                                         percentage=0.01, min_total_counts=1) ## otherwise might have problem...
    spark_obj@lib_size = apply(spark_obj@counts, 2, sum)
    spark_obj = SPARK::spark.vc(spark_obj, 
                       covariates = NULL, 
                       lib_size = spark_obj@lib_size, 
                       num_core = 4,
                       verbose = F,
                       fit.model = 'gaussian')
    spark_obj = SPARK::spark.test(spark_obj, check_positive = T, verbose = F)
    significant_peak_idx = which(spark_obj@res_mtest$adjusted_pvalue <= 0.5)
    SVGnames = rownames(spark_obj@res_mtest[significant_peak_idx, ])
    cat("=== Number of features selected:", length(SVGnames), '\n')

    spatialPCA_obj = CreateSpatialPCAObject(counts=GetAssayData(data, 'counts'),
                                            location=xy_coords, project=sample,
                                            customGenelist=SVGnames, min.loctions=20, min.features=20)
    ## estimate spatial PCs
    spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ", bandwidth.set.by.user=NULL)
    ## ground truth information
    truth = data@meta.data$clusters[match(colnames(spatialPCA_obj@normalized_expr),colnames(GetAssayData(data, 'counts')))]
    plot_cluster(location=spatialPCA_obj@location, clusterlabel=truth, 
                 pointsize=1.5, text_size=20, 
                 title_in=paste("Ground Truth", sample),
                 color_in=rainbow(length(unique(data@meta.data$clusters))),
                 legend='bottom')
    ggsave(file.path(result_dir, paste0(sample, '_GroundTruth.png')))

    if (dim(spatialPCA_obj@params$expr)[1] < 20) {
        cat("! Error processing data because the feature space is too small as", dim(spatialPCA_obj@params$expr)[1], '\n')
    } else {
        spatialPCA_obj = SpatialPCA_EstimateLoading(spatialPCA_obj, fast=FALSE, SpatialPCnum=20) 
        spatialPCA_obj = SpatialPCA_SpatialPCs(spatialPCA_obj, fast=FALSE)
        clusterlabel = walktrap_clustering(clusternum=length(unique(data@meta.data$clusters)),
                                           latent_dat=spatialPCA_obj@SpatialPCs, 
                                           knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2])))
        clusterlabel_refine = refine_cluster_10x(clusterlabels=clusterlabel,location=spatialPCA_obj@location,shape="hexagon")
        plot_cluster(location=spatialPCA_obj@location, clusterlabel=clusterlabel_refine, 
                     pointsize=1.5, text_size=20, 
                     title_in=paste("SpatialPCA", sample),
                     color_in=rainbow(length(unique(data@meta.data$clusters))),
                     legend='bottom')
        ggsave(file.path(result_dir, paste0(sample, '_SpatialPCA.png')))
    
        cat(sample, 'ARI:', ARI(clusterlabel_refine, truth), '\n')
    }

    ## what if we first use TF-IDF normalization on the data and use it as input
    cat("After TF-IDF on", sample, '\n')
    TFIDF_data = RunTFIDF(GetAssayData(data, 'counts'))
    xy_coords = as.matrix(data@meta.data[, c('x', 'y')])
    spatialPCA_obj = CreateSpatialPCAObject(counts=TFIDF_data,
                                            location=xy_coords, project=sample,
                                            gene.type="spatial", sparkversion="spark", numCores_spark=4, 
                                            gene.number=10000, customGenelist=NULL, min.loctions = 20, min.features=20)
    spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ", bandwidth.set.by.user=NULL)
 
    if (dim(spatialPCA_obj@params$expr)[1] < 20) {
        cat("! Error processing data because the feature space is too small as", dim(spatialPCA_obj@params$expr)[1], '\n')
    } else {
        spatialPCA_obj = SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ", bandwidth.set.by.user=NULL)
        spatialPCA_obj = SpatialPCA_EstimateLoading(spatialPCA_obj, fast=FALSE, SpatialPCnum=20) 
        spatialPCA_obj = SpatialPCA_SpatialPCs(spatialPCA_obj, fast=FALSE)
        clusterlabel = walktrap_clustering(clusternum=length(unique(data@meta.data$clusters)),
                                           latent_dat=spatialPCA_obj@SpatialPCs, 
                                           knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2])))
        clusterlabel_refine = refine_cluster_10x(clusterlabels=clusterlabel,location=spatialPCA_obj@location,shape="hexagon")
        plot_cluster(location=spatialPCA_obj@location, clusterlabel=clusterlabel_refine, 
                     pointsize=1.5, text_size=20, 
                     title_in=paste("SpatialPCA TFIDF", sample),
                     color_in=rainbow(length(unique(data@meta.data$clusters))),
                     legend='bottom')
        ggsave(file.path(result_dir, paste0(sample, '_TFIDF_SpatialPCA.png')))
        truth = data@meta.data$clusters[match(colnames(spatialPCA_obj@normalized_expr), colnames(TFIDF_data))]
        cat('ARI (After TFIDF):', ARI(clusterlabel_refine, truth), '\n')
    }
}

