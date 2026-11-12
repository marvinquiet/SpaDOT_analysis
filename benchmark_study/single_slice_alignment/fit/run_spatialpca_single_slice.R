#!/usr/bin/env Rscript
# Independent raw-slice SpatialPCA fits based on explore_spatialPCA's native steps.
# Cross-slice Seurat integration is intentionally omitted for a single-slice fit.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("usage: run_spatialpca_single_slice.R bridge_dir output_dir K_csv")
bridge <- args[[1]]
out <- args[[2]]
kvals <- as.integer(strsplit(args[[3]], ",", fixed=TRUE)[[1]])
if (dir.exists(out)) stop("output root collision")


suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SpatialPCA))
set.seed(2024)

timepoints <- sort(list.dirs(bridge, recursive=FALSE, full.names=FALSE))
if (length(timepoints) != length(kvals)) stop("timepoint/K mismatch")
dir.create(out, recursive=TRUE)

for (i in seq_along(timepoints)) {
    tp <- timepoints[[i]]
    block <- file.path(out, tp)
    dir.create(block)
    current <- file.path(bridge, tp)
    genes <- readLines(file.path(current, "genes.txt"))
    spots <- readLines(file.path(current, "spots.txt"))
    coords <- read.csv(file.path(current, "coordinates.csv"), stringsAsFactors=FALSE)
    if (!identical(coords$spot, spots)) stop("spot/coordinate order mismatch")
    counts <- Matrix::readMM(file.path(current, "counts_genes_by_spots.mtx"))
    if (!identical(dim(counts), c(length(genes), length(spots)))) stop("matrix dimensions mismatch")
    counts <- as(counts, "dgCMatrix")
    rownames(counts) <- genes
    colnames(counts) <- spots
    location <- as.matrix(coords[, c("pixel_x", "pixel_y")])
    rownames(location) <- spots

    obj <- SpatialPCA::CreateSpatialPCAObject(
        counts=counts, location=location, project="SpatialPCA",
        gene.type="spatial", sparkversion="sparkx", numCores_spark=5,
        gene.number=3000, customGenelist=NULL,
        min.loctions=0, min.features=0)
    selected <- rownames(obj@normalized_expr)
    if (length(selected) == 0 || anyDuplicated(selected)) stop("native selected genes missing")
    obj@location <- scale(obj@location)
    obj <- SpatialPCA::SpatialPCA_buildKernel(
        obj, kerneltype="gaussian", bandwidth.set.by.user=0.1)
    obj <- SpatialPCA::SpatialPCA_EstimateLoading(
        obj, fast=TRUE, SpatialPCnum=20)
    obj <- SpatialPCA::SpatialPCA_SpatialPCs(obj, fast=TRUE)
    latent <- t(as.matrix(obj@SpatialPCs))
    fitted_spots <- rownames(obj@location)
    if (is.null(fitted_spots)) fitted_spots <- spots
    if (!all(fitted_spots %in% spots) || anyDuplicated(fitted_spots)) stop("native spot filter mismatch")
    if (nrow(latent) != length(fitted_spots) || any(!is.finite(latent))) stop("invalid native latent")
    rownames(latent) <- fitted_spots
    predicted <- stats::kmeans(latent, centers=kvals[[i]])$cluster
    if (length(unique(predicted)) != kvals[[i]]) stop("invalid cluster count")
    write.csv(data.frame(spot=fitted_spots, predicted=as.vector(predicted)),
              file.path(block, "predicted.csv"), row.names=FALSE)
    write.csv(data.frame(spot=fitted_spots, latent),
              file.path(block, "latent.csv"), row.names=FALSE)
    writeLines(selected, file.path(block, "selected_genes.txt"))
    saveRDS(obj, file.path(block, "native_spatialpca.rds"))
    cat(sprintf("COMPLETE SpatialPCA %s spots=%d native_genes=%d K=%d\n",
                tp, length(fitted_spots), length(selected), kvals[[i]]))
    flush.console()
}
writeLines("COMPLETE", file.path(out, "status.txt"))
