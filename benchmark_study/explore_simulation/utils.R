## Utilities functions for Splatter simulation

#' Simulate counts on LIBD with t0, t1, t2 timepoint
## === Settings:
## t0: 4 clusters, 1, 2, 3, 4
## t1: 5 clusters, 1, 2, 2.1, 3, 4
## t2: 7 clusters, 1, 2, 2.1, 3, 3.1, 3.2, 4
## Goal: Generate 7 clusters/groups of cells based on t2 point
## Method: Generate three batches of four cell types and subset accordingly
 
#' @param subsample
#' @param nGenes
#' @param nGroups
#' @param nBatches
#' @param seed 
#' @param ... splatter parameters

splatter_simu_LIBD = function(init_params, nCells, nGenes=5000, nGroups=4, nBatches=3, seed=1234,
			      out.prob=0, de.prop=0.5, de.facLoc=0.5, de.facScale=0.5,
			      batch.rmEffect=F, batch.facLoc=0.2, batch.facScale=0.2) {
    require(splatter)
    #N = nrow(subsample)
    #location = subsample[, 2:3]
    #label = subsample$label
    params = setParams(
        init_params,
        batchCells = rep(nCells*nGroups, nBatches),
        batch.rmEffect = batch.rmEffect,
        batch.facLoc = batch.facLoc,
        batch.facScale = batch.facScale,
        nGenes = nGenes,
        group.prob = rep(1/nGroups, nGroups),
        out.prob = out.prob,
        de.prob = de.prop,
        de.facLoc = de.facLoc,
        de.facScale = de.facScale,
        seed = seed)
    sim_groups = splatSimulate(
        params = params,
        method = "groups",
        verbose = FALSE)
    # subset certain batch and groups
    idx = which((colData(sim_groups)$Batch == 'Batch1' & colData(sim_groups)$Group %in% c('Group1', 'Group4')) |
            (colData(sim_groups)$Group == 'Group2' & colData(sim_groups)$Batch %in% c('Batch1', 'Batch2')) |
            (colData(sim_groups)$Group == 'Group3'))
    sub_sim_groups = sim_groups[, idx]
    return(sub_sim_groups)
}


#' Simulate counts on Brain
## === Settings:
# t0: 3 regions, 4 cell types; For neurons/GABA, 20% DEGs between regions; other two, no DEGs between regions
# For timepoint, when there is a region split, there is abatch effect
# --- For excitatory neurons:
# t0-forebrain: G1B1, t1-forebrain1: G1B2, t1-forebrain2: G1B3, t2-forebrain1: G1B2, t2-forebrain2: G1B3
# t0-midbrain: G2B1, t1-midbrain: G2B1, t2-midbrin: G2B1
# t0-hindbrain: G3B1, t1-hindbrain: G3B1, t2-hindbrain1: G3B2, t2-hindbrain2: G3B3
# --- For Astrocytes:
# t0-forebrain: G1B1, t1-forebrain1: G1B2, t1-forebrain2: G1B3, t2-forebrain1: G1B2, t2-forebrain2: G1B3
# t0-midbrain: G1B1, t1-midbrain: G1B1, t2-midbrin: G1B1
# t0-hindbrain: G1B1, t1-hindbrain: G1B1, t2-hindbrain1: G1B4, t2-hindbrain2: G1B5

splatter_simu_brain = function(init_params, nCells, nGenes=5000, nGroups=4, nBatches=3, seed=1234,
			      out.prob=0, de.prop=0.5, de.facLoc=0.5, de.facScale=0.5,
			      batch.rmEffect=F, batch.facLoc=0.2, batch.facScale=0.2) {
    require(splatter)
    #N = nrow(subsample)
    #location = subsample[, 2:3]
    #label = subsample$label
    params = setParams(
        init_params,
        batchCells = rep(nCells*nGroups, nBatches),
        batch.rmEffect = batch.rmEffect,
        batch.facLoc = batch.facLoc,
        batch.facScale = batch.facScale,
        nGenes = nGenes,
        group.prob = rep(1/nGroups, nGroups),
        out.prob = out.prob,
        de.prob = de.prop,
        de.facLoc = de.facLoc,
        de.facScale = de.facScale,
        seed = seed)
    sim_groups = splatSimulate(
        params = params,
        method = "groups",
        verbose = FALSE)
    return(sim_groups)
}


splatter_simu_brain_prop = function(init_params, nCells, nGenes=5000, nGroups=3, seed=1234,
			      out.prob=0, de.prop=0.5, de.facLoc=0.5, de.facScale=0.5,
			      batch.rmEffect=F, batch.facLoc=0, batch.facScale=0) {
    require(splatter)
    #N = nrow(subsample)
    #location = subsample[, 2:3]
    #label = subsample$label
    params = setParams(
        init_params,
        batchCells = rep(nCells*nGroups, 1),
        batch.rmEffect = batch.rmEffect,
        batch.facLoc = batch.facLoc,
        batch.facScale = batch.facScale,
        nGenes = nGenes,
        group.prob = rep(1/nGroups, nGroups),
        out.prob = out.prob,
        de.prob = de.prop,
        de.facLoc = de.facLoc,
        de.facScale = de.facScale,
        seed = seed)
    sim_groups = splatSimulate(
        params = params,
        method = "groups",
        verbose = FALSE)
    return(sim_groups)
}


#' Generate coordinates from an image
#' @param image_file filepath to drawn image
#' @param num_clusters pre-determined number of clusters in the image
#' @param sample_size number of locations to be subsampled from the image
image_simu_LIBD = function(image_file, num_clusters, sample_size=20000, seed=1234, Knum_start=1) {
    require(jpeg)
    require(grid)
    require(gridExtra)
    mandrill = readJPEG(image_file, native = FALSE)

    ## generate dataframe for RGB channel
    df = data.frame(
      red = matrix(mandrill[,,1], ncol=1),
      green = matrix(mandrill[,,2], ncol=1),
      blue = matrix(mandrill[,,3], ncol=1)
    )
    ### compute the k-means clustering, to obtain regions from image by color
    set.seed(seed)
    K = kmeans(df, num_clusters+1, nstart=Knum_start, iter.max=100) # one more region for white background
    df$label = K$cluster
    print('Label distribution:')
    print(table(df$label))
    ### Replace the color of each pixel in the image with the mean 
    ### R,G, and B values of the cluster in which the pixel resides:
    # get the coloring
    colors = data.frame(
        label = 1:nrow(K$centers), 
        R = K$centers[,"red"],
        G = K$centers[,"green"],
        B = K$centers[,"blue"]
    )
    # merge color codes on to df
    # IMPORTANT: we must maintain the original order of the df after the merge!
    df$order = 1:nrow(df)
    df = merge(df, colors)
    df = df[order(df$order),]
    df$order = NULL
    # get mean color channel values for each row of the df.
    R = matrix(df$R, nrow=dim(mandrill)[1])
    G = matrix(df$G, nrow=dim(mandrill)[1])
    B = matrix(df$B, nrow=dim(mandrill)[1])
    # reconstitute the segmented image in the same shape as the input image
    mandrill.segmented = array(dim=dim(mandrill))
    mandrill.segmented[,,1] = R
    mandrill.segmented[,,2] = G
    mandrill.segmented[,,3] = B
    print('Plot reconstructed figure:')
    jpeg(file.path(result_dir, paste0('reconstructed_', basename(image_file))))
    grid.raster(mandrill.segmented)
    dev.off()
    RGB_label = R*B*G
    pixel_ind = c(1:dim(df)[1])
    x_coor = rep(1:dim(mandrill.segmented)[2], each=dim(mandrill.segmented)[1])
    y_coor = -rep(1:dim(mandrill.segmented)[1], times=dim(mandrill.segmented)[2])
    data_groundtruth = data.frame(pixel_ind, x_coor, y_coor)
    data_groundtruth$label = as.integer(as.factor(c(RGB_label))) # largest one is white background because R=1, G=1, B=1
    # remove white background 
    data_groundtruth_use = data_groundtruth[-which(data_groundtruth$label==(num_clusters+1)),]

    set.seed(seed)
    # sample 10000 locations from the data
    subsample = data_groundtruth_use[sample(1:dim(data_groundtruth_use)[1], sample_size, replace=F),]
    min_ycoord = min(subsample$y_coor)
    subsample$y_coor = subsample$y_coor-min_ycoord # make coordinates >=0
    return(list(gt=data_groundtruth, gt_in_use=data_groundtruth_use, subsample=subsample))
}


#' Estimate parameters from a dataset
#' @param data_dir data directory
splatter_estimate_params = function(data_dir, cell_type='Neurons') {
    require(Matrix)
    require(readr)
    require(openxlsx)
    require(splatter)
    PFC_ref = read.table(file.path(data_dir, "GSE104276_all_pfc_2394_UMI_count_NOERCC.xls.gz"))
    PFC_ref = as.matrix(PFC_ref)
    celltype_ref = read.xlsx(file.path(data_dir, "GSE104276_readme_sample_barcode.xlsx"),5)
    PFC_ref_use = PFC_ref[,match(as.character(celltype_ref[,1]), colnames(PFC_ref))]
    raw.data = PFC_ref_use
    cell_types = celltype_ref$cell_types
    names(cell_types) = as.character(celltype_ref[,1]) # create cell_types named list
    cell_types = as.factor(cell_types) # convert to factor data type
    cnts = as.matrix(raw.data)
    init_params = splatEstimate(cnts[,which(cell_types==cell_type)])
    return(init_params)
}

#' Permute Z to change cell type proportion for each region
#' I do not have a better idea to generate this permutation now
#' @param z region label
map_z2c = function(z)
{
  dplyr::case_when(
    # first index has the first corresponding proportion
    z == "1" ~ c("1", "2", "3", "4"),
    z == "2" ~ c("2", "3", "4", "1"),
    z == "3" ~ c("3", "4", "1", "2"),
    z == "4" ~ c("4", "1", "2", "3"),
    z == "5" ~ c("5", "3", "4", "1"),
    z == "6" ~ c("6", "4", "1", "2"),
    z == "7" ~ c("7", "4", "1", "2")
    )
}

#' Map scenarios to proportion
#' @param scenario different scenarios
map_scenario2prop = function(scenario) {
    dplyr::case_when(
        scenario == 'scenario2' ~ c(0.6, 0.3, 0.05, 0.05),
        scenario == 'scenario4' ~ c(0.35, 0.30, 0.30, 0.05)
    )
}



#  1. Set parameters as to discrimate the groups and batches, 
#  where 2 and 2.1 has batch effect and 3, 3.1 and 3.2 also 
#  has batch effect
#  2. Generate 1, 4 together as two groups;
#  3. Generate 2, 2.1 as three group with fewer DE genes;
#  4. Generate 3, 3.1 and 3.2 as three group with fewer DE genes;
#  5. Concat them all together
 
splatter_simu_LIBD_failed_version = function() {
### simulate 1 and 4
#group.prob = c(0.5, 0.5)  # same probability for 7 clusters
#out.prob = 0
#de.prop = de.facLoc = de.facScale = 0.5
#batch.rmEffect = TRUE # no batch
#batch.facLoc = 0 # no batch factor location
#params = setParams(
#    init_params,
#    batchCells = rep(N*2, 1),
#    batch.rmEffect = batch.rmEffect,
#    batch.facLoc = batch.facLoc,
#    nGenes = J,
#    group.prob = group.prob,
#    out.prob = out.prob,
#    de.prob = de.prop,
#    de.facLoc = de.facLoc,
#    de.facScale = de.facScale,
#    seed = sim_seed)
#sim_1_and_4_groups = splatSimulate(
#    params = params,
#    method = "groups",
#    verbose = FALSE)
#levels(colData(sim_1_and_4_groups)$Group) = c('Group1', 'Group4')
#
### simulate 2 and 2.1
#group.prob = c(0.5, 0.5)  # same probability for 7 clusters
#out.prob = 0
#de.prop = de.facLoc = de.facScale = 0.1
#batch.rmEffect = TRUE # no batch
#batch.facLoc = 0 # no batch factor location
#params = setParams(
#    init_params,
#    batchCells = rep(N*2, 1),
#    batch.rmEffect = batch.rmEffect,
#    batch.facLoc = batch.facLoc,
#    nGenes = J,
#    group.prob = group.prob,
#    out.prob = out.prob,
#    de.prob = de.prop,
#    de.facLoc = de.facLoc,
#    de.facScale = de.facScale,
#    seed = sim_seed)
#sim_2_groups = splatSimulate(
#    params = params,
#    method = "groups",
#    verbose = FALSE)
#levels(colData(sim_2_groups)$Group) = c('Group2', 'Group2.1')
#
### simulate 3, 3.1 and 3.2
#group.prob = rep(1/3, 3)  # same probability for 7 clusters
#out.prob = 0
#de.prop = de.facLoc = de.facScale = 0.1
#batch.rmEffect = TRUE # no batch
#batch.facLoc = 0 # no batch factor location
#params = setParams(
#    init_params,
#    batchCells = rep(N*2, 1),
#    batch.rmEffect = batch.rmEffect,
#    batch.facLoc = batch.facLoc,
#    nGenes = J,
#    group.prob = group.prob,
#    out.prob = out.prob,
#    de.prob = de.prop,
#    de.facLoc = de.facLoc,
#    de.facScale = de.facScale,
#    seed = sim_seed)
#sim_3_groups = splatSimulate(
#    params = params,
#    method = "groups",
#    verbose = FALSE)
#levels(colData(sim_3_groups)$Group) = c('Group3', 'Group3.1', 'Group3.2')

## concat three objects together
# directly use cbind could cause problem

}
