#' simulate spatial RNA-seq data based on LIBD
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(splatter)
library(scater)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggpubr)
library(HDF5Array)
library(optparse)

## on Mulan
base_dir = '/net/mulan/home/wenjinma'
# project_dir = file.path(base_dir, "projects/spatialATAC-explore")
# work_dir = file.path(project_dir, 'explore_spatialPCA_simulation')
work_dir = '/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_simulation'
setwd(work_dir)
source('utils.R')

# set optparse
option_list = list(
    make_option(c("-x", "--x_DEG"), type="integer", default=20, help="x% DEGs when having domain split"),
    make_option(c("-y", "--y_DEG"), type="integer", default=50, help="x% DEGs among domains"),
    make_option(c("-s", "--seed"), type="integer", default=2024, help="random seed")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)


data_dir = '/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain'
# data_dir = file.path(base_dir, "data/spatial_simulation_brain")
# result_dir = file.path(project_dir, 'results', 'spatialPCA_simulation_brain')
# dir.create(result_dir, showWarnings = FALSE)
# after the generation, we move the data to /net/mulan/home/wenjinma/data/spatial_simulation_brain
result_dir = data_dir

seed = opt$seed
x_DEG = opt$x_DEG
y_DEG = opt$y_DEG
seed_dir = file.path(result_dir, seed)
dir.create(seed_dir, showWarnings=F)
scenario_dir = file.path(result_dir, seed, paste(x_DEG, y_DEG, sep='_'))
dir.create(scenario_dir, showWarnings=F)

#' estimate for once and use it later to do simulation
# NPCs_init_params = splatter_estimate_params('/net/zootopia/disk1/wenjinma/data/spatial_simulation', cell_type='Stem cells')
# saveRDS(NPCs_init_params, file=file.path(result_dir, "NPCs_init_params_pFC_for_LIBD.RDS"))
NPCs_init_params = readRDS(file.path(result_dir, "NPCs_init_params_pFC_for_LIBD.RDS"))
#Neurons_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='Neurons')
#saveRDS(Neurons_init_params, file = file.path(result_dir, "Neurons_init_params_pFC_for_LIBD.RDS"))
Neurons_init_params = readRDS(file.path(result_dir, "Neurons_init_params_pFC_for_LIBD.RDS"))
#GABANeurons_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='GABAergic neurons')
#saveRDS(GABANeurons_init_params, file = file.path(result_dir, "GABA_init_params_pFC_for_LIBD.RDS"))
GABANeurons_init_params = readRDS(file.path(result_dir, "GABA_init_params_pFC_for_LIBD.RDS"))
#Astrocytes_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='Astrocytes')
#saveRDS(Astrocytes_init_params, file = file.path(result_dir, "Astrocytes_init_params_pFC_for_LIBD.RDS"))
Astrocytes_init_params = readRDS(file.path(result_dir, "Astrocytes_init_params_pFC_for_LIBD.RDS"))
#OPC_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='OPC')
#saveRDS(OPC_init_params, file = file.path(result_dir, "OPC_init_params_pFC_for_LIBD.RDS"))
OPC_init_params = readRDS(file.path(result_dir, "OPC_init_params_pFC_for_LIBD.RDS"))

#' generate locations for image t0, t1, t2 
tp_clusters = c(3, 4, 5)
names(tp_clusters) = c('t0', 't1', 't2')
Knum_start = c(10, 10, 1)
names(Knum_start) = c('t0', 't1', 't2')
location_res = list()
for (timepoint in names(tp_clusters)) {
    #image_file = file.path(data_dir, paste0('brain_develop_', timepoint, '.jpg'))
    #num_clusters = unname(tp_clusters[timepoint])
    #res = image_simu_LIBD(image_file, num_clusters=num_clusters, sample_size=3000*num_clusters, seed=seed, 
    #                      Knum_start=unname(Knum_start[timepoint]))
    #subsample = res[['subsample']]
    #subsample$label = as.factor(subsample$label)
    ### test subsampled results
    #ggscatter(subsample, x='x_coor', y='y_coor', color='label', size=1)
    #subsample$label = as.factor(subsample$label)
    #ggsave(file.path(result_dir, paste0(timepoint, '_simulated_image.png')))
    # save results to list
    #saveRDS(res, file.path(result_dir, paste0(timepoint, '_simulated_image_data.RDS')))
    res = readRDS(file.path(result_dir, paste0(timepoint, '_simulated_image_data.RDS')))
    location_res[[timepoint]] = res
}
### find cluster correspondence among timepoints -> manually TODO: how to make it automatically corresponded?
#label_dict = list()
#label_dict[['t0']] = c("1"=1, "2"=2, "3"=3)
## label_dict[['t1']] = c("1"=3, "2"=1, "3"=4, "4"=2)
#label_dict[['t1']] = c("3"=1, "1"=2, "4"=3, "2"=4)
## label_dict[['t2']] = c("1"=4, "2"=1, "3"=2, "4"=3, "5"=5)
#label_dict[['t2']] = c("4"=1, "1"=2, "2"=3, "3"=4, "5"=5)
#color_palettes = RColorBrewer::brewer.pal(5, "Set1")
#names(color_palettes) = 1:5
#for (timepoint in names(tp_clusters)) {
#    subsample = location_res[[timepoint]][['subsample']]
#    subsample$curated_label = label_dict[[timepoint]][as.character(subsample$label)]
#    subsample$curated_label = as.factor(subsample$curated_label)
#    location_res[[timepoint]][['subsample']] = subsample
#    ggscatter(subsample, x='x_coor', y='y_coor', color='curated_label', palette=color_palettes)
#    ggsave(file.path(result_dir, paste0(timepoint, '_simulated_image_curated.png')))
#}

## generate simulated data for t0 with four cell types in each brain region
simulated_neurons = splatter_simu_brain(Neurons_init_params, nCells=5000, nGenes=5000, nGroups=3, nBatches=3, 
                              seed=seed, out.prob=0, de.prop=y_DEG/100, # de.facLoc=y_DEG/100, #de.facScale=y_DEG/100,
                              batch.rmEffect=F, batch.facLoc=x_DEG/100, #batch.facScale=x_DEG/100
                              )
colData(simulated_neurons)$cell_type = 'Neurons'
#saveHDF5SummarizedExperiment(simulated_neurons, scenario_dir, prefix="simulated_neurons_")
simulated_GABA = splatter_simu_brain(GABANeurons_init_params, nCells=5000, nGenes=5000, nGroups=3, nBatches=3,
                              seed=seed, out.prob=0, de.prop=y_DEG/100, #de.facLoc=y_DEG/100, #de.facScale=y_DEG/100,
                              batch.rmEffect=F, batch.facLoc=x_DEG/100, #batch.facScale=x_DEG/100
                              )
colData(simulated_GABA)$cell_type = 'GABA'
#saveHDF5SummarizedExperiment(simulated_GABA, scenario_dir, prefix="simulated_GABA_")
simulated_Astrocytes = splatter_simu_brain(Astrocytes_init_params, nCells=1000, nGenes=5000, nGroups=1, nBatches=5, 
                              seed=seed, out.prob=0, batch.rmEffect=F, batch.facLoc=x_DEG/100, #batch.facScale=x_DEG/100
                              )
colData(simulated_Astrocytes)$cell_type = 'Astrocytes'
colData(simulated_Astrocytes)$Group = 'Group1'
#saveHDF5SummarizedExperiment(simulated_Astrocytes, scenario_dir, prefix="simulated_Astrocytes_")
simulated_OPC = splatter_simu_brain(OPC_init_params, nCells=1000, nGenes=5000, nGroups=1, nBatches=5, 
                              seed=seed, out.prob=0, batch.rmEffect=F, batch.facLoc=x_DEG/100, # batch.facScale=x_DEG/100
                               )
colData(simulated_OPC)$cell_type = 'OPC'
colData(simulated_OPC)$Group = 'Group1'
#saveHDF5SummarizedExperiment(simulated_OPC, scenario_dir, prefix="simulated_OPC_")

#simulated_neurons = logNormCounts(simulated_neurons)
#simulated_neurons = runUMAP(simulated_neurons)
#g = plotUMAP(simulated_neurons, colour_by = "Group", shape_by = "Batch")
#ggsave(file.path(result_dir, 'simulated_neurons_t0_three_regions_UMAP.png'))
#
#simulated_Astrocytes = logNormCounts(simulated_Astrocytes)
#simulated_Astrocytes = runUMAP(simulated_Astrocytes)
#g = plotUMAP(simulated_Astrocytes, colour_by = "Group", shape_by = "Batch")
#ggsave(file.path(result_dir, 'simulated_Astrocytes_t0_three_regions_UMAP.png'))
#
#simulated_GABA = logNormCounts(simulated_GABA)
#simulated_GABA = runUMAP(simulated_GABA)
#simulated_OPC = logNormCounts(simulated_OPC)
#simulated_OPC = runUMAP(simulated_OPC)

rowData(simulated_neurons) = rowData(simulated_neurons)[, 'Gene', drop=F]
colnames(simulated_neurons) = paste0(colnames(simulated_neurons), '_', colData(simulated_neurons)$cell_type)
rowData(simulated_GABA) = rowData(simulated_GABA)[, 'Gene', drop=F]
colnames(simulated_GABA) = paste0(colnames(simulated_GABA), '_', colData(simulated_GABA)$cell_type)
rowData(simulated_Astrocytes) = rowData(simulated_Astrocytes)[, 'Gene', drop=F]
colnames(simulated_Astrocytes) = paste0(colnames(simulated_Astrocytes), '_', colData(simulated_Astrocytes)$cell_type)
rowData(simulated_OPC) = rowData(simulated_OPC)[, 'Gene', drop=F]
colnames(simulated_OPC) = paste0(colnames(simulated_OPC), '_', colData(simulated_OPC)$cell_type)

idx = colData(simulated_neurons)$Group == 'Group2' & colData(simulated_neurons)$Batch %in% c('Batch2', 'Batch3')
sub_simulated_neurons = simulated_neurons[, !idx]
idx = colData(simulated_GABA)$Group == 'Group2' & colData(simulated_GABA)$Batch %in% c('Batch2', 'Batch3')
sub_simulated_GABA = simulated_GABA[, !idx]
simulated_sce = cbind(sub_simulated_neurons, sub_simulated_GABA, simulated_Astrocytes, simulated_OPC)

simulated_sce = logNormCounts(simulated_sce)
simulated_sce = runUMAP(simulated_sce)
colData(simulated_sce)$label = paste0(colData(simulated_sce)$Group, '-', colData(simulated_sce)$Batch)
g = plotUMAP(simulated_sce, colour_by = "cell_type", shape_by = "Group")
ggsave(file.path(scenario_dir, 'simulated_sce_UMAP_Group.png'))
g = plotUMAP(simulated_sce, colour_by = "cell_type", shape_by = "Batch")
ggsave(file.path(scenario_dir, 'simulated_sce_UMAP_Batch.png'))
# saveHDF5SummarizedExperiment(simulated_sce, scenario_dir, prefix="simulated_sce")

## remove cells with no genes expressed
#idx_zerosize = apply(counts(simulated_sce), MARGIN = 2, sum) == 0
#simulated_sce = simulated_sce[, !idx_zerosize]
#colData(simulated_sce)$label = paste0(colData(simulated_sce)$Group, '-', colData(simulated_sce)$Batch)
#saveRDS(simulated_sce, file=file.path(scenario_dir, paste0('LIBD_res_seed_', seed, '.RDS')))
#simulated_sce = readRDS(file.path(scenario_dir, paste0('LIBD_res_seed_', seed, '.RDS')))

colData(simulated_sce)$label = paste0(colData(simulated_sce)$Group, '-', colData(simulated_sce)$Batch)
## assign count data
scenario_props = c(0.45, 0.45, 0.05, 0.05)
celltype_labels = c('Neurons', 'GABA', 'Astrocytes', 'OPC')
names(scenario_props) = celltype_labels
for (timepoint in names(tp_clusters)) {
    metadata = as.data.frame(colData(simulated_sce))
    if (timepoint == 't0') {
        subsample = location_res[[timepoint]][['subsample']]
        subsample$celltype = NA
        cnt = array(NA, c(nrow(simulated_sce), nrow(subsample)))
        # 1 -> forebraini (G1B1), 2 -> midbrain (G2B1 Neurons, G1B1 non-neurons), 3 -> hindbrain (G3B1 Neurons, G1B1 non-neurons)
        for (label in unique(subsample$label)) {
            label_idx = subsample$label == label
            num_cells = sum(label_idx)
            if (num_cells == 0) next
            final_cells_selected = c()
            sampled_number = 0 # sampled cell numbers
            for (celltype_idx in 1:length(unique(metadata$cell_type))) {
                celltype = unique(metadata$cell_type)[celltype_idx]
                if (celltype_idx != length(unique(metadata$cell_type))) {
                    prop_sampled_number = round(num_cells * unname(scenario_props[celltype]))
                    sampled_number = sampled_number + prop_sampled_number
                } else {
                    prop_sampled_number = num_cells - sampled_number
                }
                if (label == 1) {
                    cells_idx =  which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch1')
                }
                if (label == 2) {
                    if (celltype %in% c('Neurons', 'GABA')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group2-Batch1')
                    }
                    if (celltype %in% c('Astrocytes', 'OPC')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch1')
                    }
                }
                if (label == 3) {
                    if (celltype %in% c('Neurons', 'GABA')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group3-Batch1')
                    }
                    if (celltype %in% c('Astrocytes', 'OPC')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch1')
                    }
                }
                cells_selected = sample(cells_idx, prop_sampled_number, replace=F)
                final_cells_selected = c(final_cells_selected, cells_selected)
            }
            cnt[, label_idx] = as.matrix(counts(simulated_sce))[, final_cells_selected]
            subsample[label_idx, 'celltype'] = colData(simulated_sce[, final_cells_selected])$cell_type
            }
        }

    if (timepoint == 't1') {
        subsample = location_res[[timepoint]][['subsample']]
        subsample$celltype = NA
        cnt = array(NA, c(nrow(simulated_sce), nrow(subsample)))
        # 1 -> midbrain (G2B1 Neurons, G1B1 Non-Neurons); 2 -> forebrain2 (G1B3), 3 -> forebrain1 (G1B2), 4 -> hindbrain (G3B1 Neurons, G1B1 non-neurons)
        for (label in unique(subsample$label)) {
            label_idx = subsample$label == label
            num_cells = sum(label_idx)
            if (num_cells == 0) next
            final_cells_selected = c()
            sampled_number = 0 # sampled cell numbers
            for (celltype_idx in 1:length(unique(metadata$cell_type))) {
                celltype = unique(metadata$cell_type)[celltype_idx]
                if (celltype_idx != length(unique(metadata$cell_type))) {
                    prop_sampled_number = round(num_cells * unname(scenario_props[celltype]))
                    sampled_number = sampled_number + prop_sampled_number
                } else {
                    prop_sampled_number = num_cells - sampled_number
                }
                cat(label, celltype, prop_sampled_number, '\n')
                if (label == 1) {
                    if (celltype %in% c('Neurons', 'GABA')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group2-Batch1')
                    }
                    if (celltype %in% c('Astrocytes', 'OPC')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch1')
                    }
                }
                if (label == 2) {
                    cells_idx =  which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch3')
                }
                if (label == 3) {
                    cells_idx =  which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch2')
                }
                if (label == 4) {
                    if (celltype %in% c('Neurons', 'GABA')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group3-Batch1')
                    }
                    if (celltype %in% c('Astrocytes', 'OPC')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch1')
                    }
                }
                cells_selected = sample(cells_idx, prop_sampled_number, replace=F)
                final_cells_selected = c(final_cells_selected, cells_selected)
            }
            cnt[, label_idx] = as.matrix(counts(simulated_sce))[, final_cells_selected]
            subsample[label_idx, 'celltype'] = colData(simulated_sce[, final_cells_selected])$cell_type
            }
        }

    if (timepoint == 't2') {
        subsample = location_res[[timepoint]][['subsample']]
        subsample$celltype = NA
        cnt = array(NA, c(nrow(simulated_sce), nrow(subsample)))
        # 1 -> midbrain (G2B1 Neurons, G1B1 Non-Neurons); 2 -> hindbrain1 (G3B2 Neurons, G1B4 non-neurons), 
        # 3 -> forebrain2 (G1B3), 4 -> forebrain1 (G1B2), 5 -> hindbrain2 (G3B3 Neurons, G1B5 non-neurons)
        for (label in unique(subsample$label)) {
            label_idx = subsample$label == label
            num_cells = sum(label_idx)
            if (num_cells == 0) next
            final_cells_selected = c()
            sampled_number = 0 # sampled cell numbers
            for (celltype_idx in 1:length(unique(metadata$cell_type))) {
                celltype = unique(metadata$cell_type)[celltype_idx]
                if (celltype_idx != length(unique(metadata$cell_type))) {
                    prop_sampled_number = round(num_cells * unname(scenario_props[celltype]))
                    sampled_number = sampled_number + prop_sampled_number
                } else {
                    prop_sampled_number = num_cells - sampled_number
                }
                if (label == 1) {
                    if (celltype %in% c('Neurons', 'GABA')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group2-Batch1')
                    }
                    if (celltype %in% c('Astrocytes', 'OPC')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch1')
                    }
                }
                if (label == 2) {
                    if (celltype %in% c('Neurons', 'GABA')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group3-Batch2')
                    }
                    if (celltype %in% c('Astrocytes', 'OPC')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch4')
                    }
                }
                if (label == 3) {
                    cells_idx =  which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch3')
                }
                if (label == 4) {
                    cells_idx =  which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch2')
                }
                if (label == 5) {
                    if (celltype %in% c('Neurons', 'GABA')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group3-Batch3')
                    }
                    if (celltype %in% c('Astrocytes', 'OPC')) {
                        cells_idx = which(colData(simulated_sce)$cell_type == celltype & colData(simulated_sce)$label == 'Group1-Batch5')
                    }
                }
                cells_selected = sample(cells_idx, prop_sampled_number, replace=F)
                final_cells_selected = c(final_cells_selected, cells_selected)
            }
            cnt[, label_idx] = as.matrix(counts(simulated_sce))[, final_cells_selected]
            subsample[label_idx, 'celltype'] = colData(simulated_sce[, final_cells_selected])$cell_type
            }
        }
 
    colnames(cnt) = paste0('Cell', 1:nrow(subsample))
    rownames(cnt) = paste0('Gene', 1:nrow(simulated_sce))
    rownames(subsample) = colnames(cnt)
    subsample$label = as.factor(subsample$label)
    ggplot(subsample, aes(x=x_coor, y=y_coor, color=label)) + geom_point()
    ggsave(file.path(scenario_dir, paste0(timepoint, '_region_labels.png')))
    #for (region in unique(subsample$curated_label)) {
    #    region_subsample = subsample[subsample$curated_label == region, ]
    #    ggplot(region_subsample, aes(x=x_coor, y=y_coor, color=celltype)) + geom_point(alpha=0.5, size=1)
    #    ggsave(file.path(scenario_dir, paste0(timepoint, '_', region, '_celltype_labels_per_region.png')))
    #}
    ggplot(subsample, aes(x=x_coor, y=y_coor, color=celltype)) + geom_point()
    ggsave(file.path(scenario_dir, paste0(timepoint, '_celltype_labels.png')))
    writeMM(Matrix(cnt, sparse=T), file.path(scenario_dir, paste0(timepoint, '.mtx')))
    write.csv(subsample, file.path(scenario_dir, paste0(timepoint, '_metadata.csv')))
    # validate in Seurat
    obj = CreateSeuratObject(cnt, meta.data=subsample, assay='RNA')
    # but in SeuratDisk Convert(), which would turn categorical into numeric
    obj@meta.data$label = as.character(obj@meta.data$label)
    obj@meta.data$celltype = as.character(obj@meta.data$celltype)
    #obj[["RNA3"]] <- as(object = obj[["RNA"]], Class = "Assay")
    #DefaultAssay(obj) <- "RNA3"
    #obj[["RNA"]] <- NULL
    #obj <- RenameAssays(object = obj, RNA3 = 'RNA')
    #saveRDS(obj, file.path(scenario_dir, paste0(timepoint, "_brain_seed_", seed, ".RDS")))
    #SaveH5Seurat(obj, filename=file.path(scenario_dir, paste0(timepoint, "_brain.h5seurat")), overwrite=TRUE)
    #Convert(file.path(scenario_dir, paste0(timepoint, "_brain.h5seurat")), dest = "h5ad", overwrite=TRUE)
    ## test generated data
    obj = NormalizeData(obj)
    obj = FindVariableFeatures(obj)
    obj = ScaleData(obj)
    obj = RunPCA(obj)
    # for visualization purpose
    obj@meta.data$celltype = as.factor(obj@meta.data$celltype)
    DimPlot(obj, reduction='pca', group.by='celltype')
    ggsave(file.path(scenario_dir, paste0(timepoint, '_generated_SeuratObj_PCA.png')))
    obj = RunUMAP(obj, dims = 1:10)
    DimPlot(obj, reduction='umap', group.by='celltype')
    ggsave(file.path(scenario_dir, paste0(timepoint, '_generated_SeuratObj_UMAP.png')))
}



