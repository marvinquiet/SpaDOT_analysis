#' simulate proportion changes during timepoint split, with region-specific NPCs in each region
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
work_dir = '/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_simulation'
setwd(work_dir)
source('utils.R')

# set optparse
option_list = list(
    make_option(c("-s", "--seed"), type="integer", default=2024, help="random seed"),
    make_option(c("-p", "--prop"), type="double", default=0.5, help="proportion of major cell type in each region"),
    make_option(c("--DEG"), type="integer", default=500, help="number of DE genes for each cell type")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)

data_dir = '/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain'
param_dir = '/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop'
result_dir = paste0('/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_regionNPCs_DEG', opt$DEG)
dir.create(result_dir, showWarnings = FALSE)

seed = opt$seed
prop = opt$prop
seed_dir = file.path(result_dir, seed, paste0('prop_', prop))
dir.create(seed_dir, showWarnings=F, recursive = TRUE)

#' estimate for once and use it later to do simulation
# NPCs_init_params = splatter_estimate_params('/net/zootopia/disk1/wenjinma/data/spatial_simulation', cell_type='Stem cells')
# saveRDS(NPCs_init_params, file=file.path(param_dir, "NPCs_init_params_pFC_for_LIBD.RDS"))
NPCs_init_params = readRDS(file.path(param_dir, "NPCs_init_params_pFC_for_LIBD.RDS"))
#Neurons_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='Neurons')
#saveRDS(Neurons_init_params, file = file.path(param_dir, "Neurons_init_params_pFC_for_LIBD.RDS"))
Neurons_init_params = readRDS(file.path(param_dir, "Neurons_init_params_pFC_for_LIBD.RDS"))
#GABANeurons_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='GABAergic neurons')
#saveRDS(GABANeurons_init_params, file = file.path(param_dir, "GABA_init_params_pFC_for_LIBD.RDS"))
GABANeurons_init_params = readRDS(file.path(param_dir, "GABA_init_params_pFC_for_LIBD.RDS"))
#Astrocytes_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='Astrocytes')
#saveRDS(Astrocytes_init_params, file = file.path(param_dir, "Astrocytes_init_params_pFC_for_LIBD.RDS"))
Astrocytes_init_params = readRDS(file.path(param_dir, "Astrocytes_init_params_pFC_for_LIBD.RDS"))
#OPC_init_params = splatter_estimate_params(file.path(base_dir, "data/spatial_simulation"), cell_type='OPC')
#saveRDS(OPC_init_params, file = file.path(param_dir, "OPC_init_params_pFC_for_LIBD.RDS"))
OPC_init_params = readRDS(file.path(param_dir, "OPC_init_params_pFC_for_LIBD.RDS"))
# Micro_init_params = splatter_estimate_params('/net/zootopia/disk1/wenjinma/data/spatial_simulation', cell_type='Microglia')
# saveRDS(Micro_init_params, file=file.path(param_dir, "Micro_init_params_pFC_for_LIBD.RDS"))
Micro_init_params = readRDS(file.path(param_dir, "Micro_init_params_pFC_for_LIBD.RDS"))

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
    #ggsave(file.path(param_dir, paste0(timepoint, '_simulated_image.png')))
    # save results to list
    #saveRDS(res, file.path(param_dir, paste0(timepoint, '_simulated_image_data.RDS')))
    res = readRDS(file.path(param_dir, paste0(timepoint, '_simulated_image_data.RDS')))
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

## generate simulated data for t0 with one cell type - stem cells in each brain region
# simulated_NPCs = splatter_simu_brain(NPCs_init_params, nCells=20000, nGenes=5000, seed=seed, out.prob=0, nGroups=1, nBatches=1)
# colData(simulated_NPCs)$cell_type = 'NPC'
# simulated_OPCs = splatter_simu_brain(OPC_init_params, nCells=5000, nGenes=5000, seed=seed, out.prob=0, nGroups=1, nBatches=1)
# colData(simulated_OPCs)$cell_type = 'OPC'
# simulated_Astrocytes = splatter_simu_brain(Astrocytes_init_params, nCells=5000, nGenes=5000, seed=seed, out.prob=0, nGroups=1, nBatches=1)
# colData(simulated_Astrocytes)$cell_type = 'Astrocytes'
# simulated_neurons = splatter_simu_brain(Neurons_init_params, nCells=5000, nGenes=5000, seed=seed, out.prob=0, nGroups=1, nBatches=1)
# colData(simulated_neurons)$cell_type = 'Neurons'
# simulated_GABA = splatter_simu_brain(GABANeurons_init_params, nCells=5000, nGenes=5000, seed=seed, out.prob=0, nGroups=1, nBatches=1)
# colData(simulated_GABA)$cell_type = 'GABA'
# simualted_Micro = splatter_simu_brain(Micro_init_params, nCells=15000, nGenes=5000, seed=seed, out.prob=0, nGroups=1, nBatches=1)
# colData(simualted_Micro)$cell_type = 'Microglia'

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

# rowData(simulated_NPCs) = rowData(simulated_NPCs)[, 'Gene', drop=F]
# colnames(simulated_NPCs) = paste0(colnames(simulated_NPCs), '_', colData(simulated_NPCs)$cell_type)
# rowData(simulated_neurons) = rowData(simulated_neurons)[, 'Gene', drop=F]
# colnames(simulated_neurons) = paste0(colnames(simulated_neurons), '_', colData(simulated_neurons)$cell_type)
# rowData(simulated_GABA) = rowData(simulated_GABA)[, 'Gene', drop=F]
# colnames(simulated_GABA) = paste0(colnames(simulated_GABA), '_', colData(simulated_GABA)$cell_type)
# rowData(simulated_Astrocytes) = rowData(simulated_Astrocytes)[, 'Gene', drop=F]
# colnames(simulated_Astrocytes) = paste0(colnames(simulated_Astrocytes), '_', colData(simulated_Astrocytes)$cell_type)
# rowData(simulated_OPCs) = rowData(simulated_OPCs)[, 'Gene', drop=F]
# colnames(simulated_OPCs) = paste0(colnames(simulated_OPCs), '_', colData(simulated_OPCs)$cell_type)
# rowData(simualted_Micro) = rowData(simualted_Micro)[, 'Gene', drop=F]
# colnames(simualted_Micro) = paste0(colnames(simualted_Micro), '_', colData(simualted_Micro)$cell_type)


# --- total 5,000 genes among 6 cell types, but 100 genes each for different cell types
simulated_sce = splatter_simu_brain(NPCs_init_params, nCells=65000, nGenes=5000, seed=seed, out.prob=0, nGroups=1, nBatches=1)
colData(simulated_sce)$cell_type = 'None'
simulated_NPCs = splatter_simu_brain(NPCs_init_params, nCells=10000, nGenes=opt$DEG, seed=seed, out.prob=0, nGroups=3, nBatches=1)
colData(simulated_NPCs)$cell_type = 'NPC'
simulated_OPCs = splatter_simu_brain(OPC_init_params, nCells=5000, nGenes=opt$DEG, seed=seed, out.prob=0, nGroups=1, nBatches=1)
colData(simulated_OPCs)$cell_type = 'OPC'
simulated_Astrocytes = splatter_simu_brain(Astrocytes_init_params, nCells=5000, nGenes=opt$DEG, seed=seed, out.prob=0, nGroups=1, nBatches=1)
colData(simulated_Astrocytes)$cell_type = 'Astrocytes'
simulated_neurons = splatter_simu_brain(Neurons_init_params, nCells=5000, nGenes=opt$DEG, seed=seed, out.prob=0, nGroups=1, nBatches=1)
colData(simulated_neurons)$cell_type = 'Neurons'
simulated_GABA = splatter_simu_brain(GABANeurons_init_params, nCells=5000, nGenes=opt$DEG, seed=seed, out.prob=0, nGroups=1, nBatches=1)
colData(simulated_GABA)$cell_type = 'GABA'
simualted_Micro = splatter_simu_brain(Micro_init_params, nCells=15000, nGenes=opt$DEG, seed=seed, out.prob=0, nGroups=1, nBatches=1)
colData(simualted_Micro)$cell_type = 'Microglia'


# 1:200 NPC specific genes, 20000 cells
counts(simulated_sce)[1:opt$DEG, 1:30000] = as.matrix(counts(simulated_NPCs))
simulated_sce$cell_type[1:30000] = paste0(simulated_NPCs$cell_type, '-', simulated_NPCs$Group)
# 501:1000 Neurons specific genes
counts(simulated_sce)[(opt$DEG+1):(2*opt$DEG), 30001:35000] = as.matrix(counts(simulated_neurons))
simulated_sce$cell_type[30001:35000] = simulated_neurons$cell_type
# 1001:1500 GABA specific genes
counts(simulated_sce)[(2*opt$DEG+1):(3*opt$DEG), 35001:40000] = as.matrix(counts(simulated_GABA))
simulated_sce$cell_type[35001:40000] = simulated_GABA$cell_type
# 1501:2000 Astrocytes specific genes
counts(simulated_sce)[(3*opt$DEG+1):(4*opt$DEG), 40001:45000] = as.matrix(counts(simulated_Astrocytes))
simulated_sce$cell_type[40001:45000] = simulated_Astrocytes$cell_type
# 2001:2500 OPC specific genes
counts(simulated_sce)[(4*opt$DEG+1):(5*opt$DEG), 45001:50000] = as.matrix(counts(simulated_OPCs))
simulated_sce$cell_type[45001:50000] = simulated_OPCs$cell_type
# 2501:3000 Microglia specific genes
counts(simulated_sce)[(5*opt$DEG+1):(6*opt$DEG), 50001:65000] = as.matrix(counts(simualted_Micro))
simulated_sce$cell_type[50001:65000] = simualted_Micro$cell_type

# simulated_sce = cbind(simulated_NPCs, simulated_neurons, simulated_GABA, simulated_Astrocytes, simulated_OPCs, simualted_Micro)
simulated_sce = logNormCounts(simulated_sce)
simulated_sce = runPCA(simulated_sce)
# u = uwot::umap(t(logcounts(simulated_sce)), n_neighbors = 2)
# reducedDim(simulated_sce, "UMAP_uwot") = u
g = plotPCA(simulated_sce, colour_by = "cell_type")
ggsave(file.path(seed_dir, 'simulated_sce_PCA_Group.png'))

## remove cells with no genes expressed
idx_zerosize = apply(counts(simulated_sce), MARGIN = 2, sum) == 0
simulated_sce = simulated_sce[, !idx_zerosize]
#saveRDS(simulated_sce, file=file.path(seed_dir, paste0('LIBD_res_seed_', seed, '.RDS')))
#simulated_sce = readRDS(file.path(seed_dir, paste0('LIBD_res_seed_', seed, '.RDS')))

# --- renamed cell_type
## assign count data
celltype_labels = c('Neurons', 'GABA', 'Astrocytes', 'OPC', 'NPC-Group1', 'NPC-Group2', 'NPC-Group3', 'Microglia')
for (timepoint in names(tp_clusters)) {
    metadata = as.data.frame(colData(simulated_sce))
    if (timepoint == 't0') {
        subsample = location_res[[timepoint]][['subsample']]
        subsample$celltype = NA
        cnt = array(NA, c(nrow(simulated_sce), nrow(subsample)))
        # 1 -> forebraini (NPC, OPC, Astrocytes), 2 -> midbrain (NPC), 3 -> hindbrain (NPC, Neurons, GABA)
        for (label in unique(subsample$label)) {
            label_idx = subsample$label == label
            num_cells = sum(label_idx)
            if (num_cells == 0) next
            final_cells_selected = c()
            sampled_number = 0 # sampled cell numbers
            for (celltype_idx in 1:length(unique(metadata$cell_type))) {
                celltype = unique(metadata$cell_type)[celltype_idx]
                if (label == 1) {
                    scenario_props = c('NPC-Group1'=0.5, 'NPC-Group2'=0, 'NPC-Group3'=0, 'OPC'=0.25, 'Astrocytes'=0.25, 'Neurons'=0, 'GABA'=0, 'Microglia'=0)
                }
                if (label == 2) {
                    scenario_props = c('NPC-Group1'=0, 'NPC-Group2'=1-prop, 'NPC-Group3'=0, 'Neurons'=0, 'GABA'=0, 'Astrocytes'=0, 'OPC'=0, 'Microglia'=prop)
                }
                if (label == 3) {
                    scenario_props = c('NPC-Group1'=0, 'NPC-Group2'=0, 'NPC-Group3'=0.5, 'Neurons'=0.25, 'GABA'=0.25, 'Astrocytes'=0, 'OPC'=0, 'Microglia'=0)
                }
                cells_idx =  which(colData(simulated_sce)$cell_type == celltype)
                if (celltype_idx != length(unique(metadata$cell_type))) {
                    prop_sampled_number = round(num_cells * unname(scenario_props[celltype]))
                    sampled_number = sampled_number + prop_sampled_number
                    if (sampled_number > num_cells) {
                        prop_sampled_number = prop_sampled_number - (sampled_number - num_cells)
                        sampled_number = num_cells
                    }
                } else {
                    prop_sampled_number = num_cells - sampled_number
                }
                if (prop_sampled_number > 0) {
                    cells_selected = sample(cells_idx, prop_sampled_number, replace=F)
                    final_cells_selected = c(final_cells_selected, cells_selected)
                }
            }
            cnt[, label_idx] = as.matrix(counts(simulated_sce))[, final_cells_selected]
            subsample[label_idx, 'celltype'] = colData(simulated_sce[, final_cells_selected])$cell_type
        }
    }
    if (timepoint == 't1') {
        subsample = location_res[[timepoint]][['subsample']]
        subsample$celltype = NA
        cnt = array(NA, c(nrow(simulated_sce), nrow(subsample)))
        # 1 -> midbrain (NPCs); 2 -> forebrain2 (NPC, Astrocytes), 3 -> forebrain1 (NPC, OPC), 4 -> hindbrain (NPC, Neurons, GABA)
        for (label in unique(subsample$label)) {
            label_idx = subsample$label == label
            num_cells = sum(label_idx)
            if (num_cells == 0) next
            final_cells_selected = c()
            sampled_number = 0 # sampled cell numbers
            for (celltype_idx in 1:length(unique(metadata$cell_type))) {
                celltype = unique(metadata$cell_type)[celltype_idx]
                if (label == 1) {
                    scenario_props = c('NPC-Group1'=0, 'NPC-Group2'=1-prop, 'NPC-Group3'=0, 'Neurons'=0, 'GABA'=0, 'Astrocytes'=0, 'OPC'=0, 'Microglia'=prop)
                }
                if (label == 2) {
                    scenario_props = c('NPC-Group1'=1-prop, 'NPC-Group2'=0, 'NPC-Group3'=0, 'OPC'=0, 'Astrocytes'=prop, 'Neurons'=0, 'GABA'=0, 'Microglia'=0)
                }
                if (label == 3) {
                    scenario_props = c('Neurons'=0, 'GABA'=0, 'NPC-Group1'=1-prop, 'NPC-Group2'=0, 'NPC-Group3'=0, 'Astrocytes'=0, 'OPC'=prop, 'Microglia'=0)
                }
                if (label == 4) {
                    scenario_props = c('NPC-Group1'=0, 'NPC-Group2'=0, 'NPC-Group3'=0.5, 'Neurons'=0.25, 'GABA'=0.25, 'Astrocytes'=0, 'OPC'=0, 'Microglia'=0)
                }
                cells_idx =  which(colData(simulated_sce)$cell_type == celltype)
                if (celltype_idx != length(unique(metadata$cell_type))) {
                    prop_sampled_number = round(num_cells * unname(scenario_props[celltype]))
                    sampled_number = sampled_number + prop_sampled_number
                    if (sampled_number > num_cells) {
                        prop_sampled_number = prop_sampled_number - (sampled_number - num_cells)
                        sampled_number = num_cells
                    }
                } else {
                    prop_sampled_number = num_cells - sampled_number
                }
                if (prop_sampled_number > 0) {
                    cells_selected = sample(cells_idx, prop_sampled_number, replace=F)
                    final_cells_selected = c(final_cells_selected, cells_selected)
                }
            }
            cnt[, label_idx] = as.matrix(counts(simulated_sce))[, final_cells_selected]
            subsample[label_idx, 'celltype'] = colData(simulated_sce[, final_cells_selected])$cell_type
        }
    }
    if (timepoint == 't2') {
        subsample = location_res[[timepoint]][['subsample']]
        subsample$celltype = NA
        cnt = array(NA, c(nrow(simulated_sce), nrow(subsample)))
        # 1 -> midbrain (NPCs); 2 -> hindbrain1 (NPC, Neurons), 
        # 3 -> forebrain2 (NPC, Astrocytes), 4 -> forebrain1 (NPC, OPC), 5 -> hindbrain2 (NPC, GABA)
        for (label in unique(subsample$label)) {
            label_idx = subsample$label == label
            num_cells = sum(label_idx)
            if (num_cells == 0) next
            final_cells_selected = c()
            sampled_number = 0 # sampled cell numbers
            for (celltype_idx in 1:length(unique(metadata$cell_type))) {
                celltype = unique(metadata$cell_type)[celltype_idx]
                if (label == 1) {
                    scenario_props = c('NPC-Group1'=0, 'NPC-Group2'=1-prop, 'NPC-Group3'=0, 'Neurons'=0, 'GABA'=0, 'Astrocytes'=0, 'OPC'=0, 'Microglia'=prop)
                }
                if (label == 2) {
                    scenario_props = c('NPC-Group1'=0, 'NPC-Group2'=0, 'NPC-Group3'=1-prop, 'OPC'=0, 'Astrocytes'=0, 'Neurons'=prop, 'GABA'=0, 'Microglia'=0)
                }
                if (label == 3) {
                    scenario_props = c('Neurons'=0, 'GABA'=0, 'NPC-Group1'=1-prop, 'NPC-Group2'=0, 'NPC-Group3'=0, 'Astrocytes'=prop, 'OPC'=0, 'Microglia'=0)
                }
                if (label == 4) {
                    scenario_props = c('NPC-Group1'=1-prop, 'NPC-Group2'=0, 'NPC-Group3'=0, 'Neurons'=0, 'GABA'=0, 'Astrocytes'=0, 'OPC'=prop, 'Microglia'=0)
                }
                if (label == 5) {
                    scenario_props = c('NPC-Group1'=0, 'NPC-Group2'=0, 'NPC-Group3'=1-prop, 'Neurons'=0, 'GABA'=prop, 'Astrocytes'=0, 'OPC'=0, 'Microglia'=0)
                }
                cells_idx =  which(colData(simulated_sce)$cell_type == celltype)
                if (celltype_idx != length(unique(metadata$cell_type))) {
                    prop_sampled_number = round(num_cells * unname(scenario_props[celltype]))
                    sampled_number = sampled_number + prop_sampled_number
                    if (sampled_number > num_cells) {
                        prop_sampled_number = prop_sampled_number - (sampled_number - num_cells)
                        sampled_number = num_cells
                    }
                } else {
                    prop_sampled_number = num_cells - sampled_number
                }
                if (prop_sampled_number > 0) {
                    cells_selected = sample(cells_idx, prop_sampled_number, replace=F)
                    final_cells_selected = c(final_cells_selected, cells_selected)
                }
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
    ggsave(file.path(seed_dir, paste0(timepoint, '_region_labels.png')))
    for (region in unique(subsample$label)) {
       region_subsample = subsample[subsample$label == region, ]
       ggplot(region_subsample, aes(x=x_coor, y=y_coor, color=celltype)) + geom_point(alpha=0.5, size=1)
       ggsave(file.path(seed_dir, paste0(timepoint, '_', region, '_celltype_labels_per_region.png')))
    }
    ggplot(subsample, aes(x=x_coor, y=y_coor, color=celltype)) + geom_point()
    ggsave(file.path(seed_dir, paste0(timepoint, '_celltype_labels.png')))
    writeMM(Matrix(cnt, sparse=T), file.path(seed_dir, paste0(timepoint, '.mtx')))
    write.csv(subsample, file.path(seed_dir, paste0(timepoint, '_metadata.csv')))
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
    ggsave(file.path(seed_dir, paste0(timepoint, '_generated_SeuratObj_PCA.png')))
    obj = RunUMAP(obj, dims = 1:10)
    DimPlot(obj, reduction='umap', group.by='celltype')
    ggsave(file.path(seed_dir, paste0(timepoint, '_generated_SeuratObj_UMAP.png')))
}
