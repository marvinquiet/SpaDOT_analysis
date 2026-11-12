## conda activate /media/sda/wenjinma/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(PRECAST)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)
library(optparse)
library(Matrix)

## on Mulan
work_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_PRECAST"
setwd(work_dir)

# set optparse
option_list = list(
    make_option(c("-s", "--seed"), type="integer", default=2024, help="random seed"),
    make_option(c("--prop"), type="double", default=0.8, help="Major cell type proportion"),
    make_option(c("--DEG"), type="integer", default=500, help="Number of DEGs per cell type"),
    make_option(c("--mixedNPCs"), action="store_true", default=FALSE, help="Whether the NPCs is a mixture")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)

seed = opt$seed
prop = opt$prop
DEG = opt$DEG
if (opt$mixedNPCs) {
    scenario_dir = file.path(paste0("/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_mixedNPCs_DEG", DEG), seed, paste('prop', prop, sep='_'))
    result_dir = file.path(work_dir, paste0('simulation_prop_results_mixedNPCs_DEG', DEG), seed, paste('prop', prop, sep='_'))
} else {
    scenario_dir = file.path(paste0("/net/zootopia/disk1/wenjinma/data/spatial_simulation_brain_prop_DEG", DEG), seed, paste('prop', prop, sep='_'))
    result_dir = file.path(work_dir, paste0('simulation_prop_results_DEG', DEG), seed, paste('prop', prop, sep='_'))
}
dir.create(result_dir, recursive=T)

# load data
t0_mtx = readMM(file.path(scenario_dir, 't0.mtx'))
t0_metadata = read.csv(file.path(scenario_dir, 't0_metadata.csv'), header=T, row.names=1)
t0_metadata['timepoint'] = 't0'
t0_metadata = t0_metadata %>% rename(row=x_coor, col=y_coor)
t0_obj = CreateSeuratObject(t0_mtx, meta.data=t0_metadata)
colnames(t0_obj) = paste0('t0_', colnames(t0_obj))
t1_mtx = readMM(file.path(scenario_dir, 't1.mtx'))
t1_metadata = read.csv(file.path(scenario_dir, 't1_metadata.csv'), header=T, row.names=1)
t1_metadata['timepoint'] = 't1'
t1_metadata = t1_metadata %>% rename(row=x_coor, col=y_coor)
t1_obj = CreateSeuratObject(t1_mtx, meta.data=t1_metadata)
colnames(t1_obj) = paste0('t1_', colnames(t1_obj))
t2_mtx = readMM(file.path(scenario_dir, 't2.mtx'))
t2_metadata = read.csv(file.path(scenario_dir, 't2_metadata.csv'), header=T, row.names=1)
t2_metadata['timepoint'] = 't2'
t2_metadata = t2_metadata %>% rename(row=x_coor, col=y_coor)
t2_obj = CreateSeuratObject(t2_mtx, meta.data=t2_metadata)
colnames(t2_obj) = paste0('t2_', colnames(t2_obj))
# sample 6000 cells to accelerate results
set.seed(seed)
sample_idx = sample(c(colnames(t0_obj), colnames(t1_obj), colnames(t2_obj)), 6000, replace=F)
seuList = list()
seuVec = list(t0_obj, t1_obj, t2_obj)
seuMetaList = list()
for (i in 1:length(seuVec)) {
    common_cells = intersect(colnames(seuVec[[i]]), sample_idx)
    sub_obj = seuVec[[i]][, common_cells]
    seuList[[i]] = sub_obj
    seuMetaList[[i]] = sub_obj@meta.data
}
seuMeta = do.call(rbind, seuMetaList)
# run PRECAST
start_time = proc.time()
preobj = CreatePRECASTObject(seuList = seuList, selectGenesMethod = "HVGs", gene.number = 2000) 
PRECASTObj = AddAdjList(preobj, platform = "Visium")
PRECASTObj = AddParSetting(PRECASTObj, Sigma_equal = TRUE, coreNum = 1, maxIter = 30, verbose = TRUE)
PRECASTObj = PRECAST(PRECASTObj, K = 5) # set 5 for now because in total there will be 5 clusters
end_time = proc.time()
elapsed_time = end_time - start_time
write(unname(elapsed_time["elapsed"]), file.path(result_dir, 'PRECAST_running_time.txt'))
## backup the fitting results in resList
resList = PRECASTObj@resList
PRECASTObj = SelectModel(PRECASTObj)
seuInt = IntegrateSpaData(PRECASTObj, species='Unknown')
seuInt = AddTSNE(seuInt, n_comp = 2) 
seuInt@meta.data = seuMeta[rownames(seuInt@meta.data), ]
# perform KMeans clustering on the embeddings
sample_res_all = NULL
ARI_res = c()
for (tp in unique(seuInt$timepoint)) {
    tp_seu = subset(seuInt, subset=timepoint == tp)
    tp_seu_embed = Embeddings(tp_seu, reduction='PRECAST')
    tp_kmeans = kmeans(tp_seu_embed, centers=length(unique(tp_seu$label)))
    tp_seu$cluster = tp_kmeans$cluster
    tp_ari = ARI(tp_seu$label, tp_seu$cluster)
    cat(tp, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # plot embeddings
    cols_cluster = chooseColors(palettes_name = 'Nature 10', n_colors = length(unique(tp_seu$label)), plot_colors = F)
    p1 = dimPlot(tp_seu, item='celltype', point_size = 1,  font_family='serif')
    p2 = dimPlot(tp_seu, item='label', point_size = 1,  font_family='serif')
    p3 = dimPlot(tp_seu, item='cluster', font_family='serif', cols=cols_cluster) # Times New Roman
    drawFigs(list(p1, p2, p3), common.legend=FALSE, align='hv', layout.dim=c(3,1)) 
    ggsave(file.path(result_dir, paste0(tp, '_tSNE.png')))
    # save sample res
    sample_res = cbind(tp_seu_embed, tp_seu@meta.data)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'PRECAST_results.csv'))
names(ARI_res) = unique(seuInt$timepoint)
write.csv(data.frame(ARI_res), file.path(result_dir, 'PRECAST_ARI.csv'))

