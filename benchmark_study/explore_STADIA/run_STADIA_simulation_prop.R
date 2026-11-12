## conda activate ~/envs/spatialATAC
.libPaths("/net/mulan/home/wenjinma/Rlib")

library(stadia)
library(dplyr)
library(Seurat)
library(aricode)
library(ggplot2)
library(optparse)
library(Matrix)

## on Mulan
work_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/benchmark_study/explore_STADIA"
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
# sample 9000 cells to accelerate results
sample_idx = sample(c(colnames(t0_obj), colnames(t1_obj), colnames(t2_obj)), 3000*3)
seuList = list()
seuVec = list(t0_obj, t1_obj, t2_obj)
seuMetaList = list()
for (i in 1:length(seuVec)) {
    common_cells = intersect(colnames(seuVec[[i]]), sample_idx)
    sub_obj = seuVec[[i]][, common_cells]
    sub_assay_v3 = CreateAssayObject(counts=sub_obj[['RNA']]$'counts')
    sub_obj[['RNA']] = sub_assay_v3
    seuList[[i]] = sub_obj
    seuMetaList[[i]] = sub_obj@meta.data
}
seuMeta = do.call(rbind, seuMetaList)
## set parameters for STADIA model
etas = 0.15
K = 5
hyper = HyperParameters(seuList, eta = etas)
start_time = proc.time()
stadia = stadia(seuList, hyper, n_cluster = K, min.features = 0, min.spots = 20, # do not remove cells based on feature expression
                platform = "others", em.maxiter = 30)
end_time = proc.time()
elapsed_time = end_time - start_time
write(unname(elapsed_time["elapsed"]), file.path(result_dir, 'STADIA_running_time.txt'))
save(stadia, file=file.path(result_dir, 'stadia.RData'))

# perform KMeans clustering on the embeddings
sample_res_all = NULL
ARI_res = c()
for (tp in unique(seuMeta$timepoint)) {
    tp_meta = seuMeta[seuMeta$timepoint == tp, ]
    common_cells = intersect(rownames(tp_meta), colnames(stadia$factors))
    tp_embed = t(stadia$factors)[common_cells, ]
    tp_meta = tp_meta[common_cells, ]
    tp_kmeans = kmeans(tp_embed, centers=length(unique(tp_meta$label)))
    tp_meta$cluster = tp_kmeans$cluster
    tp_ari = ARI(tp_meta$label, tp_meta$cluster)
    cat(tp, 'ARI:', tp_ari, '\n')
    ARI_res = c(ARI_res, tp_ari)
    # visualizing results
    g = ggplot(tp_meta, aes(row, col)) + 
        geom_point(size = .3, shape = 21, color = "#00000000", aes(fill = factor(cluster))) +
        labs(subtitle = paste0("ARI:", round(tp_ari, 2)))
    ggsave(file.path(result_dir, paste0(tp, '_results.png')))
    # save embeddings
    sample_res = cbind(tp_embed, tp_meta)
    if (is.null(sample_res_all)) {
        sample_res_all = sample_res
    } else {
        sample_res_all = rbind(sample_res_all, sample_res)
    }
}
write.csv(sample_res_all, file.path(result_dir, 'STADIA_results.csv'))
names(ARI_res) = unique(seuMeta$timepoint)
write.csv(data.frame(ARI_res), file.path(result_dir, 'STADIA_ARI.csv'))

