# library(devtools)
# install_github('theislab/kBET')
# install_github("immunogenomics/lisi")

library(kBET)
library(lisi)
library(optparse)

# set optparse
option_list = list(
    make_option(c("--result_dir"), type='character', default=NULL, help="Result directory."),
    make_option(c("--method"), type="character", default=NULL, 
        help="Spatial clustering method")
)
parser = OptionParser(option_list=option_list)
opt = parse_args(parser)
result_dir = opt$result_dir
method = opt$method

latent_data = read.csv(file.path(result_dir, paste0(method, '_latent.csv')), row.names=1, header=T)
batch_data = read.csv(file.path(result_dir, paste0(method, '_batch.csv')), row.names=1, header=T)
# random sample at 10%, 20%, ..., 50% of the data to compute rejection rate
percents = c(0.1, 0.2, 0.3, 0.4, 0.5)
set.seed(1234)
kBET_res = c()
for (pt in percents) {
    cat('Downsample for kBET:', pt, '\n')
    sample_idx = sample(1:nrow(latent_data), round(nrow(latent_data)*pt))
    latent_data_sample = latent_data[sample_idx, ]
    batch_data_sample = batch_data[sample_idx, 1]
    kBET_res = c(kBET_res, 1-kBET(as.matrix(latent_data_sample), batch_data_sample, plot=FALSE)$average.pval) # acceptance rate, higher the better
}
kBET_res_df = data.frame(acceptance_rate=kBET_res, sample_precent=percents)
write.csv(kBET_res_df, file.path(result_dir, paste0(method, '_kBET.csv')), row.names=T)

# compute lisi
colnames(batch_data) = "timepoint"
lisi_res = compute_lisi(as.matrix(latent_data), batch_data, c('timepoint'))
write(median(lisi_res$timepoint), file.path(result_dir, paste0(method, '_lisi.txt')))

