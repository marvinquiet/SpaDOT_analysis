args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[1]
suppressPackageStartupMessages(library(lisi))
options(digits = 17)
stopifnot(formals(lisi::compute_lisi)$perplexity == 30, formals(lisi::compute_lisi)$nn_eps == 0)
d <- read.csv(file.path(out, "labels.csv"), stringsAsFactors = FALSE)
tp_counts <- c(t0=2397, t1=2312, t2=3546)
stopifnot(nrow(d) == sum(tp_counts), !anyNA(d))
reps <- c("GPVAE10", "GAT10", "combined20")
r <- do.call(rbind, lapply(names(tp_counts), function(tp) {
    dat <- d[d$timepoint == tp, ]
    stopifnot(nrow(dat) == tp_counts[[tp]])
    # Same rigid orientation used before compute_lisi in the Mouse downstream R script.
    loc <- data.frame(x = if (tp == "t1") dat$pixel_y else -dat$pixel_y, y = -dat$pixel_x)
    sapply(reps, function(rep) {
        md <- data.frame(label = as.character(dat[[rep]]))
        vals <- lisi::compute_lisi(loc, md, "label")$label
        stopifnot(length(vals) == nrow(dat), all(is.finite(vals)), all(vals >= 1 - 1e-10), all(vals <= 5 + 1e-10))
        median(vals)
    })
}))
base <- read.csv(file.path(out, "pas_chaos_ari.csv"), stringsAsFactors = FALSE)
base$median_LISI <- r[cbind(match(base$timepoint, names(tp_counts)), match(base$representation, reps))]
stopifnot(nrow(base) == 9L, all(is.finite(as.matrix(base[,c("median_LISI", "PAS", "CHAOS", "ARI")]))))
metrics <- c("median_LISI", "PAS", "CHAOS", "ARI")
avg <- aggregate(base[, metrics], list(representation = base$representation), mean)
write.csv(base[,c("representation", "timepoint", "stage", "n_spots", metrics)], file.path(out, "metrics_per_timepoint.csv"), row.names = FALSE)
write.csv(avg, file.path(out, "mean_metrics.csv"), row.names = FALSE)
capture.output(list(lisi = packageVersion("lisi"), RANN = packageVersion("RANN"), defaults = formals(lisi::compute_lisi)),
               file = file.path(out, "r_metric_versions.txt"))
