# Recompute median spatial LISI with the original benchmark package and defaults.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
suppressPackageStartupMessages(library(lisi))
options(digits = 17)
stopifnot(as.numeric(formals(lisi::compute_lisi)$perplexity) == 30,
          as.numeric(formals(lisi::compute_lisi)$nn_eps) == 0)
labels <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
methods <- setdiff(names(labels), c("spot", "timepoint", "pixel_x", "pixel_y", "true_labels"))
stages <- c(t0 = 1477L, t1 = 2929L, t2 = 4410L)
stopifnot(!anyDuplicated(labels$spot),
          length(methods) == 13L, nrow(labels) == sum(stages), !anyNA(labels))
rows <- lapply(methods, function(method) {
    values <- vapply(names(stages), function(stage) {
        data <- labels[labels$timepoint == stage, ]
        stopifnot(nrow(data) == stages[[stage]])
        xy <- data[, c("pixel_x", "pixel_y")]
        pred <- data.frame(kmeans = as.character(data[[method]]))
        score <- lisi::compute_lisi(xy, pred, "kmeans")$kmeans
        recoded <- data.frame(domain = as.character(as.integer(factor(pred$kmeans))))
        check <- lisi::compute_lisi(xy, recoded, "domain")$domain
        stopifnot(length(score) == nrow(data), all(is.finite(score)),
                  all(score >= 1 - 1e-6), all(score <= length(unique(pred$kmeans)) + 1e-6),
                  max(abs(score - check)) < 1e-10)
        median(score)
    }, numeric(1))
    data.frame(index = "median-LISI", timepoint = c(names(stages), "avg"),
               value = c(values, mean(values)), method = method)
})
write.csv(do.call(rbind, rows), args[2], row.names = FALSE)
cat("LISI", as.character(packageVersion("lisi")), "verified for", length(methods) * length(stages), "method-stage pairs\n")
