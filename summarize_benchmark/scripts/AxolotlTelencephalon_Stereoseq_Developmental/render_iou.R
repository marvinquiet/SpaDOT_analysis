#!/usr/bin/env Rscript
# Figure 5e. Compute overlap and IoU for every predicted domain, then use
# the historical Hungarian assignment maximizing intersection / predicted size.
#
# Rscript render_iou.R PREDICTIONS DOMAIN_KEY OUTPUT_DIR
# PREDICTIONS is collect_predictions.py's predictions_lisi.csv, with all
# 13 methods on matched spots. DOMAIN_KEY is the stage-specific SpaDOT display
# key used by the same collector. The key changes only the displayed label.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop("Usage: Rscript render_iou.R PREDICTIONS DOMAIN_KEY OUTPUT_DIR")
}
data <- read.csv(args[1], check.names = FALSE, colClasses = "character")
key <- read.csv(args[2], check.names = FALSE, colClasses = "character")
out <- args[3]
if (dir.exists(out)) stop("Use a new output directory.")

palette <- c(SpaDOT = "#c1121f", spatiAlign = "#ded6c1",
             spatialPCA = "#d3b300", spaVAE = "#0082c0", SPIRAL = "#fb6107",
             SEDR = "#f09232", STADIA = "#ccd5ae", PRECAST = "#68a756",
             spaGCN = "#4c2623", STAligner = "#1d1815", GraphST = "#165451",
             STAMP = "#ABDDDE", STIFT = "#6F67AB")
methods <- names(palette)
stages <- c("t0", "t1", "t2")
targets <- c("t0 - MP&DP&LP", "t1 - MP&DP", "t1 - LP",
             "t2 - MP", "t2 - DP", "t2 - LP")
stopifnot(all(c("timepoint", "true_labels", methods) %in% names(data)),
          all(c("domain_id", "paper_label") %in% names(key)),
          !anyDuplicated(key$domain_id), setequal(data$SpaDOT, key$domain_id),
          setequal(data$timepoint, stages))
data$SpaDOT <- key$paper_label[match(data$SpaDOT, key$domain_id)]
stopifnot(!anyNA(data[c("timepoint", "true_labels", methods)]))

score_domains <- function(data, method, stage) {
    predicted <- data[[method]]
    truth <- data$true_labels
    # Preserve first-occurrence ordering for the same assignment tie rule as
    # the original analyze_Axolotl_downstream.R.
    intersections <- table(factor(predicted, levels = unique(predicted)),
                           factor(truth, levels = unique(truth)))
    overlap <- intersections / rowSums(intersections)
    union <- outer(rowSums(intersections), colSums(intersections), "+") - intersections
    iou <- intersections / union
    assignment <- as.integer(clue::solve_LSAP(1 - overlap))
    rows <- seq_len(nrow(intersections))
    full <- expand.grid(pred_labels = rownames(intersections),
                        annotation = colnames(intersections), stringsAsFactors = FALSE)
    full$method <- method
    full$timepoint <- stage
    full$overlap <- as.vector(overlap)
    full$IoU <- as.vector(iou)
    full <- full[c("method", "pred_labels", "overlap", "IoU", "timepoint", "annotation")]

    matched <- full[rows + nrow(intersections) * (assignment - 1L), , drop = FALSE]
    matched$assigned_annotation <- colnames(intersections)[assignment]
    matched$timepoint_annotation <- paste(stage, matched$annotation, sep = " - ")
    list(full = full, matched = matched)
}

results <- unlist(lapply(stages, function(stage) {
    stage_data <- data[data$timepoint == stage, , drop = FALSE]
    lapply(methods, function(method) score_domains(stage_data, method, stage))
}), recursive = FALSE)
all_pairs <- do.call(rbind, lapply(results, function(result) result$full))
matched <- do.call(rbind, lapply(results, function(result) result$matched))
shown <- matched[matched$timepoint_annotation %in% targets, ]
shown$timepoint_annotation <- factor(shown$timepoint_annotation, levels = targets)
shown$method <- factor(shown$method, levels = methods)
stopifnot(all(table(shown$timepoint_annotation, shown$method) == 1L))

dir.create(out, recursive = TRUE)
write.csv(all_pairs, file.path(out, "IoU_summary.csv"), row.names = FALSE)
write.csv(shown, file.path(out, "IoU_MPDPLP_matched.csv"), row.names = FALSE)

spadot <- shown[shown$method == "SpaDOT", ]
other <- shown[shown$method != "SpaDOT", ]
p <- ggplot() +
    geom_point(data = other, aes(x = overlap, y = IoU, color = method, shape = timepoint),
               size = 2) +
    geom_vline(data = spadot, aes(xintercept = overlap), colour = palette["SpaDOT"],
               linetype = "dashed", alpha = 0.5) +
    geom_hline(data = spadot, aes(yintercept = IoU), colour = palette["SpaDOT"],
               linetype = "dashed", alpha = 0.5) +
    geom_point(data = spadot, aes(x = overlap, y = IoU, color = method, shape = timepoint),
               size = 4) +
    facet_wrap(~timepoint_annotation, nrow = 2, scales = "free") +
    scale_color_manual(values = palette) +
    theme_pubr() + theme(legend.position = "none") +
    labs(x = "Overlap Percentage", y = "IoU")
ggsave(file.path(out, "IoU_vs_Percentage.pdf"), plot = p, width = 8, height = 4)
