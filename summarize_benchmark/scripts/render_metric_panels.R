# Update the Figure 3B dot plot with matched per-stage domain scores.
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L, args[4] %in% c("MouseEmbryo", "ChickenHeart", "Axolotl"))
scores <- read.csv(args[1], stringsAsFactors = FALSE)
archived <- read.csv(args[2], stringsAsFactors = FALSE)
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Keep the original palette assignment (CSV order) and display order (factor levels).
original_colors <- c("#c1121f", "#d3b300", "#165451", "#f09232", "#4c2623", "#0082c0",
                     "#1d1815", "#68a756", "#ded6c1", "#fb6107", "#ccd5ae", "#ABDDDE")
original_methods <- c("SpaDOT", "spatiAlign", "spatialPCA", "spaVAE", "SPIRAL", "SEDR",
                      "STADIA", "PRECAST", "spaGCN", "STAligner", "GraphST", "STAMP")
archive_order <- unique(archived$method)
stopifnot(setequal(archive_order, original_methods), length(archive_order) == length(original_colors))
colors <- c(setNames(original_colors, archive_order), STIFT = "#6A4C93")
methods <- c(original_methods, "STIFT")
stages <- sort(setdiff(unique(scores$timepoint), "avg"))
expected_stages <- if (args[4] == "ChickenHeart") 4L else 3L
stopifnot(identical(stages, paste0("t", seq_len(expected_stages) - 1L)))
stage_shapes <- setNames(c(16, 17, 15, 18)[seq_along(stages)], stages)
metrics <- c("ARI", "AMI", "H", "median-LISI", "PAS", "CHAOS")
scores <- scores[scores$index %in% metrics, c("index", "timepoint", "value", "method")]
scores$value <- as.numeric(scores$value)
stopifnot(setequal(scores$method, methods), setequal(scores$timepoint, c(stages, "avg")),
          setequal(scores$index, metrics), all(is.finite(scores$value)),
          !anyDuplicated(scores[c("index", "timepoint", "method")]),
          nrow(scores) == length(methods) * (length(stages) + 1L) * length(metrics))
scores$method <- factor(scores$method, levels = methods)
scores$timepoint <- factor(scores$timepoint, levels = c(stages, "avg"))

summaries <- scores[scores$timepoint == "avg", ]
means <- aggregate(value ~ index + method, scores[scores$timepoint != "avg", ], mean)
means <- merge(means, summaries[c("index", "method", "value")], by = c("index", "method"))
stopifnot(nrow(means) == length(metrics) * length(methods),
          max(abs(means$value.x - means$value.y)) < 2e-6)

draw_metric <- function(metric) {
    stage_values <- scores[scores$index == metric & scores$timepoint != "avg", ]
    mean_values <- scores[scores$index == metric & scores$timepoint == "avg", ]
    reference <- mean_values$value[mean_values$method == "SpaDOT"]
    ggplot() +
        geom_hline(yintercept = reference, linetype = "dashed", color = colors[["SpaDOT"]],
                   linewidth = 0.65) +
        geom_point(data = stage_values, aes(method, value, shape = timepoint),
                   size = 1.65, alpha = 0.30) +
        geom_point(data = mean_values, aes(method, value, color = method),
                   size = 2.75, alpha = 0.80) +
        scale_color_manual(values = colors, breaks = methods) +
        scale_shape_manual(values = stage_shapes, breaks = stages) +
        scale_x_discrete(drop = FALSE, expand = expansion(add = 0.55)) +
        scale_y_continuous(expand = expansion(mult = c(0.07, 0.08))) +
        labs(x = NULL, y = metric) +
        theme_pubr(base_size = 8) +
        theme(axis.text.x = element_blank(),
              axis.title.y = element_text(size = 8, margin = margin(r = 2)),
              plot.margin = margin(7, 2.5, 0, 1.5, unit = "pt"),
              legend.position = "bottom",
              legend.box = "horizontal",
              legend.margin = margin(0),
              legend.spacing.x = grid::unit(0.8, "mm"),
              legend.box.spacing = grid::unit(0, "mm"),
              legend.key.width = grid::unit(3.1, "mm"),
              legend.key.height = grid::unit(3.25, "mm"),
              legend.text = element_text(size = 6.7, margin = margin(r = 1, unit = "pt"))) +
        guides(color = guide_legend(nrow = 3, title = NULL, order = 1,
                                    override.aes = list(size = 2.8, alpha = 0.8)),
               shape = guide_legend(nrow = 3, title = NULL, order = 2,
                                    override.aes = list(size = 2.1, alpha = 0.35)))
}

# Keep SpaDOT centered vertically beside a four-column comparator legend.
legend_entries <- function(labels, x, y, fills, shapes = 16, size = 2.5) {
    grid::grobTree(
        grid::pointsGrob(grid::unit(x, "npc"), grid::unit(y, "npc"),
                         pch = shapes, size = grid::unit(size, "mm"),
                         gp = grid::gpar(col = fills)),
        grid::textGrob(labels, grid::unit(x, "npc") + grid::unit(5.3, "pt"),
                       grid::unit(y, "npc"), just = "left",
                       gp = grid::gpar(fontfamily = "sans", fontsize = 6.7))
    )
}

comparators <- setdiff(methods, "SpaDOT")
stopifnot(length(comparators) == 12L, "STIFT" %in% comparators)
legend_rows <- c(0.78, 0.50, 0.22)
stage_rows <- if (expected_stages == 4L) c(0.86, 0.62, 0.38, 0.14) else legend_rows
legend <- grid::grobTree(
    legend_entries("SpaDOT", 0.118, 0.50, grDevices::adjustcolor(colors[["SpaDOT"]], alpha.f = 0.8), size = 3.5),
    grid::segmentsGrob(x0 = 0.245, x1 = 0.245, y0 = 0.20, y1 = 0.80,
                       gp = grid::gpar(col = "#C9C9C9", lwd = 0.65)),
    legend_entries(comparators, rep(c(0.280, 0.424, 0.568, 0.712), each = 3L),
                   rep(legend_rows, 4L), grDevices::adjustcolor(colors[comparators], alpha.f = 0.8)),
    legend_entries(stages, rep(0.866, expected_stages), stage_rows, rep("#A5A5A5", expected_stages),
                   shapes = unname(stage_shapes), size = 1.85)
)

draw_triplet <- function(chosen, filename) {
    panels <- lapply(chosen, function(metric) draw_metric(metric) + theme(legend.position = "none"))
    strip <- cowplot::plot_grid(plotlist = panels, nrow = 1, align = "h", axis = "tb")
    combined <- cowplot::plot_grid(strip, legend, ncol = 1, rel_heights = c(1, 0.29))
    ggsave(file.path(out, filename), combined, width = 14, height = 6.45, units = "cm",
           device = "pdf", bg = "white")
}

if (args[4] != "MouseEmbryo") {
    draw_triplet(c("ARI", "AMI", "H"), paste0(args[4], "_ARI_AMI_H.pdf"))
}
draw_triplet(c("median-LISI", "PAS", "CHAOS"), paste0(args[4], "_spatial_continuity.pdf"))
