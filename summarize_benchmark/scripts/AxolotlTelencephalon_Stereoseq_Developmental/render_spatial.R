# Figure 5A uses the original F5 palette and axes. Each stage shares one
# coordinate window and physical panel size across annotation and prediction.
suppressPackageStartupMessages(library(ggplot2))
library(grid)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
data <- read.csv(args[1], stringsAsFactors = FALSE, colClasses = c("character", "numeric", "numeric", "character", "character", "character"))
stages <- c(t0 = "St.44", t1 = "St.54", t2 = "St.57")
expected <- c(t0 = 1477, t1 = 2929, t2 = 4410)
gt_colors <- c("MP&DP&LP" = "#292444", "VZ" = "#ECCBAE", "CP" = "#EE7E81",
               "VLMC" = "#971E24", "Septum" = "#F3DF6C", "Striatum" = "#F98400",
               "MP&DP" = "#55597A", "LP" = "#A6C5D0", "MP" = "#567B87", "DP" = "#74AAB8")
pred_colors <- c("1" = "#EE7E81", "2" = "#F3DF6C", "3" = "#971E24",
                 "4" = "#292444", "4-1" = "#A6C5D0", "4-2" = "#55597A",
                 "4-2-1" = "#567B87", "4-2-2" = "#74AAB8", "5" = "#F98400", "6" = "#ECCBAE")
stopifnot(!anyNA(data), all(data$true_labels %in% names(gt_colors)),
          all(data$paper_label %in% names(pred_colors)), nrow(data) == sum(expected))
make_plot <- function(points, field, palette, xlim, ylim) {
    ggplot(points, aes(pixel_x, pixel_y, color = .data[[field]])) +
        geom_point(size = 0.29) + scale_color_manual(values = palette) +
        coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
        theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
}
panels <- lapply(names(stages), function(stage) {
    d <- data[data$timepoint == stage, ]
    stopifnot(nrow(d) == expected[[stage]])
    span <- max(diff(range(d$pixel_x)), diff(range(d$pixel_y)))
    delta <- c(-1, 1) * 0.528 * span
    xlim <- mean(range(d$pixel_x)) + delta
    ylim <- mean(range(d$pixel_y)) + delta
    gt <- make_plot(d, "true_labels", gt_colors, xlim, ylim)
    pred <- make_plot(d, "paper_label", pred_colors, xlim, ylim)
    g <- ggplot_build(gt)
    p <- ggplot_build(pred)
    g_coords <- g$data[[1]][c("x", "y")]
    p_coords <- p$data[[1]][c("x", "y")]
    gx <- g$layout$panel_params[[1]]$x.range
    gy <- g$layout$panel_params[[1]]$y.range
    stopifnot(identical(g_coords, p_coords), identical(gx, p$layout$panel_params[[1]]$x.range),
              identical(gy, p$layout$panel_params[[1]]$y.range), isTRUE(all.equal(diff(gx), diff(gy))))
    list(gt = ggplotGrob(gt), pred = ggplotGrob(pred), stage = stages[[stage]])
})
draw_at <- function(row, col, drawing) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    force(drawing)
    popViewport()
}
draw_legend <- function(top, bottom, palette) {
    gp <- gpar(fontsize = 7.5)
    n <- length(top)
    stopifnot(n == length(bottom))
    legend_width <- vapply(seq_len(n), function(i) {
        max(convertWidth(grobWidth(textGrob(top[i], gp = gp)), "cm", valueOnly = TRUE),
            convertWidth(grobWidth(textGrob(bottom[i], gp = gp)), "cm", valueOnly = TRUE)) + 0.42
    }, numeric(1))
    left <- (convertWidth(unit(1, "npc"), "cm", valueOnly = TRUE) - sum(legend_width)) / 2
    for (i in seq_len(n)) {
        for (k in 1:2) {
            label <- c(top[i], bottom[i])[k]
            y <- c(0.73, 0.27)[k]
            grid.circle(x = unit(left + 0.08, "cm"), y = unit(y, "npc"), r = unit(0.059, "cm"),
                        gp = gpar(col = NA, fill = unname(palette[[label]])))
            grid.text(label, x = unit(left + 0.21, "cm"), y = unit(y, "npc"), just = "left", gp = gp)
        }
        left <- left + legend_width[i]
    }
}
# Fit the two rows into the same slide footprint as the updated Figure 3A.
heights <- c(0.16, 0.48, 0.33, 3.78, 0.81, 0.54, 0.48, 0.33, 3.78, 0.81, 0.36)
stopifnot(isTRUE(all.equal(sum(heights[1:5]), 5.56)), isTRUE(all.equal(sum(heights[6:11]), 6.30)))
grDevices::pdf(args[2], width = 14 / 2.54, height = 11.86 / 2.54, useDingbats = FALSE,
               title = "Figure 5A | Axolotl telencephalon spatial domains")
grid.newpage()
pushViewport(viewport(width = unit(13.52, "cm"), height = unit(sum(heights), "cm"),
                      layout = grid.layout(length(heights), 3L, heights = unit(heights, "cm"))))
draw_at(2, 1:3, grid.text("Region Annotation", gp = gpar(fontsize = 11.5)))
draw_at(7, 1:3, grid.text("SpaDOT Prediction", gp = gpar(fontsize = 11.5)))
for (i in seq_along(panels)) {
    for (row in c(3, 8)) draw_at(row, i, grid.text(panels[[i]]$stage, gp = gpar(fontsize = 8.7, fontface = "bold")))
    draw_at(4, i, grid.draw(panels[[i]]$gt))
    draw_at(9, i, grid.draw(panels[[i]]$pred))
}
draw_at(5, 1:3, draw_legend(c("CP", "LP", "MP&DP", "Septum", "VLMC"),
                           c("DP", "MP", "MP&DP&LP", "Striatum", "VZ"), gt_colors))
draw_at(10, 1:3, draw_legend(c("1", "3", "4-1", "4-2-1", "5"),
                            c("2", "4", "4-2", "4-2-2", "6"), pred_colors))
popViewport()
invisible(dev.off())
cat("F5A coordinates, axis ranges, aspect ratios and all spots verified\n")
