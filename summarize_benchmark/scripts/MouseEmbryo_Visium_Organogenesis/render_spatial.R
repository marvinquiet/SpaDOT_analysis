# Draw freshly clustered SpaDOT 1.0.7 domains with shared 1:1 geometry.
suppressPackageStartupMessages(library(ggplot2))
library(grid)
args <- commandArgs(trailingOnly = TRUE)
input <- read.csv(args[1], colClasses = c("character", "numeric", "numeric", "character", "character"))
source_file <- parse(args[2])
extract_assignment <- function(symbol) {
    entries <- Filter(function(expr) {
        is.call(expr) && length(expr) == 3L && identical(expr[[1]], as.name("=")) &&
            identical(expr[[2]], as.name(symbol))
    }, as.list(source_file))
    stopifnot(length(entries) == 1L)
    eval(entries[[1]][[3]], envir = baseenv())
}

stages <- c(t0 = "E12.5", t1 = "E13.5", t2 = "E15.5")
pred_palette <- extract_assignment("tp_color_palette")
truth_palette <- c("CNS" = "#F8766D", "Heart" = "#A3A500", "Liver" = "#00BF7D",
                   "Men/PNS" = "#00B0F6", "Mesenchyme" = "#E76BF3")

# Preserve the stage-specific orientation in analyze_ATACvisium_noremoval_downstream.R.
input$x <- ifelse(input$timepoint == "t1", input$pixel_y, -input$pixel_y)
input$y <- -input$pixel_x

plot_stage <- function(points, field, palette, xlim, ylim) {
    ggplot(points, aes(x = x, y = y, color = .data[[field]])) +
        geom_point(size = 0.5) +
        scale_color_manual(values = palette) +
        coord_fixed(ratio = 1, xlim = xlim, ylim = ylim, expand = FALSE) +
        theme_void(base_family = "sans") +
        theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
}

panels <- vector("list", length(stages))
audit <- vector("list", length(stages))
for (i in seq_along(stages)) {
    stage <- names(stages)[i]
    points <- input[input$timepoint == stage, ]
    stopifnot(!anyNA(points$prediction), all(points$prediction %in% names(pred_palette)))
    span <- max(diff(range(points$x)), diff(range(points$y)))
    halfspan <- span * 0.525
    xlim <- mean(range(points$x)) + c(-halfspan, halfspan)
    ylim <- mean(range(points$y)) + c(-halfspan, halfspan)
    gt <- plot_stage(points, "true_labels", truth_palette, xlim, ylim)
    pred <- plot_stage(points, "prediction", pred_palette, xlim, ylim)
    gtb <- ggplot_build(gt)
    predb <- ggplot_build(pred)
    gp <- gtb$layout$panel_params[[1]]
    pp <- predb$layout$panel_params[[1]]
    stopifnot(identical(gtb$data[[1]][c("x", "y")], predb$data[[1]][c("x", "y")]),
              identical(gp$x.range, pp$x.range), identical(gp$y.range, pp$y.range),
              isTRUE(all.equal(diff(gp$x.range), diff(gp$y.range))), nrow(gtb$data[[1]]) == nrow(points))
    panels[[i]] <- list(gt = ggplotGrob(gt), pred = ggplotGrob(pred))
    audit[[i]] <- data.frame(stage = stages[[i]], spots = nrow(points),
                            raw_width = diff(range(points$x)), raw_height = diff(range(points$y)),
                            shared_x_min = xlim[1], shared_x_max = xlim[2],
                            shared_y_min = ylim[1], shared_y_max = ylim[2],
                            coordinate_ratio = 1, identical_spot_coordinates = TRUE)
}

draw_at <- function(row, col, expr) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    force(expr)
    popViewport()
}

draw_legend <- function(labels, palette) {
    text_gp <- gpar(fontfamily = "sans", fontsize = 12)
    widths <- vapply(labels, function(label) convertWidth(grobWidth(textGrob(label, gp = text_gp)), "in", valueOnly = TRUE), numeric(1))
    gap <- 0.18
    spacing <- 0.13
    total <- sum(widths + spacing) + (length(labels) - 1) * gap
    cursor <- (convertWidth(unit(1, "npc"), "in", valueOnly = TRUE) - total) / 2
    for (i in seq_along(labels)) {
        grid.circle(x = unit(cursor + 0.035, "in"), r = unit(0.038, "in"),
                    gp = gpar(col = NA, fill = unname(palette[[labels[i]]])))
        grid.text(labels[i], x = unit(cursor + spacing, "in"), just = "left", gp = text_gp)
        cursor <- cursor + widths[i] + spacing + gap
    }
}

heights <- c(0.30, 0.23, 2.18, 0.34, 0.08, 0.30, 0.23, 2.18, 0.34)
grDevices::pdf(args[3], width = 8.15, height = sum(heights) + 0.22, useDingbats = FALSE,
               title = "Figure 3A | Mouse Embryo spatial domains")
grid.newpage()
pushViewport(viewport(width = unit(7.85, "in"), height = unit(sum(heights), "in"),
                      layout = grid.layout(length(heights), 3L, heights = unit(heights, "in"))))
draw_at(1L, 1:3, grid.text("Region Annotation", gp = gpar(fontsize = 16)))
draw_at(6L, 1:3, grid.text("SpaDOT Prediction", gp = gpar(fontsize = 16)))
draw_at(1L, 1L, grid.text("a.", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 16)))
for (i in seq_along(stages)) {
    for (row in c(2L, 7L)) draw_at(row, i, grid.text(stages[[i]], gp = gpar(fontsize = 13, fontface = "bold")))
    draw_at(3L, i, grid.draw(panels[[i]]$gt))
    draw_at(8L, i, grid.draw(panels[[i]]$pred))
}
draw_at(4L, 1:3, draw_legend(names(truth_palette), truth_palette))
draw_at(9L, 1:3, draw_legend(names(pred_palette), pred_palette))
popViewport()
invisible(dev.off())
write.csv(do.call(rbind, audit), args[4], row.names = FALSE)
