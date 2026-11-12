# Updated Chicken Heart spatial domains on the manuscript's paired coordinate frame.
suppressPackageStartupMessages(library(ggplot2))
library(grid)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
input <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
key <- read.csv(args[2], colClasses = "character")
pred_palette <- setNames(key$color[!duplicated(key$paper_label)], key$paper_label[!duplicated(key$paper_label)])
pred_palette <- pred_palette[c("1","1-1","2","3","3-1","3-2","3-3","4","5","5-1","5-2","6")]
truth_palette <- c("Epicardium- like"="#000000", "Ventricle"="#00b4d8", "Valves"="#ffb3c6",
                   "Outflow tract"="#fca311", "Atria"="#a3b18a", "Epicardium"="#000000",
                   "Endothelium"="#344e41", "Trabecular LV and endocardium"="#caf0f8",
                   "Right ventricle"="#0077b6", "Compact LV and inter-ventricular septum"="#90e0ef",
                   "Endocardium"="#48cae4")
truth_text <- setNames(names(truth_palette), names(truth_palette))
truth_text["Epicardium- like"] <- "Epicardium-like"
truth_text["Trabecular LV and endocardium"] <- "Trabecular LV and\nendocardium"
truth_text["Compact LV and inter-ventricular septum"] <- "Compact LV and\ninter-ventricular septum"
stages <- c(t0="Day 4", t1="Day 7", t2="Day 10", t3="Day 14")
counts <- c(t0=747L, t1=1966L, t2=1916L, t3=1967L)
stopifnot(nrow(input) == sum(counts), !anyDuplicated(input$spot),
          setequal(unique(input$SpaDOT), key$domain_id))
input$x <- input$pixel_y
input$y <- input$pixel_x
input$truth <- trimws(gsub("[[:space:]]+", " ", input$true_labels))
# Display the Day 14 annotation exactly as the original F4A plotting block does.
input$truth[input$timepoint == "t3" & input$truth == "Trabecular LV and endocardium"] <- "Endocardium"
input$pred <- key$paper_label[match(input$SpaDOT, key$domain_id)]
stopifnot(all(input$truth %in% names(truth_palette)), all(input$pred %in% names(pred_palette)))

panel <- function(data, field, colors, x_limits, y_limits) {
    ggplot(data, aes(x, y, colour = .data[[field]])) +
        geom_point(size = 0.50) +
        scale_color_manual(values = colors) +
        coord_fixed(ratio = 1, xlim = x_limits, ylim = y_limits, expand = FALSE) +
        theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
}
stage_data <- lapply(names(stages), function(t) {
    data <- input[input$timepoint == t, ]
    stopifnot(nrow(data) == counts[[t]], length(unique(data$pred)) == sum(key$stage == stages[[t]]))
    half <- max(diff(range(data$x)), diff(range(data$y))) * 0.522
    x_limits <- mean(range(data$x)) + c(-half, half)
    y_limits <- mean(range(data$y)) + c(-half, half)
    truth <- panel(data, "truth", truth_palette, x_limits, y_limits)
    pred <- panel(data, "pred", pred_palette, x_limits, y_limits)
    tr <- ggplot_build(truth)
    pr <- ggplot_build(pred)
    stopifnot(identical(tr$data[[1]][c("x", "y")], pr$data[[1]][c("x", "y")]),
              nrow(tr$data[[1]]) == counts[[t]])
    list(truth = ggplotGrob(truth), pred = ggplotGrob(pred), truth_labels = unique(data$truth),
         pred_labels = unique(data$pred),
         audit = data.frame(timepoint=t, stage=stages[[t]], spots=nrow(data), clusters=length(unique(data$pred)),
                            xmin=x_limits[1], xmax=x_limits[2], ymin=y_limits[1], ymax=y_limits[2],
                            coordinate_ratio=1, paired_coordinates_identical=TRUE))
})
put <- function(row, col, content) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    force(content)
    popViewport()
}
draw_stage_legend <- function(labels, palette, display, cols, text_size, row_height) {
    ordered <- names(palette)[names(palette) %in% labels]
    rows <- ceiling(length(ordered) / cols)
    row_id <- (seq_along(ordered)-1L) %% rows
    col_id <- (seq_along(ordered)-1L) %/% rows
    total_width <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE)
    widths <- vapply(seq_len(cols)-1L, function(j) {
        items <- ordered[col_id == j]
        if (!length(items)) return(0)
        max(vapply(items, function(x) {
            convertWidth(grobWidth(textGrob(display[[x]], gp=gpar(fontsize=text_size))),
                         "in", valueOnly=TRUE)
        }, numeric(1))) + 0.13
    }, numeric(1))
    spacing <- 0.08
    block_width <- sum(widths) + spacing*(cols-1L)
    stopifnot(block_width <= total_width)
    left <- (total_width-block_width)/2 + cumsum(c(0, head(widths, -1L) + spacing))
    y <- unit(0.5, "npc") + unit((rows-1)/2*row_height - row_id*row_height, "in")
    for (i in seq_along(ordered)) {
        label <- ordered[i]
        x <- left[col_id[i]+1L]
        grid.circle(x=unit(x+0.032,"in"), y=y[i], r=unit(0.032,"in"),
                    gp=gpar(col=NA, fill=palette[[label]]))
        grid.text(display[[label]], x=unit(x+0.095,"in"), y=y[i], just="left",
                  gp=gpar(fontsize=text_size, lineheight=0.90))
    }
}

heights <- c(0.30, 0.24, 1.75, 0.83, 0.10, 0.30, 0.24, 1.75, 0.40)
grDevices::pdf(args[3], width = 10.0, height = sum(heights)+0.18, useDingbats = FALSE,
               title="Chicken Heart spatial domains | SpaDOT 1.0.7")
grid.newpage()
pushViewport(viewport(width=unit(9.70,"in"), height=unit(sum(heights),"in"),
                      layout=grid.layout(length(heights), 4L, heights=unit(heights,"in"))))
put(1L, 1:4, grid.text("Region Annotation", gp=gpar(fontsize=15)))
put(6L, 1:4, grid.text("SpaDOT Prediction", gp=gpar(fontsize=15)))
pred_text <- setNames(names(pred_palette), names(pred_palette))
for (i in seq_along(stages)) {
    d <- stage_data[[i]]
    for (r in c(2L,7L)) put(r, i, grid.text(stages[[i]], gp=gpar(fontsize=11)))
    put(3L,i,grid.draw(d$truth))
    put(4L,i,draw_stage_legend(d$truth_labels, truth_palette, truth_text, 2L, 8.0, 0.18))
    put(8L,i,grid.draw(d$pred))
    put(9L,i,draw_stage_legend(d$pred_labels, pred_palette, pred_text, 4L, 9.0, 0.17))
}
popViewport()
invisible(dev.off())
write.csv(do.call(rbind,lapply(stage_data,function(x)x$audit)), args[4], row.names=FALSE)
