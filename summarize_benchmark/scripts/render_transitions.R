# Run the existing dataset-specific plot_OT_heatmap() helpers for Figures 3c–5c.
# Rscript render_transitions.R INPUT OUTPUT [BENCHMARK_SCRIPTS]
# INPUT is the paper OT command's output/inputs directory, including display keys.
# BENCHMARK_SCRIPTS defaults to the existing dataset scripts next to this file.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
    stop("Usage: Rscript render_transitions.R INPUT OUTPUT [BENCHMARK_SCRIPTS]")
}
inputs <- normalizePath(args[1])
output <- args[2]
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
scripts <- if (length(args) == 3L) args[3] else dirname(script_file)
scripts <- normalizePath(scripts)
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(reshape2)
})
invisible(Sys.setlocale("LC_COLLATE", "C"))

find_assignment <- function(expressions, variable) {
    found <- Filter(function(expr) {
        is.call(expr) && length(expr) == 3L &&
            (identical(expr[[1]], as.name("=")) || identical(expr[[1]], as.name("<-"))) &&
            identical(expr[[2]], as.name(variable))
    }, as.list(expressions))
    stopifnot(length(found) == 1L)
    found[[1]]
}

read_matrix <- function(path) {
    as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
}

datasets <- data.frame(
    dataset = c("MouseEmbryo", "ChickenHeart", "AxolotlTelencephalon"),
    directory = c("MouseEmbryo_Visium_Organogenesis", "ChickenHeart_Visium_Cardiogenesis",
                  "AxolotlTelencephalon_Stereoseq_Developmental"),
    helper = c("ATACVisium_utils.R", "ChickenHeart_utils.R", "Axolotl_utils.R"),
    downstream = c("analyze_ATACvisium_noremoval_downstream.R", "analyze_ChickenHeart_downstream.R",
                   "analyze_Axolotl_downstream.R"),
    figure = c("F3", "F4", "F5"),
    stages = c(3L, 4L, 3L)
)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
plot_inputs <- file.path(output, "plot_inputs")
dir.create(plot_inputs, showWarnings = FALSE)
audit <- list()
for (i in seq_len(nrow(datasets))) {
    dataset <- datasets[i, ]
    helper_dir <- file.path(scripts, dataset$directory)
    env <- new.env()
    eval(find_assignment(parse(file.path(helper_dir, dataset$helper)), "plot_OT_heatmap"), env)
    for (variable in c("timepoints_dict", "cluster_mapping")) {
        eval(find_assignment(as.list(body(env$plot_OT_heatmap)), variable), env)
    }

    input_dir <- file.path(inputs, dataset$dataset)
    plot_dir <- file.path(plot_inputs, dataset$dataset)
    dir.create(plot_dir, showWarnings = FALSE)
    key <- read.csv(file.path(input_dir, "domain_key.csv"), colClasses = "character")
    palette <- setNames(key$color, key$paper_label)
    stopifnot(all(vapply(split(key$color, key$paper_label), function(x) length(unique(x)) == 1L, logical(1))))
    palette <- palette[!duplicated(names(palette))]
    if (dataset$dataset != "ChickenHeart") {
        eval(find_assignment(parse(file.path(helper_dir, dataset$downstream)), "tp_color_palette"), env)
        stopifnot(identical(tolower(unname(palette)), tolower(unname(env$tp_color_palette[names(palette)]))))
        palette <- env$tp_color_palette
    }
    times <- paste0("t", seq_len(dataset$stages) - 1L)
    names_for_stage <- unname(env$timepoints_dict[times])
    permutation <- setNames(vector("list", length(times)), times)
    paper_labels <- setNames(vector("list", length(times)), times)
    for (j in seq_along(times)) {
        stage <- names_for_stage[j]
        original <- if (times[j] %in% names(env$cluster_mapping)) {
            env$cluster_mapping[[times[j]]]
        } else {
            env$cluster_mapping[[stage]]
        }
        original <- original[as.character(seq_along(original) - 1L)]
        current <- key[key$stage == stage, ]
        stopifnot(setequal(unname(original), current$paper_label))
        permutation[[j]] <- current$domain_id[match(original, current$paper_label)]
        paper_labels[[j]] <- unname(original)
        stopifnot(!anyNA(palette[paper_labels[[j]]]))
    }

    expected <- list()
    for (j in seq_len(length(times) - 1L)) {
        suffix <- paste(j - 1L, j, sep = "_")
        for (direction in c("row", "col")) {
            input_name <- if (direction == "row") "forward_conditional" else "reverse_conditional"
            mat <- read_matrix(file.path(input_dir, paste0(input_name, "_", suffix, ".csv")))
            mat <- mat[permutation[[j]], permutation[[j + 1L]], drop = FALSE]
            write.csv(mat, file.path(plot_dir, paste0("transition_prob_", suffix, "_", direction, "_norm.csv")))
        }
        score <- read_matrix(file.path(input_dir, paste0("association_score_", suffix, ".csv")))
        score <- score[permutation[[j]], permutation[[j + 1L]], drop = FALSE]
        dimnames(score) <- paper_labels[c(j, j + 1L)]
        expected[[paste(times[j], times[j + 1L], sep = "_")]] <- score
    }

    # Capture the original pair plots and legends without modifying their appearance.
    env$plots <- list()
    env$ggsave <- eval(quote(function(filename, plot = ggplot2::last_plot(), ...) {
        plots[[basename(filename)]] <<- plot
        invisible(filename)
    }), env)
    env$plot_OT_heatmap(times, "spaDOT", plot_dir, tp_color_palette = palette, prob_cutoff = 0.20)
    stopifnot(length(env$plots) == 2L * length(expected))
    panels <- legends <- list()
    for (j in seq_along(expected)) {
        pair <- names(expected)[j]
        plot_file <- names(env$plots)[startsWith(names(env$plots), paste0(pair, "_OT_heatmap"))]
        panel <- env$plots[[plot_file[!grepl("legend", plot_file)]]]
        legend <- env$plots[[plot_file[grepl("legend", plot_file)]]]
        frame <- panel$data
        score <- expected[[j]]
        values <- score[cbind(match(frame[[times[j]]], rownames(score)),
                              match(frame[[times[j + 1L]]], colnames(score)))]
        build <- ggplot_build(panel)
        axes <- build$layout$panel_params[[1]]
        stopifnot(nrow(frame) == length(score), !anyDuplicated(frame[1:2]),
                  length(values) == nrow(frame), max(abs(values - frame$value)) < 2e-15,
                  identical(is.na(frame$value_edited), values < 0.20),
                  identical(as.character(axes$y$get_labels()), sort(rownames(score))),
                  identical(as.character(axes$x$get_labels()), sort(colnames(score))),
                  identical(unname(panel$theme$axis.text.x$colour), unname(palette[sort(colnames(score))])),
                  identical(unname(panel$theme$axis.text.y$colour), unname(palette[sort(rownames(score))])))
        audit[[length(audit) + 1L]] <- data.frame(
            dataset = dataset$dataset, pair = pair,
            source = names_for_stage[j], target = names_for_stage[j + 1L],
            rows_bottom_to_top = paste(sort(rownames(score)), collapse = " "),
            columns_left_to_right = paste(sort(colnames(score)), collapse = " "),
            cells = length(score), highlighted = sum(values >= 0.20),
            max_abs_error = max(abs(values - frame$value))
        )
        panels[[j]] <- panel
        legends[[j]] <- legend
    }
    panels <- cowplot::align_plots(plotlist = panels, align = "h", axis = "tb")
    pairs <- Map(function(panel, legend) {
        cowplot::plot_grid(panel, legend, nrow = 1, rel_widths = c(4, 2))
    }, panels, legends)
    combined <- cowplot::plot_grid(plotlist = pairs, nrow = 1)
    ggsave(file.path(output, paste0(dataset$figure, "_", dataset$dataset, "_transition.pdf")),
           combined, width = 6 * length(pairs), height = 4, bg = "white")
}
write.csv(do.call(rbind, audit), file.path(output, "plot_audit.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "R_session.txt"))
