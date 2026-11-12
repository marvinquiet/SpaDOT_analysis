# Fig. S42b. Usage: Rscript render_s42b.R <prepared-input-dir> <output-dir>
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('Pass prepared input and output directories')
input <- args[[1]]
output <- args[[2]]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
metrics_df <- read.csv(file.path(input, 'real_transition_scores.csv'))
metrics_df$method = factor(metrics_df$method, levels=c('SpaDOT','PRECAST','SEDR','SpaGCN','GraphST (positional 3000)','SpatialPCA','spatiAlign','spaVAE'))
levels(metrics_df$method)[levels(metrics_df$method) == "GraphST (positional 3000)"] <- "GraphST"
metrics_df$Data = factor(
  metrics_df$Data,
  levels = c("Mouse Embryo", 
             "Chicken Heart", 
             "Axolotl Telencephalon"),
  labels = c("Mouse\nOrganogenesis",
             "Chicken\nCardiogenesis",
             "Axolotl Telencephalon\nNeurogenesis")
)


gap_df = metrics_df %>%
  group_by(Data, Day) %>%
  summarise(
    SpaDOT_value = naive_accuracy[method == "SpaDOT"],
    best_other = max(naive_accuracy[method != "SpaDOT"], na.rm = TRUE),
    log2FC = log2(SpaDOT_value / best_other),
    .groups = "drop"
  )

plot_df = metrics_df %>%
  left_join(gap_df, by = c("Data", "Day"))

transition_labels <- data.frame(
  Data = c(
    rep("Mouse\nOrganogenesis", 2),
    rep("Chicken\nCardiogenesis", 3),
    rep("Axolotl Telencephalon\nNeurogenesis", 2)
  ),
  Day = c("0->1", "1->2", "0->1", "1->2", "2->3", "0->1", "1->2"),
  label = c(
    "E12.5->E13.5", "E13.5->E15.5",
    "D4->D7", "D7->D10", "D10->D14",
    "St.44->St.54", "St.54->St.57"
  )
)

plot_df <- plot_df %>%
  left_join(transition_labels, by = c("Data", "Day"))

plot_df$label = factor(plot_df$label, levels=c(
    "E12.5->E13.5", "E13.5->E15.5",
    "D4->D7", "D7->D10", "D10->D14",
    "St.44->St.54", "St.54->St.57"
  ))

plot_df$Data = factor(plot_df$Data, levels=c("Mouse\nOrganogenesis", "Chicken\nCardiogenesis", "Axolotl Telencephalon\nNeurogenesis"))


method_palette <- read.csv(file.path(input, 'method_palette.csv'), stringsAsFactors=FALSE)
method_colors <- setNames(method_palette$hex,method_palette$method)
method_colors['SpaGCN'] <- method_colors['spaGCN']
method_colors['SpatialPCA'] <- method_colors['spatialPCA']
stopifnot(all(levels(plot_df$method) %in% names(method_colors)))

p <- ggplot(plot_df, aes(x = label, y = naive_accuracy, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = .9), na.rm = TRUE) +
  geom_point(aes(y = ifelse(status == "failed", .025, NA_real_),
                 group = method, shape = "Failed"),
             position = position_dodge(width = .9), size = 3, na.rm = TRUE) +
  scale_shape_manual(name = NULL, values = c(Failed = 4)) +
  scale_fill_manual(values=method_colors) +
  geom_text(
      data = subset(plot_df, method == "SpaDOT"),
      aes(label = paste0("log2FC=", sprintf("%.2f", log2FC)),
          color = log2FC > 0),
      vjust = -0.6,
      size = 4,
      fontface = "bold"
    ) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red"), guide = "none")  + 
  facet_grid(~ Data, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) + 
  labs(
    x = "Transitions",
    y = "Transition Accuracy",
    fill = "Methods"
)


ggsave(file.path(output, 'S42b_real_transition_accuracy.pdf'), plot=p, width=16, height=5.5)
