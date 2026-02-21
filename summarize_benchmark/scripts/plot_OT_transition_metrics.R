library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

datasets = c('MouseEmbryo_Visium_Organogenesis', 'ChickenHeart_Visium_Cardiogenesis', 'AxolotlTelencephalon_Stereoseq_Developmental')
metrics_df_list = list()
for (dataset in datasets) {
    metrics_df_list[[dataset]] = read.csv(file.path(dataset, 'OT_transition_comparisons', 'OT_method_comparison.csv'), header=T)
}
metrics_df = do.call(rbind, metrics_df_list)
metrics_df$method = factor(metrics_df$method, levels=c('SpaDOT', 'DeST-OT', 'Moscot', 'TOAST'))
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
    best_other = max(naive_accuracy[method != "SpaDOT"]),
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


ggplot(plot_df, aes(x = label, y = naive_accuracy, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
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
    legend.text = element_text(size = 12),
  ) + 
  labs(
    x = "Transitions",
    y = "Transition Accuracy",
    fill = "Methods"
)
ggsave('OT_transition_naive_accuracy.pdf', width=9, height=5)


ggplot(metrics_df, aes(x = Day, y = Accuracy, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Data, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +    
  theme_bw() +
  labs(
    x = "Transitions",
    y = "Transition Accuracy (curated)"
)
ggsave('OT_transition_curated_accuracy.pdf', width=6, height=4)
