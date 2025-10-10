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
metrics_df$Data = factor(metrics_df$Data, levels=c('Mouse Embryo', 'Chicken Heart', 'Axolotl Telencephalon'))

ggplot(metrics_df, aes(x = Day, y = naive_accuracy, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Data, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +    
  theme_bw() +
  labs(
    x = "Transitions",
    y = "Transition Accuracy"
)
ggsave('OT_transition_naive_accuracy.pdf', width=6, height=4)


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