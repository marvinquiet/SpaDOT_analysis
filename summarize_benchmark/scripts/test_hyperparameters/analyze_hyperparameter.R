# conda activate ~/envs/spatialATAC
Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)

library(ggplot2)
library(patchwork)
library(ggpubr)
library(glue)
library(dplyr)
library(tidyr)


# ---
quick_search_version = 'quick_search_20260121'
result_dir = '~/projects/SpaDOT/examples'
map_vector = c('ChickenHeart'='ChickenHeart', 'MouseEmbryo'='ATAC-Visium_noremoval', 'AxolotlTelencephalon'='Axolotl')
result_list = list()
for (dataset in c('ChickenHeart', 'MouseEmbryo', 'AxolotlTelencephalon')) {
    quick_search_df = read.csv(glue('{result_dir}/{dataset}_data/hyperparam_results/{quick_search_version}/quick_search_results.csv'), header=T, row.names=1)
    df_prepped = quick_search_df %>%
      mutate(
        mean_ari_clean = if_else(status == "success", mean_ari, NA_real_)
      )
    default_res_df = read.csv(glue('/net/mulan/home/wenjinma/projects/SpaDOT_analysis/summarize_benchmark/figures/{unname(map_vector[dataset])}_figures/final_{unname(map_vector[dataset])}_domain_metrics.csv'), header=T)
    df_prepped = df_prepped[, c('mean_ari_clean'), drop=FALSE]
    default_ari = default_res_df %>%
        filter(method == 'SpaDOT', index == 'ARI', timepoint == 'avg') %>%
        pull(value)
    df_prepped['default', ] = c(default_ari)
    df_prepped$mean_ari_clean = as.numeric(df_prepped$mean_ari_clean)
    df_prepped = df_prepped %>%
          mutate(delta_mean_ari = mean_ari_clean - default_ari)
    df_prepped$hyperparameters = rownames(df_prepped)
    df_prepped_long = df_prepped %>%
        pivot_longer(
            cols = c(mean_ari_clean, delta_mean_ari),
            names_to = "metric",
            values_to = "value"
    )
    df_prepped_long$dataset = dataset
    result_list[[dataset]] = df_prepped_long
}
result_all = bind_rows(result_list)


result_all <- result_all %>%
  mutate(hyperparameters = factor(
    hyperparameters,
    levels = c("default", setdiff(unique(hyperparameters), "default"))
  ))

ggplot(result_all %>% filter(metric == "mean_ari_clean"),
       aes(x = hyperparameters, y = dataset, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = "Mean ARI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Mean_ARI_per_dataset.pdf', height=2, width=8)


ggplot(result_all %>% filter(metric == "delta_mean_ari", hyperparameters != 'default'),
       aes(x = hyperparameters, y = dataset, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "white", name = "Delta ARI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Delta_ARI_per_dataset.pdf', height=2, width=8)







































































































































































































































































































































































































































