# conda activate ~/envs/spatialATAC
Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)

# --- calculate concordance index and top-K overlap --- 
c_index <- function(x, y) {
  n <- length(x)
  concordant <- 0
  total <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (x[i] == x[j] || y[i] == y[j]) next
      total <- total + 1
      concordant <- concordant + as.numeric((x[i] - x[j]) * (y[i] - y[j]) > 0)
    }
  }
  concordant / total
}

topk_overlap <- function(x, y, k) {
  # get indices of top-k (highest values)
  topx <- order(x, decreasing = TRUE)[1:k]
  topy <- order(y, decreasing = TRUE)[1:k]
  overlap <- length(intersect(topx, topy)) / k
  return(overlap)
}

final_result_dir = "/net/mulan/home/wenjinma/projects/SpaDOT_analysis/summarize_benchmark/figures"
MouseOrganogenesis_final = read.csv(file.path(final_result_dir, 'ATAC-Visium_noremoval_figures/final_ATAC-Visium_noremoval_domain_metrics.csv'), header=T)
MouseOrganogenesis_final_ARI = MouseOrganogenesis_final %>% filter(index=='ARI', timepoint == 'avg')
MouseOrganogenesis_final_ARI$index = NULL
MouseOrganogenesis_final_ARI$X = NULL
rownames(MouseOrganogenesis_final_ARI) = MouseOrganogenesis_final_ARI$method

ChickenHeart_final = read.csv(file.path(final_result_dir, 'ChickenHeart_figures/final_ChickenHeart_domain_metrics.csv'), header=T)
ChickenHeart_final_ARI = ChickenHeart_final %>% filter(index=='ARI', timepoint == 'avg')
ChickenHeart_final_ARI$index = NULL
ChickenHeart_final_ARI$X = NULL
rownames(ChickenHeart_final_ARI) = ChickenHeart_final_ARI$method

Axolotl_final = read.csv(file.path(final_result_dir, 'Axolotl_figures/final_Axolotl_domain_metrics.csv'), header=T)
Axolotl_final_ARI = Axolotl_final %>% filter(index=='ARI', timepoint == 'avg')
Axolotl_final_ARI$index = NULL
Axolotl_final_ARI$X = NULL
rownames(Axolotl_final_ARI) = Axolotl_final_ARI$method

# load simulation results
simulation_result_dir = '/net/mulan/home/wenjinma/projects/SpaDOT_analysis/summarize_benchmark/figures/simulation_figures/simulation_downstream'
simulation_final = read.csv(file.path(simulation_result_dir, 'final_simulation_domain_metrics.csv'), header=T)
simulation_final_ARI = simulation_final %>% filter(index=='ARI', timepoint == 'avg')
simulation_final_ARI_avg = simulation_final_ARI %>%
    group_by(method, simu_comb) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>% as.data.frame()
methods = unique(simulation_final_ARI_avg$method)
for (simu_comb in unique(simulation_final_ARI_avg$simu_comb)) {
    simu_subset =  simulation_final_ARI_avg[simulation_final_ARI_avg$simu_comb == simu_comb, ]
    rownames(simu_subset) = simu_subset$method
    print(paste0("Spearman correlation for simulation combination: ", simu_comb))

    res_mouse = cor.test(simu_subset[methods, ]$value,
                        MouseOrganogenesis_final_ARI[methods, ]$value,
                        method = "spearman")
    res_chicken = cor.test(simu_subset[methods, ]$value,
                            ChickenHeart_final_ARI[methods, ]$value,
                            method = "spearman")
    res_axolotl = cor.test(simu_subset[methods, ]$value,
                            Axolotl_final_ARI[methods, ]$value,
                            method = "spearman")
    # Print formatted results
    cat("MouseOrganogenesis:  rho =", round(res_mouse$estimate, 3),
        ", p =", signif(res_mouse$p.value, 3), "\n")

    cat("ChickenHeart:        rho =", round(res_chicken$estimate, 3),
        ", p =", signif(res_chicken$p.value, 3), "\n")

    cat("Axolotl:             rho =", round(res_axolotl$estimate, 3),
        ", p =", signif(res_axolotl$p.value, 3), "\n")

    # Printed concordance index
    cat("MouseOrganogenesis Concordance Index:", round(c_index(simu_subset[methods, ]$value,
                                                                MouseOrganogenesis_final_ARI[methods, ]$value), 3), "\n")
    cat("ChickenHeart Concordance Index:", round(c_index(simu_subset[methods, ]$value,
                                                        ChickenHeart_final_ARI[methods, ]$value), 3), "\n")
    cat("Axolotl Concordance Index:", round(c_index(simu_subset[methods, ]$value,
                                                Axolotl_final_ARI[methods, ]$value), 3), "\n")

    # Printed top-5 overlap
    k = 5
    cat("MouseOrganogenesis Top-", k, " Overlap:", round(topk_overlap(simu_subset[methods, ]$value,
                                                                    MouseOrganogenesis_final_ARI[methods, ]$value, k), 5), "\n")
    cat("ChickenHeart Top-", k, " Overlap:", round(topk_overlap(simu_subset[methods, ]$value,
                                                        ChickenHeart_final_ARI[methods, ]$value, k), 5), "\n")
    cat("Axolotl Top-", k, " Overlap:", round(topk_overlap(simu_subset[methods, ]$value,
                                                Axolotl_final_ARI[methods, ]$value, k), 5), "\n")
}




# load simulation prop results
