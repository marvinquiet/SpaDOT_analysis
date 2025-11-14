# conda activate ~/envs/spatialATAC
Sys.setenv(RETICULATE_PYTHON = "/net/mulan/home/wenjinma/envs/spatialATAC/bin/python")
.libPaths("/net/mulan/home/wenjinma/Rlib")
library(reticulate)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)

result_dir = '/net/mulan/home/wenjinma/projects/SpaDOT_analysis/explore_STAMP/ATAC-Visium_noremoval_results'
STAMP_df = read.csv(file.path(result_dir, 'STAMP_5topic_proportions.csv'), header=T, row.names=1)
STAMP_df$pixel_x_bak = STAMP_df$pixel_x
STAMP_df$pixel_y_bak = STAMP_df$pixel_y


STAMP_df = STAMP_df %>%
    mutate(pixel_x = case_when(
        timepoint == 't0' ~ -pixel_y_bak,
        timepoint == 't2' ~ -pixel_y_bak,
        timepoint == 't1' ~ pixel_y_bak
    )) %>%
    mutate(pixel_y = -pixel_x_bak)

STAMP_df = STAMP_df %>%
    mutate(timepoint_original = case_when(
        timepoint == 't0' ~ "E12.5",
        timepoint == 't1' ~ "E13.5",
        timepoint == 't2' ~ "E15.5"
    ))

for (topic in paste0('Topic', 1:5)) {
    p = ggplot(STAMP_df, aes_string(x='pixel_x', y='pixel_y', color=topic)) +
        geom_point(size=1) +
        scale_color_viridis_c() +
        facet_wrap(~timepoint_original, nrow=1, scale='free') +
        theme_bw() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size=14),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10)
        ) +
        ggtitle(paste0('STAMP - ', topic))
    ggsave(filename=file.path(result_dir, paste0('STAMP_ATAC-Visium_noremoval_', topic, '_5topics.png')),
           plot=p,
           width=9, height=3)
}

