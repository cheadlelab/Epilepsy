rm(list = ls(all = TRUE))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(dplyr)
  library(Matrix)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(scales)
  library(ggplot2)
  library(ggrepel)
})
# source("//grid/cheadle_home/qianyu/Epilepsy/figures/manuscript/VolcanoPlot.R")

##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))

exc_cell_types <- c("Exc_L23IT", "Exc_L4IT", "Exc_L5IT_GRIN3A", "Exc_L6IT")
celltypes <- exc_cell_types

fc_thr <- log2(1.5)
eps    <- 1e-150

deg_list <- list()

for (ct in celltypes) {
  message("DE for: ", ct)
  
  deg <- final_results[[ct]][["Focal_vs_Distal"]][["deg"]]
  deg$gene  <- rownames(deg)
  deg$Group <- ct
  deg_list[[ct]] <- deg
}

diff_sc <- bind_rows(deg_list)

diff_sc <- diff_sc %>%
  mutate(
    log10fdr   = -log10(p_val_adj + eps),
    abs_log2FC = abs(avg_log2FC)
  )
top10 <- diff_sc %>%
  group_by(Group) %>%
  slice_max(order_by = abs_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

mycol <- hue_pal()(length(celltypes))

p <- ggplot() +
  geom_point(data = diff_sc,
             aes(x = avg_log2FC, y = log10fdr),
             size = 0.7, color = 'grey70', alpha = 0.6) +
  geom_point(data = top10,
             aes(x = avg_log2FC, y = log10fdr, color = Group),
             size = 1.2) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), size = 0.3, color = "grey50", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.3, color = "grey50") +
  coord_flip() +
  facet_grid(. ~ Group, scales = "free") +
  scale_color_manual(values = mycol[seq_along(unique(top10$Group))]) +
  labs(x = "avg_log2FC (epileptic vs healthy)", y = expression(-log[10]("adjusted p"))) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = 'none',
    panel.grid      = element_blank(),
    axis.text.x     = element_text(angle = 45, vjust = 0.8),
    strip.text.x    = element_text(size = 10, face = 'bold'),
    plot.margin     = margin(5,5,5,5)
  )

pt6 <- 6
pt6_geom <- 6/2.845

p1 <- p +
  theme(
    text         = element_text(size = pt6),
    axis.text    = element_text(size = pt6),
    axis.title   = element_text(size = pt6),
    strip.text.x = element_text(size = pt6, face = "bold"),
    legend.text  = element_text(size = pt6),
    legend.title = element_text(size = pt6)
  ) +
  geom_text_repel(
    data = top10,
    aes(x = avg_log2FC, y = log10fdr, label = gene, color = Group),
    size = pt6_geom,
    seed = 233,
    direction = "y",
    box.padding = 0.15,
    point.padding = 0.1,
    min.segment.length = 0,
    segment.size = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  )

pdf("//grid/cheadle_home/qianyu/Epilepsy/output/Figure.Focal.volcanoplot.pdf", width=3.8, height=2)
p1
dev.off()

