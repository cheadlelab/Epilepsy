rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(ggplot2)

##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"
# load(paste0(file_path, "/step2.1.Comparisons.sc.downsampling.RData"))
load(paste0(file_path, "/step2.1.DEG.singlecell.RData"))


celltypes <- names(patient_regressed)

cell_types <- c(
  "CD14_Mono",
  "CD16_Mono",
  "cDC1",
  "cDC2",
  "pDC",
  "EMP",
  "Plasmablast",
  "B_naive",
  "B_memory",
  "B_intermediate",
  "CD4_Naive",
  "CD8_Naive",
  "CD4_Memory",
  "CD8_Memory",
  "CD4_CTL",
  "Treg",
  "MAIT",
  "MKI67+_NKT",
  "gdT",
  "NK",
  "NK_CD56bright",
  "Platelets"
)

plot_df <- lapply(celltypes, function(ct) {
  res <- patient_regressed[[ct]]
  if (!is.null(res) && nrow(res) > 0) {
    res <- res %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > log2(1.2)) %>%
      mutate(Direction = ifelse(avg_log2FC > 0, "Up", "Down"))
    up_n <- sum(res$Direction == "Up")
    down_n <- sum(res$Direction == "Down")
  } else {
    up_n <- 0; down_n <- 0
  }
  tibble(CellType = ct,
         Up = up_n,
         Down = down_n)
}) %>% bind_rows()

plot_df_long <- plot_df %>%
  tidyr::pivot_longer(cols = c("Up", "Down"), names_to = "Direction", values_to = "Count") %>%
  filter(CellType %in% cell_types) %>%
  mutate(CellType = factor(CellType, levels = cell_types),
         Direction = factor(Direction, levels = c("Up", "Down")))

dir_colors <- c("Up" = "#E41A1C", "Down" = "#377EB8")

p <- ggplot(plot_df_long, aes(x = CellType, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = dir_colors) +
  theme_bw() +
  labs(title = "DEG Count (PBMC)",
       x = "Cell Type", y = "DEG Count") +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    plot.title = element_text(hjust = 0.5, size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.position = "top"
  )
print(p)



pdf(file = "C:/Users/thech/Desktop/Figure.tmp.pdf", width = 3.5, height = 2.8)
print(p)
dev.off() # Close the PDF device
