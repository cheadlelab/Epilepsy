rm(list = ls(all = TRUE))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(ggpubr)   # for significance annotation on the plot
})

## ------------------ Load ------------------
file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output/Blood"
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"

load(paste0(file_path, "/step1.2.umap.fine.RData"))
Idents(test) <- "condition"

## ------------------ Filter & table ------------------
# Look at distribution before filtering
table(test$Annotation.fine, test$patient)

# Remove cell types whose names contain "/"
test <- subset(test, subset = !grepl("/", Annotation.fine))

# Look at distribution after filtering
table(test$Annotation.fine, test$patient)

# Cross table of cell types vs patients
# Normalize by columns (per patient) to get percentages
tab      <- table(test$Annotation.fine, test$patient)
tab_prop <- prop.table(tab, margin = 2) * 100  # percentage (0-100)

## ------------------ Choose celltype ------------------
celltype <- "Treg"   # <<< change this to the cell type of interest

stopifnot(celltype %in% rownames(tab_prop))

# Extract the percentage of this cell type for each patient
pct_vec <- as.numeric(tab_prop[celltype, , drop = TRUE])

df_box <- data.frame(
  Patient   = colnames(tab_prop),
  Percent   = pct_vec,
  Condition = ifelse(grepl("^Epileptic", colnames(tab_prop)), "Epileptic", "Healthy"),
  stringsAsFactors = FALSE
)

# Remove possible NA values (in case some patients had 0 cells)
df_box <- df_box %>% filter(!is.na(Percent))

## ------------------ Plot with significance ------------------
# Comparison between Healthy and Epileptic
cmp <- list(c("Healthy", "Epileptic"))

p <- ggplot(df_box, aes(x = Condition, y = Percent)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.75, size = 2) +
  labs(
    title = paste0("Per-patient percent of ", celltype),
    x = NULL, y = "Percent of cells in patient (%)"
  ) +
  theme_bw(base_size = 12) +
  stat_compare_means(
    comparisons = cmp,
    method = "wilcox.test",
    label = "p.signif"       # show significance stars; use "p.format" to show exact p-value
  )

print(p)

## ------------------ Wilcoxon test (text output) ------------------
if (all(table(df_box$Condition) >= 2)) {
  w <- wilcox.test(Percent ~ Condition, data = df_box, exact = FALSE)
  cat("\nWilcoxon test on per-patient percent (", celltype, "):\n", sep = "")
  print(w)
} else {
  cat("\n[Warn] Not enough patients per group for Wilcoxon test.\n")
}

## ------------------ Optional: save the figure ------------------
# ggsave(paste0("box_percent_", celltype, ".png"), p, width = 6, height = 4, dpi = 300)
