rm(list = ls(all = TRUE))

# Load necessary libraries
library(Seurat)
library(dplyr)
library(egg)
library(ggplot2)
library(scales)

##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))
load(paste0(file_path, "/step2.1.Comparisons.sc.RData"))

# Convert to character if necessary
test@meta.data$Annotation.fine

# Immediate Early Genes (IEG)
IEG <- c("FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "ARC", "NR4A1", "NR4A2", "NPAS4", "ATF3", "DUSP1", "DUSP6")
HSP <- c("HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1")
genes_of_interest <- HSP

# Define cell types and conditions
cell_types <- c("Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1", 
                "Exc_L5IT_GRIN3A", "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L6IT",
                "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier",
                "Astrocyte", "Microglia", "OPC", "Olg")
cell_types <- rev(cell_types)

# ---------- Pseudobulk aggregation (same as before) ----------
pseudo_test <- AggregateExpression(
  test,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("Annotation.fine", "condition")
)

# Unify annotation names (avoid underscore/dash issues)
pseudo_test$Annotation.fine <- gsub("-", "_", pseudo_test$Annotation.fine)
pseudo_test$Annotation.fine <- gsub("Exc_L5IT_GRIN3A_", "Exc_L5IT_GRIN3A-", pseudo_test$Annotation.fine)

# Define group order
conditions <- c("Distal", "Stimulated", "Focal")
group_levels <- c(paste0(cell_types, "_Distal"),
                  paste0(cell_types, "_Stimulated"),
                  paste0(cell_types, "_Focal"))

meta_pb <- pseudo_test@meta.data
meta_pb$group <- paste(meta_pb$Annotation.fine, meta_pb$condition, sep = "_")
meta_pb$group <- factor(meta_pb$group, levels = group_levels)
meta_pb$Condition <- factor(sub(".*_(Focal|Stimulated|Distal)$", "\\1", meta_pb$group),
                            levels = conditions)

# ---------- Extract pseudobulk matrix ----------
# Use log-normalized "data" (Seurat default is natural log after NormalizeData)
mat <- Seurat::GetAssayData(pseudo_test, assay = "RNA", slot = "data")
genes_keep <- intersect(genes_of_interest, rownames(mat))
mat <- mat[genes_keep, , drop = FALSE]

# Reorder columns according to defined group order
ord <- order(meta_pb$group)
mat <- mat[, ord, drop = FALSE]
meta_pb <- meta_pb[ord, , drop = FALSE]

# ---------- Compute log2FC vs Distal per cell type ----------
# On the log (natural) scale: lnFC = ln(expr_cond) - ln(expr_distal); convert to log2 by dividing ln(2)
ln2 <- log(2)
groups_vec <- as.character(meta_pb$group)
celltype_of_group <- sub("_(Focal|Stimulated|Distal)$", "", groups_vec)
condition_of_group <- as.character(meta_pb$Condition)

mat_l2fc <- mat * 0  # initialize

for (ct in unique(celltype_of_group)) {
  # columns for this cell type
  cols_ct <- which(celltype_of_group == ct)
  # distal column for this cell type
  col_dist <- which(groups_vec == paste0(ct, "_Distal"))
  if (length(col_dist) != 1) {
    # if baseline missing, set NA
    mat_l2fc[, cols_ct] <- NA_real_
    next
  }
  # lnFC relative to distal
  lnFC <- sweep(mat[, cols_ct, drop = FALSE], 1, mat[, col_dist, drop = FALSE], FUN = "-")
  # convert to log2FC
  mat_l2fc[, cols_ct] <- lnFC / ln2
}

# For the Distal baseline itself, log2FC is 0 by definition
# (already true due to difference against itself)

# Optional clipping for visualization
mat_l2fc[mat_l2fc >  2.5] <-  2.5
mat_l2fc[mat_l2fc < -2.5] <- -2.5

# ---------- Significance mask (unchanged logic) ----------
sig_hits <- dplyr::tibble()
comparisons <- c("Focal_vs_Distal", "Stimulated_vs_Distal")
for (celltype in cell_types) {
  for (comparison in comparisons) {
    deg <- final_results[[celltype]][[comparison]]$deg
    if (!is.null(deg) && nrow(deg) > 0) {
      deg$symbol <- rownames(deg)
      add <- deg %>%
        dplyr::filter(p_val_adj < 0.05,
                      abs(avg_log2FC) > log2(1.5),
                      symbol %in% genes_keep) %>%
        dplyr::mutate(
          features.plot = symbol,
          condition = ifelse(grepl("Focal", comparison), "Focal", "Stimulated"),
          group = paste(celltype, condition, sep = "_")
        ) %>%
        dplyr::select(features.plot, group) %>%
        dplyr::distinct()
      sig_hits <- dplyr::bind_rows(sig_hits, add)
    }
  }
}
sig_hits <- dplyr::distinct(sig_hits)

# ---------- Long format for ggplot ----------
library(tidyr)
df_hm <- as.data.frame(mat_l2fc)
df_hm$Gene <- rownames(df_hm)
df_hm <- tidyr::pivot_longer(df_hm, cols = -Gene, names_to = "Col", values_to = "l2fc")

df_hm$group <- meta_pb$group[match(df_hm$Col, colnames(mat_l2fc))]
df_hm$Condition <- meta_pb$Condition[match(df_hm$Col, colnames(mat_l2fc))]

df_hm$Gene   <- factor(df_hm$Gene, levels = genes_keep)
df_hm$group  <- factor(df_hm$group, levels = group_levels)
df_hm$Condition <- factor(df_hm$Condition, levels = conditions)

df_hm$Significant <- paste(df_hm$Gene, df_hm$group) %in% paste(sig_hits$features.plot, sig_hits$group)

# ---------- Heatmap (colored by log2FC vs Distal) ----------
expression_colors <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
)

p_both <- ggplot(df_hm, aes(x = Gene, y = group, fill = l2fc)) +
  geom_tile() +
  geom_text(
    data = subset(df_hm, Significant),
    aes(label = "*"),
    size = 2.2, color = "black", na.rm = TRUE
  ) +
  scale_fill_gradientn(
    colours = expression_colors,
    values = scales::rescale(c(-2, -1, 0, 1, 2)),
    limits = c(-2.5, 2.5),
    name = "log2FC vs Distal"
  ) +
  scale_x_discrete(position = "top") +
  facet_grid(rows = vars(Condition), scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, face = "italic", color = "black", size = 7),
    axis.text.y = element_text(color = "black", size = 6),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.margin = margin(3, 3, 3, 3, "pt")
  ) +
  labs(x = "Genes", y = "Cell type × Condition")



########################################################################################
# ---------- Keep only "Stimulated" ----------
df_plot <- subset(df_hm, Condition == "Stimulated")

# Remove the "_Stimulated" suffix from y-axis labels
df_plot$group_clean <- sub("_Stimulated$", "", as.character(df_plot$group))
df_plot$group_clean <- factor(df_plot$group_clean, levels = cell_types)

pp <- ggplot(df_plot, aes(x = Gene, y = group_clean, fill = l2fc)) +
  geom_tile() +
  geom_text(
    data = subset(df_plot, Significant),
    aes(label = "*", y = group_clean),
    size = 2.2, color = "black", na.rm = TRUE
  ) +
  scale_fill_gradientn(
    colours = expression_colors,
    values = scales::rescale(c(-2, -1, 0, 1, 2)),
    limits = c(-2.5, 2.5),
    name = "log2FC vs Distal"
  ) +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, face = "italic", color = "black", size = 7),
    axis.text.y = element_text(color = "black", size = 6),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 7),
    panel.grid = element_blank(),
    plot.margin = margin(3, 3, 3, 3, "pt")
  ) +
  labs(x = "Genes", y = "Cell type")









# ---------- Save PDF ----------
pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V3/Figure.Stimulated.HSP.pdf", width = 2.5, height = 2.5)
print(pp)
dev.off()
