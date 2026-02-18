rm(list = ls(all = TRUE))

#########################
## 0. Load packages and merged object
#########################

library(Seurat)
library(dplyr)
library(jsonlite)
library(arrow)
library(spacexr)
library(ggplot2)
library(scales)
library(tibble)
library(ggpubr)

options(future.globals.maxSize = 128 * 1024^3)

# Detect OS and load merged Xenium object
os_name <- Sys.info()[["sysname"]]
save.path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/qianyu/Epilepsy/Xenium/output/",
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/",
                    stop(sprintf("Unsupported OS: %s", os_name))
)

# This should load a Seurat object named "test"
load(paste0(save.path, "step1.merged.RData"))


#########################
## 1. Define target cell type to analyze
#########################

celltype <- c("Microglia")
# You can change this to any predicted.celltype, e.g.:
# celltype <- "Astrocyte"
# celltype <- "Olg"
# celltype <- "OPC"


#########################
## 2. Define Focal vs Distal information by sample
#########################

sample_info <- tibble(
  orig.ident = c("1C SA", "4E", "1C SB", "5C", "2E",
                 "3E SB", "1E", "4C", "5E", "2C", "3E SA"),
  exp.ident  = c("Non-epileptiform", "Epileptiform", "Non-epileptiform",
                 "Non-epileptiform", "Epileptiform", "Epileptiform",
                 "Epileptiform", "Non-epileptiform", "Epileptiform",
                 "Non-epileptiform", "Epileptiform")
) %>%
  mutate(region_type = if_else(exp.ident == "Epileptiform", "Focal", "Distal"))


#########################
## 3. Build meta table and keep only selected samples
#########################

meta <- test@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  left_join(sample_info, by = "orig.ident")

# Keep only these samples: 1E, 4C, 2C, 1C SA, 4E, 2E, 3E SB
keep_samples <- c("1E", "4C", "2C", "1C SA", "4E", "2E", "3E SB")
meta <- meta %>% filter(orig.ident %in% keep_samples)

# (Optional) quick check
# table(meta$orig.ident, meta$region_type)


#########################
## 4. Compute per-sample cell-type proportions per layer
#########################
# For each sample (orig.ident) and each anatomical layer, compute:
# - total number of cells
# - number of cells belonging to the target celltype
# - proportion = celltype_cells / total_cells

layer_by_sample <- meta %>%
  filter(!is.na(layer)) %>%
  group_by(orig.ident, region_type, layer) %>%
  summarise(
    total_cells    = n(),
    celltype_cells = sum(predicted.celltype %in% celltype),
    prop_celltype  = celltype_cells / total_cells,
    .groups = "drop"
  )

# Exclude "Other" layer
layer_df <- layer_by_sample %>%
  filter(layer != "Other")


#########################
## 5. Statistical testing (Wilcoxon rank-sum test)
#########################
# Compare Focal vs Distal within each layer (unpaired Wilcoxon)

p_to_sig <- function(p){
  if (p <= 0.0001) return("****")
  else if (p <= 0.001) return("***")
  else if (p <= 0.01) return("**")
  else if (p <= 0.05) return("*")
  else return("ns")
}

layer_stats <- layer_df %>%
  group_by(layer) %>%
  summarise(
    p_value      = wilcox.test(prop_celltype ~ region_type)$p.value,
    focal_mean   = mean(prop_celltype[region_type == "Focal"]),
    distal_mean  = mean(prop_celltype[region_type == "Distal"]),
    significance = p_to_sig(p_value),
    .groups = "drop"
  )

layer_stats   # summary table of p-values and means


#########################
## 6. Visualization: Focal vs Distal for each layer
#########################

# Optional: enforce a consistent layer order
layer_df$layer <- factor(
  layer_df$layer,
  levels = c("Layer-Exc23", "Layer-Exc456", "Layer-NVU", "Layer-Olg")
)

# Basic boxplot
p_box <- ggplot(layer_df,
                aes(x = layer, y = prop_celltype,
                    fill = region_type, color = region_type)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.9) +
  ylab(paste0(celltype, " proportion within layer")) +
  xlab("Layer") +
  ggtitle(paste0(celltype, " proportion across layers (Focal vs Distal)")) +
  theme_classic() +
  scale_color_manual(values = c("Focal" = "#e31a1c",   # red-ish
                                "Distal" = "#1f78b4")) +  # blue-ish
  scale_fill_manual(values  = c("Focal" = "#e31a1c",
                                "Distal" = "#1f78b4"))
p_box
#########################
## 7. Visualization with p-values on the plot
#########################
library(ggpubr)

dodge_width <- 1

layer_df <- layer_df %>%
  mutate(
    region_type = factor(region_type, levels = c("Distal", "Focal")),
    layer = factor(layer, levels = c("Layer-Exc23", "Layer-Exc456", "Layer-NVU", "Layer-Olg"))
  )

cols <- c("Distal" = "#1cbdc3", "Focal" = "#f3766e")


p_box_sig <- ggplot(layer_df,
                    aes(x = layer, y = prop_celltype,
                        fill = region_type, color = region_type)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA,
               position = position_dodge(width = 1)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 1),
             size = 2, alpha = 0.9) +
  ylab(paste0(celltype, " proportion within layer")) +
  xlab("Layer") +
  ggtitle(paste0(celltype, " proportion across layers (Focal vs Distal)")) +
  theme_classic() +
  ggpubr::stat_compare_means(
    aes(group = region_type),
    method = "wilcox.test",
    label  = "p.format"
  ) +
  scale_fill_manual(values = cols, drop = FALSE) +
  scale_color_manual(values = cols, drop = FALSE)


p_box_sig # 6*3



layer_df <- layer_df %>%
  filter(layer == "Layer-Exc23")

p_box_sig <- ggplot(layer_df,
                    aes(x = layer, y = prop_celltype,
                        fill = region_type, color = region_type)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA,
               position = position_dodge(width = 1)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 1),
             size = 2, alpha = 0.9) +
  ylab(paste0(celltype, " proportion within layer")) +
  xlab("Layer") +
  ggtitle(paste0(celltype, " proportion in Layer-Exc23 (Focal vs Distal)")) +
  theme_classic() +
  ggpubr::stat_compare_means(
    aes(group = region_type),
    method = "wilcox.test",
    label  = "p.format"
  ) +
  scale_fill_manual(values = cols, drop = FALSE) +
  scale_color_manual(values = cols, drop = FALSE)

p_box_sig # 6*3





