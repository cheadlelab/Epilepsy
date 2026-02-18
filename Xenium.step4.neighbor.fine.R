# https://divingintogeneticsandgenomics.com/post/neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data-in-seurat-and-bioconductor/
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ
rm(list = ls(all = TRUE))

#load packages
library(Seurat)
library(dplyr)
library(jsonlite)
library(dplyr)
library(arrow)
library(spacexr)
library(ggplot2)
library(dplyr)

source("//grid/cheadle_home/qianyu/Epilepsy/Xenium/ReadXenium2.R")
# source("/grid/cheadle/home/qianyu/Epilepsy/Xenium/ReadXenium2.R")

save.path <- "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/"


for (i in 1:11) {
  load(paste0(save.path, "step1.sample_", i, ".RData"))
  sample_name <- paste0("sample_", i)
  xenium.obj <- get(sample_name)
  
  xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "predicted.celltype")
  xenium.obj@meta.data$layer <- xenium.obj@meta.data$niches
  
  
  assign(sample_name, xenium.obj)
}
sample_list <- list(
  "sample1" = sample_1,
  "sample2" = sample_2,
  "sample3" = sample_3,
  "sample4" = sample_4,
  "sample5" = sample_5,
  "sample6" = sample_6,
  "sample7" = sample_7,
  "sample8" = sample_8,
  "sample9" = sample_9,
  "sample10" = sample_10,
  "sample11" = sample_11
)




# =========================
# Xenium neighborhood analysis: Focal vs Distal (no layers)
# Full script (orders fixed on both axes, remove "Other" and NA)
# =========================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# 2) Per-sample summary
summary_list <- setNames(vector("list", length(sample_list)), names(sample_list))
for (nm in names(sample_list)) {
  x <- sample_list[[nm]]
  
  niche_mat <- GetAssayData(x, assay = "niche", slot = "data")
  cs <- colSums(niche_mat)
  keep <- which(cs > 0)
  if (!length(keep)) next
  niche_mat  <- niche_mat[, keep, drop = FALSE]
  cs         <- cs[keep]
  niche_prop <- sweep(niche_mat, 2, cs, "/")
  
  meta <- x@meta.data[colnames(niche_prop),
                      c("layer","predicted.celltype","orig.ident"),
                      drop = FALSE]
  meta$cell <- rownames(meta)
  
  prop_df <- as.data.frame(t(niche_prop))
  prop_df$cell <- rownames(prop_df)
  
  long_df <- prop_df %>%
    pivot_longer(-cell, names_to = "neighbor_celltype", values_to = "prop") %>%
    left_join(meta, by = "cell")
  
  summary_df <- long_df %>%
    group_by(layer, predicted.celltype, neighbor_celltype) %>%
    summarise(mean_prop = mean(prop, na.rm = TRUE),
              n_cells   = dplyr::n_distinct(cell),
              .groups   = "drop")
  
  oi <- unique(x@meta.data$orig.ident)[1]
  condition <- ifelse(grepl("C", oi), "Distal", "Focal")
  patient   <- as.integer(gsub("[^0-9].*", "", oi))
  
  summary_df$orig.ident <- oi
  summary_df$condition  <- condition
  summary_df$patient    <- patient
  summary_df <- relocate(summary_df, orig.ident, condition, patient, .before = layer)
  
  summary_list[[nm]] <- summary_df
}

# 3) Combine
summary_all <- bind_rows(summary_list)

# 4) Patient-level (drop "Other" on BOTH axes)
patient_df <- summary_all %>%
  filter(!is.na(patient), !is.na(condition),
         predicted.celltype != "Other",
         neighbor_celltype != "Other") %>%
  group_by(patient, condition, predicted.celltype, neighbor_celltype) %>%
  summarise(val = weighted.mean(mean_prop, w = n_cells, na.rm = TRUE),
            .groups = "drop")

# 5) All combos present
comb <- tidyr::expand_grid(
  predicted.celltype    = unique(patient_df$predicted.celltype),
  neighbor_celltype = unique(patient_df$neighbor_celltype)
)

# 6) F/D ratio + t-test
res_df <- pmap_dfr(
  comb,
  function(predicted.celltype, neighbor_celltype) {
    sub <- patient_df %>%
      filter(predicted.celltype == !!predicted.celltype,
             neighbor_celltype == !!neighbor_celltype)
    
    foc <- sub %>% filter(condition == "Focal")  %>% distinct(patient, .keep_all = TRUE) %>% pull(val)
    dis <- sub %>% filter(condition == "Distal") %>% distinct(patient, .keep_all = TRUE) %>% pull(val)
    
    eps <- 1e-6
    mF <- mean(foc, na.rm = TRUE); mD <- mean(dis, na.rm = TRUE)
    effect <- (mF + eps) / (mD + eps)
    
    pval <- NA_real_
    if (length(foc) >= 2 && length(dis) >= 2) {
      pval <- tryCatch(stats::t.test(foc, dis)$p.value, error = function(e) NA_real_)
    }
    
    tibble(predicted.celltype, neighbor_celltype, effect, p = pval)
  }
) %>%
  mutate(sig = ifelse(!is.na(p) & p < 0.05, "*", ""))

# 7) Fix BOTH axis orders and drop NA columns/rows
cell_order <- c(
  # Excitatory neurons
  "Exc_L23IT", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L234IT_ERBB4",
  "Exc_L5IT_GRIN3A", "Exc_L5IT_GRIN3A-", "Exc_L5ET", "Exc_L56IT_CAR3", 
  "Exc_L56NP", "Exc_L6IT", "Exc_L6CT", "Exc_L6b",
  
  # Inhibitory neurons
  "Inh_PVALB", "Inh_SST", "Inh_VIP", "Inh_LAMP5", "Inh_LAMP5_LHX6", 
  "Inh_CXCL14", "Chandelier",
  
  # Glia
  "Astrocyte", "Microglia", "OPC", "Olg",
  
  # Immune / vascular
  "Macrophage", "NK-T", "Vascular"
)

res_df$neighbor_celltype <- gsub("-(?!$)", "_", res_df$neighbor_celltype, perl = TRUE)
res_df$neighbor_celltype <- gsub("^NK_T$", "NK-T", res_df$neighbor_celltype)


res_df <- res_df %>%
  filter(!is.na(predicted.celltype), !is.na(neighbor_celltype)) %>%
  mutate(
    predicted.celltype    = factor(predicted.celltype,    levels = cell_order),
    neighbor_celltype = factor(neighbor_celltype, levels = cell_order)
  )
res_df$log_effect <- log2(res_df$effect)

# 8) Heatmap
p_global <- ggplot(res_df, aes(x = neighbor_celltype, y = predicted.celltype, fill = log_effect)) +
  geom_tile() +
  geom_text(aes(label = sig), vjust = 0.5, size = 2) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,                      # log(1) = 0
    name = "log2(Focal/Distal)"
  ) +
  labs(
    title = "log2(Focal / Distal) ratio",
    x = "Neighbor cell type", y = "Cell type"
  ) +
  theme_minimal(base_size = 6) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y  = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    plot.title   = element_text(size = 6),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6),
    panel.grid   = element_blank()
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 0.3,
      barheight = 2
    )
  )

p_global