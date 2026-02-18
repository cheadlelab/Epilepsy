# =========================
# Xenium neighborhood analysis: Focal vs Distal (no layers)
# Full script (orders fixed on both axes, remove "Other" and NA)
# =========================

rm(list = ls(all = TRUE))

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# 1) Load Seurat list object
load("//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/step1.list.RData")
# expecting: sample_list

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
                      c("layer","celltype.corse","orig.ident"),
                      drop = FALSE]
  meta$cell <- rownames(meta)
  
  prop_df <- as.data.frame(t(niche_prop))
  prop_df$cell <- rownames(prop_df)
  
  long_df <- prop_df %>%
    pivot_longer(-cell, names_to = "neighbor_celltype", values_to = "prop") %>%
    left_join(meta, by = "cell")
  
  summary_df <- long_df %>%
    group_by(layer, celltype.corse, neighbor_celltype) %>%
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
         celltype.corse != "Other",
         neighbor_celltype != "Other") %>%
  group_by(patient, condition, celltype.corse, neighbor_celltype) %>%
  summarise(val = weighted.mean(mean_prop, w = n_cells, na.rm = TRUE),
            .groups = "drop")

# 5) All combos present
comb <- tidyr::expand_grid(
  celltype.corse    = unique(patient_df$celltype.corse),
  neighbor_celltype = unique(patient_df$neighbor_celltype)
)

# 6) F/D ratio + t-test
res_df <- pmap_dfr(
  comb,
  function(celltype.corse, neighbor_celltype) {
    sub <- patient_df %>%
      filter(celltype.corse == !!celltype.corse,
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
    
    tibble(celltype.corse, neighbor_celltype, effect, p = pval)
  }
) %>%
  mutate(sig = ifelse(!is.na(p) & p < 0.05, "*", ""))

# 7) Fix BOTH axis orders and drop NA columns/rows
cell_order <- c("Exc_L23","Exc_L456",
                "Inh",
                "Astrocyte","Microglia","OPC","Olg",
                "Macrophage",
                "Vascular")

res_df$neighbor_celltype <- gsub("-", "_", res_df$neighbor_celltype)



res_df <- res_df %>%
  filter(!is.na(celltype.corse), !is.na(neighbor_celltype)) %>%
  mutate(
    celltype.corse    = factor(celltype.corse,    levels = cell_order),
    neighbor_celltype = factor(neighbor_celltype, levels = cell_order)
  )

# 8) Heatmap
p_global <- ggplot(res_df, aes(x = neighbor_celltype, y = celltype.corse, fill = effect)) +
  geom_tile() +
  geom_text(aes(label = sig), vjust = 0.5, size = 2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 1, name = "Focal/Distal") +
  labs(title = "Focal / Distal ratio",
       x = "Neighbor cell type", y = "Cell type") +
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

pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V3/Figure.Xenium.pdf", width = 2.3, height = 1.6)
p_global
dev.off() # Close the PDF device