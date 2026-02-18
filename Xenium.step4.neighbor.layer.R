rm(list = ls(all = TRUE))

library(Seurat)
library(dplyr)
library(jsonlite)
library(arrow)
library(spacexr)
library(ggplot2)
library(tidyr)

load("//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/step1.list.RData")

# Prepare an empty list to hold results for each sample
summary_list <- setNames(vector("list", length(sample_list)), names(sample_list))

for (nm in names(sample_list)) {
  xenium.obj <- sample_list[[nm]]

  # Inspect a small slice of the niche assay (first 5 cells)
  mat <- GetAssayData(xenium.obj, assay = "niche", slot = "data")[, 1:5]
  as.data.frame(mat[apply(mat, 1, sum) > 0, ])
  
  # 1) Extract the niche assay (rows = neighbor cell types, cols = cells)
  niche_mat <- GetAssayData(xenium.obj, assay = "niche", slot = "data")
  
  # 2) Normalize each cell's neighborhood counts to proportions
  col_sums <- colSums(niche_mat)
  keep_cols <- which(col_sums > 0)                        # drop cells with zero neighbors
  niche_mat <- niche_mat[, keep_cols, drop = FALSE]
  col_sums  <- col_sums[keep_cols]
  niche_prop <- sweep(niche_mat, 2, col_sums, "/")        # divide each column by its total
  
  # 3) Add metadata (layer and celltype.corse) and reshape to long format
  meta <- xenium.obj@meta.data[colnames(niche_prop), c("layer","celltype.corse"), drop = FALSE]
  meta$cell <- rownames(meta)
  
  prop_df <- as.data.frame(t(niche_prop))                 # transpose: now rows = cells
  prop_df$cell <- rownames(prop_df)
  
  long_df <- prop_df |>
    pivot_longer(cols = -cell,
                 names_to = "neighbor_celltype",
                 values_to = "prop") |>
    left_join(meta, by = "cell")
  
  # 4) Summarize by layer × celltype.corse × neighbor_celltype
  #    -> mean proportion across cells and number of contributing cells
  summary_df <- long_df |>
    group_by(layer, celltype.corse, neighbor_celltype) |>
    summarise(
      mean_prop = mean(prop, na.rm = TRUE),
      n_cells   = dplyr::n_distinct(cell),
      .groups   = "drop"
    ) |>
    arrange(layer, celltype.corse, desc(mean_prop))
  
  head(summary_df)                       # quick peek at results
  summary_list[[nm]] <- summary_df       # save into list for this sample
}






## ---- Append orig.ident / condition / patient ----

parse_orig_ident <- function(x) {
  m <- regexec("^\\s*(\\d+)([CE])", x)
  g <- regmatches(x, m)[[1]]
  if (length(g) == 0) return(list(patient = NA_integer_, condition = NA_character_))
  patient <- as.integer(g[2])
  condition <- ifelse(g[3] == "C", "Distal", "Focal")
  list(patient = patient, condition = condition)
}

for (nm in names(summary_list)) {
  df <- summary_list[[nm]]
  if (is.null(df) || nrow(df) == 0) next
  
  # pick orig.ident (first unique value if multiple)
  oi_vals <- unique(sample_list[[nm]]@meta.data$orig.ident)
  oi_vals <- oi_vals[!is.na(oi_vals)]
  oi_val  <- if (length(oi_vals)) oi_vals[1] else NA_character_
  
  info <- parse_orig_ident(oi_val)
  
  df$orig.ident <- oi_val
  df$condition  <- info$condition
  df$patient    <- info$patient
  
  # put them before layer
  df <- dplyr::relocate(df, orig.ident, condition, patient, .before = layer)
  
  summary_list[[nm]] <- df
}

# Combine into one big table
summary_all <- dplyr::bind_rows(summary_list)



###########################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# 1) Collapse sample-level summary_all to patient-level means (weighted by n_cells)
patient_df <- summary_all %>%
  filter(!is.na(patient), !is.na(condition)) %>%
  group_by(patient, condition, layer, celltype.corse, neighbor_celltype) %>%
  summarise(
    val = weighted.mean(mean_prop, w = n_cells, na.rm = TRUE),  # patient-level mean (weighted)
    .groups = "drop"
  )

# 2) For each layer: compute F/D ratio and significance across patients
layers <- sort(unique(patient_df$layer))
plots  <- setNames(vector("list", length(layers)), layers)

for (lay in layers) {
  dat_l <- patient_df %>% filter(layer == lay)
  
  # All (focal celltype, neighbor celltype) combinations observed in this layer
  comb <- tidyr::expand_grid(
    celltype.corse    = sort(unique(dat_l$celltype.corse)),
    neighbor_celltype = sort(unique(dat_l$neighbor_celltype))
  )
  
  # For each combination, collect patient-level vectors for Focal and Distal,
  # compute F/D ratio (using epsilon to avoid divide-by-zero),
  # and test significance across patients (Welch t-test).
  res <- pmap_dfr(
    comb,
    function(celltype.corse, neighbor_celltype) {
      sub <- dat_l %>%
        filter(celltype.corse == !!celltype.corse,
               neighbor_celltype == !!neighbor_celltype)
      
      foc <- sub %>% filter(condition == "Focal")  %>% distinct(patient, .keep_all = TRUE) %>% pull(val)
      dis <- sub %>% filter(condition == "Distal") %>% distinct(patient, .keep_all = TRUE) %>% pull(val)
      
      eps <- 1e-6
      mF <- mean(foc, na.rm = TRUE)
      mD <- mean(dis, na.rm = TRUE)
      effect <- (mF + eps) / (mD + eps)   # color = raw ratio F/D (no log2)
      
      pval <- NA_real_
      if (length(foc) >= 2 && length(dis) >= 2) {
        pval <- tryCatch(t.test(foc, dis)$p.value, error = function(e) NA_real_)
      }
      
      tibble(effect = effect, p = pval, nF = length(foc), nD = length(dis))
    }
  )
  
  plot_df <- bind_cols(comb, res) %>%
    mutate(sig = ifelse(!is.na(p) & p < 0.05, "*", ""))
  
  # 3) Heatmap per layer
  #    y-axis: focal cell type; x-axis: neighbor cell type; fill: F/D ratio
  p <- ggplot(plot_df, aes(x = neighbor_celltype, y = celltype.corse, fill = effect)) +
    geom_tile() +
    geom_text(aes(label = sig), vjust = 0.5, size = 3) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 1,  # midpoint=1 => no change
      name = "F/D"
    ) +
    labs(
      title = paste0(lay, " — Focal/Distal ratio"),
      x = "Neighbor cell type",
      y = "Focal cell type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank()
    )
  
  plots[[lay]] <- p
}

# Example to view:
# plots$`Layer-Exc23`
# plots$`Layer-Exc456`
# plots$`Layer-NVU`
# plots$`Layer-Olg`



plots$`Layer-Exc23`
plots$`Layer-Exc456`
plots$`Layer-NVU`
plots$`Layer-Olg`





