rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(tibble)
})

####################################################################################################
# Load
####################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))  # should provide final_results

####################################################################################################
# Params
####################################################################################################
IEG <- c("FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "ARC", "NR4A1", "NR4A2", "NPAS4", "ATF3", "DUSP1", "DUSP6")
genes_of_interest <- IEG

cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1",
  "Exc_L5IT_GRIN3A", "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L6IT",
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier",
  "Astrocyte", "Microglia", "OPC", "Olg"
)
cell_types <- rev(cell_types)

comparisons <- c("Focal_vs_Distal", "Stimulated_vs_Distal")

# heatmap colors
expression_colors <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
)

# clip range
clip_lim <- 5

# NA tile color (missing Gene×celltype -> filled as NA after complete())
na_tile_color <- "grey80"

####################################################################################################
# Helpers
####################################################################################################
# 1) Extract log2FC from final_results$deg (avg_log2FC) into a long df
get_deg_l2fc <- function(final_results, cell_types, comparisons, genes_keep) {
  out <- list()
  
  for (ct in cell_types) {
    if (is.null(final_results[[ct]])) next
    
    for (cmp in comparisons) {
      deg <- final_results[[ct]][[cmp]]$deg
      if (is.null(deg) || nrow(deg) == 0) next
      
      deg <- deg %>%
        tibble::as_tibble(rownames = "Gene") %>%
        filter(Gene %in% genes_keep) %>%
        transmute(
          Gene = Gene,
          celltype = ct,
          comparison = cmp,
          l2fc = avg_log2FC
        )
      
      out[[length(out) + 1]] <- deg
    }
  }
  
  if (length(out) == 0) {
    return(tibble(Gene = character(), celltype = character(), comparison = character(), l2fc = numeric()))
  }
  bind_rows(out)
}

# 2) Significance mask (same thresholds as you used)
get_sig_hits <- function(final_results, cell_types, comparisons, genes_keep,
                         padj_cutoff = 0.05, l2fc_cutoff = log2(1.5)) {
  sig_hits <- list()
  
  for (ct in cell_types) {
    if (is.null(final_results[[ct]])) next
    
    for (cmp in comparisons) {
      deg <- final_results[[ct]][[cmp]]$deg
      if (is.null(deg) || nrow(deg) == 0) next
      
      tmp <- deg %>%
        tibble::as_tibble(rownames = "Gene") %>%
        filter(
          Gene %in% genes_keep,
          p_val_adj < padj_cutoff,
          abs(avg_log2FC) > l2fc_cutoff
        ) %>%
        transmute(Gene = Gene, celltype = ct, comparison = cmp)
      
      sig_hits[[length(sig_hits) + 1]] <- tmp
    }
  }
  
  if (length(sig_hits) == 0) {
    return(tibble(Gene = character(), celltype = character(), comparison = character()))
  }
  bind_rows(sig_hits) %>% distinct()
}

# 3) Build focal-only heatmap df + plot (Focal_vs_Distal) with completed tiles and NA as grey
plot_focal_ieg_heatmap <- function(df_l2fc, sig_hits, genes_keep, cell_types,
                                   expression_colors, clip_lim = 2.5,
                                   na_tile_color = "grey80",
                                   add_grid = FALSE) {
  cmp_keep <- "Focal_vs_Distal"
  
  df_plot <- df_l2fc %>%
    filter(comparison == cmp_keep) %>%
    tidyr::complete(
      Gene = genes_keep,
      celltype = cell_types,
      fill = list(l2fc = NA_real_)
    ) %>%
    mutate(
      # clip (keep NA as NA)
      l2fc = ifelse(is.na(l2fc), NA_real_, pmax(pmin(l2fc, clip_lim), -clip_lim)),
      Gene = factor(Gene, levels = genes_keep),
      celltype = factor(celltype, levels = cell_types),
      comparison = cmp_keep,
      Significant = paste(Gene, celltype, comparison) %in%
        paste(sig_hits$Gene, sig_hits$celltype, sig_hits$comparison)
    )
  
  ggplot(df_plot, aes(x = Gene, y = celltype, fill = l2fc)) +
    geom_tile(
      color = if (add_grid) "grey95" else NA,
      linewidth = if (add_grid) 0.2 else 0
    ) +
    geom_text(
      data = subset(df_plot, Significant),
      aes(label = "*"),
      size = 2.2, color = "black", na.rm = TRUE
    ) +
    scale_fill_gradientn(
      colours = expression_colors,
      values  = scales::rescale(c(-2, -1, 0, 1, 2)),
      limits  = c(-clip_lim, clip_lim),
      name    = "log2FC (Focal vs Distal)",
      na.value = na_tile_color
    ) +
    scale_x_discrete(position = "bottom") +
    theme_minimal(base_size = 7) +
    theme(
      axis.title.x   = element_blank(),
      axis.title.y   = element_blank(),
      axis.text.x    = element_text(angle = 45, hjust = 1, face = "italic", color = "black", size = 7),
      axis.text.y    = element_text(color = "black", size = 6),
      legend.text    = element_text(size = 7),
      legend.title   = element_text(size = 7),
      panel.grid     = element_blank(),
      plot.margin    = margin(3, 3, 3, 3, "pt")
    )
}

####################################################################################################
# Run
####################################################################################################
genes_keep <- unique(intersect(genes_of_interest, unique(unlist(lapply(cell_types, function(ct) {
  # best-effort: find genes from any deg table (doesn't break if missing)
  g <- c()
  for (cmp in comparisons) {
    deg <- final_results[[ct]][[cmp]]$deg
    if (!is.null(deg) && nrow(deg) > 0) g <- c(g, rownames(deg))
  }
  g
})))))

# keep order as your IEG list
genes_keep <- genes_of_interest[genes_of_interest %in% genes_keep]

df_l2fc  <- get_deg_l2fc(final_results, cell_types, comparisons, genes_keep)
sig_hits <- get_sig_hits(final_results, cell_types, comparisons, genes_keep)

pp <- plot_focal_ieg_heatmap(
  df_l2fc = df_l2fc,
  sig_hits = sig_hits,
  genes_keep = genes_keep,
  cell_types = cell_types,
  expression_colors = expression_colors,
  clip_lim = clip_lim,
  na_tile_color = na_tile_color,
  add_grid = FALSE  # set TRUE if you want faint grid lines
)

####################################################################################################
# Save
####################################################################################################
pdf("//grid/cheadle_home/qianyu/Epilepsy/manuscript/V3/Figure.Focal.IEG.pdf", width = 3.5, height = 2.8)
print(pp)
dev.off()
