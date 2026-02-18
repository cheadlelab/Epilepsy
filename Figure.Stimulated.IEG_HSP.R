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
load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))  # provides final_results

####################################################################################################
# Params
####################################################################################################
# Immediate Early Genes (IEG)
IEG <- c("FOS", "JUN", "EGR1", "EGR2", "ARC", "NR4A1", "NPAS4", "DUSP6")

# Heat shock / stress genes (HSP)
HSP <- c("HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPH1")

cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1",
  "Exc_L5IT_GRIN3A", "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L6IT",
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier",
  "Astrocyte", "Microglia", "OPC", "Olg"
)
cell_types <- rev(cell_types)

cmp_keep <- "Stimulated_vs_Distal"

# heatmap colors
expression_colors <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
)

clip_lim <- 5

# significance thresholds
padj_cutoff <- 0.05
l2fc_cutoff <- log2(1.5)

# missing tile color
na_tile_color <- "grey80"

####################################################################################################
# Helpers
####################################################################################################
# Extract avg_log2FC + adjusted p from final_results for ONE comparison into long df
get_deg_onecmp <- function(final_results, cell_types, cmp_keep) {
  out <- list()
  
  for (ct in cell_types) {
    if (is.null(final_results[[ct]])) next
    if (is.null(final_results[[ct]][[cmp_keep]])) next
    
    deg <- final_results[[ct]][[cmp_keep]]$deg
    if (is.null(deg) || nrow(deg) == 0) next
    
    tmp <- deg %>%
      tibble::as_tibble(rownames = "Gene") %>%
      dplyr::transmute(
        Gene = Gene,
        celltype = ct,
        l2fc = avg_log2FC,
        padj = p_val_adj
      )
    
    out[[length(out) + 1]] <- tmp
  }
  
  if (length(out) == 0) {
    return(tibble(Gene=character(), celltype=character(), l2fc=numeric(), padj=numeric()))
  }
  dplyr::bind_rows(out)
}

# Clean facet heatmap: each panel has its OWN gene list; missing Gene×celltype tiles are grey
plot_stim_ieg_hsp_heatmap_clean <- function(df_deg, IEG, HSP, cell_types,
                                            expression_colors, clip_lim = 5,
                                            padj_cutoff = 0.05, l2fc_cutoff = log2(1.5),
                                            na_tile_color = "grey80",
                                            add_grid = FALSE) {
  
  genes_present <- unique(df_deg$Gene)
  
  IEG_keep <- IEG[IEG %in% genes_present]
  HSP_keep <- HSP[HSP %in% genes_present]
  
  # Build per-panel grids (facets contain ONLY their own genes)
  df_grid_ieg <- tidyr::expand_grid(GeneGroup = "IEG", Gene = IEG_keep, celltype = cell_types)
  df_grid_hsp <- tidyr::expand_grid(GeneGroup = "HSP", Gene = HSP_keep, celltype = cell_types)
  
  # Join separately to preserve factor levels (NO ifelse on factors!)
  df_plot_ieg <- df_grid_ieg %>%
    dplyr::left_join(
      df_deg %>% dplyr::select(Gene, celltype, l2fc, padj),
      by = c("Gene", "celltype")
    ) %>%
    dplyr::mutate(Gene = factor(Gene, levels = IEG_keep))
  
  df_plot_hsp <- df_grid_hsp %>%
    dplyr::left_join(
      df_deg %>% dplyr::select(Gene, celltype, l2fc, padj),
      by = c("Gene", "celltype")
    ) %>%
    dplyr::mutate(Gene = factor(Gene, levels = HSP_keep))
  
  df_plot <- dplyr::bind_rows(df_plot_ieg, df_plot_hsp) %>%
    dplyr::mutate(
      # clip (keep NA as NA)
      l2fc = ifelse(is.na(l2fc), NA_real_, pmax(pmin(l2fc, clip_lim), -clip_lim)),
      Significant = !is.na(padj) & (padj < padj_cutoff) & !is.na(l2fc) & (abs(l2fc) > l2fc_cutoff),
      celltype = factor(celltype, levels = cell_types),
      GeneGroup = factor(GeneGroup, levels = c("IEG", "HSP"))
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
      name    = "log2FC (Stimulated vs Distal)",
      na.value = na_tile_color
    ) +
    facet_grid(cols = vars(GeneGroup), scales = "free_x", space = "free_x") +
    scale_x_discrete(position = "bottom", drop = TRUE) +
    theme_minimal(base_size = 7) +
    theme(
      axis.title.x   = element_blank(),
      axis.title.y   = element_blank(),
      axis.text.x    = element_text(angle = 45, hjust = 1, face = "italic", color = "black", size = 7),
      axis.text.y    = element_text(color = "black", size = 6),
      legend.text    = element_text(size = 7),
      legend.title   = element_text(size = 7),
      panel.grid     = element_blank(),
      strip.text.x   = element_text(size = 7, face = "bold"),
      plot.margin    = margin(3, 3, 3, 3, "pt")
      # NOTE: I did NOT change your facet spacing; no panel.spacing.x here
    )
}

####################################################################################################
# Run
####################################################################################################
df_deg <- get_deg_onecmp(final_results, cell_types, cmp_keep)

pp <- plot_stim_ieg_hsp_heatmap_clean(
  df_deg = df_deg,
  IEG = IEG,
  HSP = HSP,
  cell_types = cell_types,
  expression_colors = expression_colors,
  clip_lim = clip_lim,
  padj_cutoff = padj_cutoff,
  l2fc_cutoff = l2fc_cutoff,
  na_tile_color = na_tile_color,
  add_grid = FALSE
)

####################################################################################################
# Save
####################################################################################################
pdf(
  file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.Stimulated.IEG_HSP.pdf",
  width = 4, height = 2.8
)
print(pp)
dev.off()
