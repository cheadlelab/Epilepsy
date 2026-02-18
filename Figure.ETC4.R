rm(list = ls(all = TRUE))
gc()

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
# two comparisons shown left-right
cmp_vec <- c("Focal_vs_Distal", "Stimulated_vs_Distal")

cell_types <- c(
  "Exc_L23IT",
  "Exc_L234IT_ERBB4",
  "Exc_L4IT",
  "Exc_L4IT_PLCH1",
  "Exc_L5ET",
  "Exc_L5IT_GRIN3A",
  "Exc_L5IT_GRIN3A-",
  "Exc_L56IT_CAR3",
  "Exc_L56NP",
  "Exc_L6CT",
  "Exc_L6IT",
  "Exc_L6b",
  "Inh_SST",
  "Inh_PVALB",
  "Inh_VIP",
  "Inh_CXCL14",
  "Inh_LAMP5",
  "Inh_LAMP5_LHX6",
  "Chandelier"
)

# Keep your original order (you reversed before)
# cell_types <- rev(cell_types)

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

# show thin grid lines?
add_grid <- FALSE

####################################################################################################
# KEGG map00190-style nuclear-encoded OXPHOS gene sets (NO MT- genes)
####################################################################################################
# Complex I (nuclear subunits)
CI <- c(
  paste0("NDUFA", c(1,3,4,11,13)),
  paste0("NDUFB", c(2,4,8)),
  paste0("NDUFS", c(2,3,5,7))
)

# Complex II
CII <- c("SDHA", "SDHB", "SDHC", "SDHD")

# Complex III
CIII <- c(
  "UQCRB", "UQCRH", "UQCR10", "UQCR11",
  "UQCRFS1", "CYC1"
)

# Complex IV (nuclear subunits only)
CIV <- c(
  "COX4I1", "COX4I2",
  "COX5A", "COX5B",
  "COX6A1", "COX6B1",  "COX6C",
  "COX7A2", "COX7A1", "COX7B", 
  "COX7C", "COX8A"
)

# Complex V (ATP synthase, nuclear subunits)
CV <- c(
  "ATP5F1B",
  "ATP5F1D", "ATP5F1E",
  "ATP5MC1", "ATP5MC3",
  "ATP5ME", "ATP5MG",
  "ATP5PF"
)

# Your axis
PLCB_ITPR2_AXIS <- c("PLCB1", "ITPR2")

# Build a long table defining gene -> group
ETC_GROUPS <- tibble::tribble(
  ~GeneGroup,            ~Gene,
  "Complex I",           CI,
  "Complex II",          CII,
  "Complex III",         CIII,
  "Complex IV",          CIV,
  "Complex V",           CV,
  "PLCB/ITPR2 AXIS",     PLCB_ITPR2_AXIS
) %>%
  tidyr::unnest(Gene) %>%
  distinct()

####################################################################################################
# Helpers
####################################################################################################
get_deg_multicmp <- function(final_results, cell_types, cmp_vec) {
  out <- list()
  
  for (ct in cell_types) {
    if (is.null(final_results[[ct]])) next
    
    for (cmp in cmp_vec) {
      if (is.null(final_results[[ct]][[cmp]])) next
      
      deg <- final_results[[ct]][[cmp]]$deg
      if (is.null(deg) || nrow(deg) == 0) next
      
      tmp <- deg %>%
        tibble::as_tibble(rownames = "Gene") %>%
        dplyr::transmute(
          Gene = Gene,
          celltype = ct,
          Comparison = cmp,
          l2fc = avg_log2FC,
          padj = p_val_adj
        )
      
      out[[length(out) + 1]] <- tmp
    }
  }
  
  if (length(out) == 0) {
    return(tibble(
      Gene = character(), celltype = character(), Comparison = character(),
      l2fc = numeric(), padj = numeric()
    ))
  }
  
  dplyr::bind_rows(out)
}

plot_etc_twocmp_heatmap <- function(df_deg, etc_groups, cell_types, cmp_vec,
                                    expression_colors, clip_lim = 5,
                                    padj_cutoff = 0.05, l2fc_cutoff = log2(1.5),
                                    na_tile_color = "grey80",
                                    add_grid = FALSE,
                                    legend_title = "log2FC") {
  
  group_levels <- c(
    "Complex I", "Complex II", "Complex III", "Complex IV", "Complex V",
    "PLCB/ITPR2 AXIS"
  )
  
  # IMPORTANT: keep the FULL gene list to align left-right comparisons (do NOT filter by df_deg genes)
  etc_groups_keep <- etc_groups %>%
    mutate(GeneGroup = factor(GeneGroup, levels = group_levels)) %>%
    arrange(GeneGroup)
  
  # Full grid: (all ETC genes) x (all cell types) x (both comparisons)
  df_grid <- etc_groups_keep %>%
    tidyr::expand_grid(
      celltype = cell_types,
      Comparison = cmp_vec
    )
  
  # Join stats + compute significance + clip
  df_plot <- df_grid %>%
    left_join(
      df_deg %>% dplyr::select(Gene, celltype, Comparison, l2fc, padj),
      by = c("Gene", "celltype", "Comparison")
    ) %>%
    mutate(
      celltype   = factor(celltype, levels = cell_types),
      Comparison = factor(Comparison, levels = cmp_vec),
      l2fc       = ifelse(is.na(l2fc), NA_real_, pmax(pmin(l2fc, clip_lim), -clip_lim)),
      Significant = !is.na(padj) & (padj < padj_cutoff) &
        !is.na(l2fc) & (abs(l2fc) > l2fc_cutoff)
    ) %>%
    group_by(GeneGroup) %>%
    mutate(Gene = factor(Gene, levels = unique(Gene))) %>%  # keep per-group input order
    ungroup()
  
  ggplot(df_plot, aes(x = celltype, y = Gene, fill = l2fc)) +
    geom_tile(
      color = if (add_grid) "grey95" else NA,
      linewidth = if (add_grid) 0.2 else 0
    ) +
    geom_text(
      data = subset(df_plot, Significant),
      aes(label = "*"),
      size = 2.0,
      color = "black",
      na.rm = TRUE
    ) +
    scale_fill_gradientn(
      colours  = expression_colors,
      values   = scales::rescale(c(-clip_lim, -clip_lim/2, 0, clip_lim/2, clip_lim)),
      limits   = c(-clip_lim, clip_lim),
      name     = legend_title,
      na.value = na_tile_color
    ) +
    facet_grid(
      rows = vars(GeneGroup),
      cols = vars(Comparison),
      scales = "free_y",
      space  = "free_y"
    ) +
    theme_minimal(base_size = 7) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        color = "black",
        size = 6
      ),
      axis.text.y = element_text(
        face = "italic",
        color = "black",
        size = 7
      ),
      
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 7),
      
      panel.grid   = element_blank(),
      strip.text.x = element_text(size = 7, face = "bold"),
      strip.text.y = element_text(size = 7, face = "bold"),
      
      plot.margin  = margin(3, 3, 3, 3, "pt")
    )
  
}

####################################################################################################
# Run
####################################################################################################
df_deg2 <- get_deg_multicmp(final_results, cell_types, cmp_vec)

pp <- plot_etc_twocmp_heatmap(
  df_deg = df_deg2,
  etc_groups = ETC_GROUPS,
  cell_types = cell_types,
  cmp_vec = cmp_vec,
  expression_colors = expression_colors,
  clip_lim = clip_lim,
  padj_cutoff = padj_cutoff,
  l2fc_cutoff = l2fc_cutoff,
  na_tile_color = na_tile_color,
  add_grid = add_grid,
  legend_title = "log2FC"
)

####################################################################################################
# Save
####################################################################################################
out_pdf <- "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V6/Figure.Focal_vs_Distal__Stimulated_vs_Distal.ETC.side_by_side.pdf"

pdf(file = out_pdf, width = 6, height = 5.2, useDingbats = FALSE)
print(pp)
dev.off()

message("Saved: ", out_pdf)
