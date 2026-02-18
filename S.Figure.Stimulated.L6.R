rm(list = ls(all = TRUE))

## Only ggplot2 (avoid dplyr conflicts)
library(ggplot2)

## ------------------------------------
## Load data
## ------------------------------------
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
load(paste0(file_path, "/step2.1.Comparisons.sc.RData"))

## Axes and contrasts
celltypes  <- c("Exc_L6CT", "Exc_L6IT", "Exc_L6b")
contrasts  <- c("Focal_vs_Distal", "Stimulated_vs_Distal")  # for completeness

## ------------------------------------
## Functional blocks (facet strips)
## Note: use line breaks to make multi-line facet titles
## ------------------------------------
gene_blocks <- list(
  # 1) Immediate early & Stress
  "Immediate early\n&\nStress" = c(
    "EGR1","JUND","RHOB","ID2","FKBP5","TMSB10",
    "HSPA8","HSP90AB1","HSP90AA1","GPX4"
  ),
  
  # 2) Synaptic & Vesicle
  "Synaptic\n&\nVesicle" = c(
    "SNAP25","SYT4","SYT11","VAMP2","RAB3A",
    "MARCKS","NPTX1","CDK5R1","CHGB","RGS4"
  ),
  
  # 3) Mito & OxPhos
  "Mito\n&\nOxPhos" = c(
    "COX6A1","COX5B","COX7C","UQCR10","NDUFA4",
    "ATP5F1D","ATP6AP1","ATP6AP2","MT-ND3","CHCHD2","CHCHD10","SOD1"
  ),
  
  # 4) Calcium & Activity modulators
  "Calcium\n&\nActivity modulators" = c(
    "CALM1","CALM2","CALM3","CAMK2N1","SLC12A5",
    "NRGN","VSNL1","CALY","HPCAL1"
  ),
  
  # 5) Cytoskeleton & Microtubule
  "Cytoskeleton\n&\nMicrotubule" = c(
    "TUBA1A","TUBA1B","TUBA4A","TUBB2A","TUBB3","TUBB4A",
    "STMN1","CFL1","ARPC5"
  ),
  
  # 6) Endo–Lysosome & Trafficking
  "Endo–Lysosome\n&\nTrafficking" = c(
    "AP2M1","ATP6V1F","ATP6V1G2","ATP6V0C","ATP6V0B",
    "TRAPPC5","TMED7","TMED2","STX17","GABARAP","MAP1LC3B"
  ),
  
  # 7) Proteostasis & Ubiquitin
  "Proteostasis\n&\nUbiquitin" = c(
    "UBB","UBC","UBQLN2","UBL5","SUMO3","STUB1","PSMB4"
  ),
  
  # 8) Axon guidance & Adhesion
  "Adhesion\n&\nAxonGuidance" = c(
    "DCC","EPHA6","SLITRK1","SLITRK4","ITGAV","THY1","OGFRL1"
  )
)

## Flatten to a character vector of genes to plot
genes_y <- unique(unlist(gene_blocks, use.names = FALSE))

## ------------------------------------
## Collect avg_log2FC AND p_val_adj per (celltype, contrast)
## ------------------------------------
pull_deg <- function(ct, contrast, genes) {
  x <- final_results[[ct]][[contrast]][["deg"]]
  if (is.null(x) || nrow(x) == 0) return(data.frame())
  x <- x[x$gene %in% genes, , drop = FALSE]
  need_cols <- c("gene","avg_log2FC","p_val_adj")
  if (!all(need_cols %in% colnames(x))) return(data.frame())
  data.frame(
    gene      = x$gene,
    avg_log2FC= x$avg_log2FC,
    p_val_adj = x$p_val_adj,
    celltype  = ct,
    contrast  = contrast,
    stringsAsFactors = FALSE
  )
}

## Build long table for 3 cell types × 2 contrasts
df_list <- list()
for (ct in celltypes) {
  for (cc in contrasts) {
    key <- paste(ct, cc, sep = "_")
    df_list[[key]] <- pull_deg(ct, cc, genes_y)
  }
}
df <- if (length(df_list)) do.call(rbind, df_list) else data.frame()

## ------------------------------------
## Complete grid (gene × celltype × contrast); merge stats
## ------------------------------------
df_full <- expand.grid(
  gene     = genes_y,
  celltype = celltypes,
  contrast = contrasts,
  stringsAsFactors = FALSE
)
if (nrow(df) > 0) {
  df_full <- merge(df_full, df, by = c("gene","celltype","contrast"), all.x = TRUE)
} else {
  df_full$avg_log2FC <- NA_real_
  df_full$p_val_adj  <- NA_real_
}

## ------------------------------------
## Map functional blocks and keep only defined blocks
## ------------------------------------
block_vec <- setNames(rep(NA_character_, length(genes_y)), genes_y)
for (blk in names(gene_blocks)) {
  gset <- intersect(genes_y, gene_blocks[[blk]])
  if (length(gset) > 0) block_vec[gset] <- blk
}
gene2block <- data.frame(gene = names(block_vec), block = as.vector(block_vec), stringsAsFactors = FALSE)
df_full <- merge(df_full, gene2block, by = "gene", all.x = TRUE)

keep_blocks <- names(gene_blocks)
df_full <- df_full[!is.na(df_full$block) & df_full$block %in% keep_blocks, , drop = FALSE]

## Remove genes that are entirely NA across all tiles (cleaner plot)
has_value <- tapply(df_full$avg_log2FC, df_full$gene, function(v) any(!is.na(v)))
df_full <- df_full[df_full$gene %in% names(has_value)[has_value], , drop = FALSE]

## ------------------------------------
## Ordering: blocks and genes within blocks
## ------------------------------------
df_full$block <- factor(df_full$block, levels = keep_blocks)

ordered_genes <- character(0)
for (blk in keep_blocks) {
  gs <- intersect(gene_blocks[[blk]], unique(df_full$gene[df_full$block == blk]))
  ordered_genes <- c(ordered_genes, gs)
}
df_full$gene <- factor(df_full$gene, levels = ordered_genes)

## ------------------------------------
## Facet columns: top headers for contrasts (left = Focal, right = Stimulated)
## x-axis will show only the cell type (rotated 45°)
## ------------------------------------
df_full$contrast_label <- ifelse(
  df_full$contrast == "Focal_vs_Distal",
  "Focal vs Distal",
  "Stimulated vs Distal"
)
df_full$contrast_label <- factor(df_full$contrast_label,
                                 levels = c("Focal vs Distal", "Stimulated vs Distal")
)

## Control the x-axis (cell types) ordering
df_full$celltype <- factor(df_full$celltype, levels = c("Exc_L6CT","Exc_L6IT","Exc_L6b"))

## ------------------------------------
## Significance stars per tile (p_val_adj < 0.05)
## ------------------------------------
df_full$star <- ""
df_full$star[!is.na(df_full$p_val_adj) & df_full$p_val_adj < 0.05] <- "*"

## ------------------------------------
## Plot
## ------------------------------------
p <- ggplot(df_full, aes(x = celltype, y = gene, fill = avg_log2FC)) +
  geom_tile() +
  geom_text(aes(label = star), size = 2.3, vjust = 0.5, hjust = 0.5) +
  facet_grid(block ~ contrast_label, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    name = "avg_log2FC",
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, na.value = "grey95"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "L6 Neurons"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),  # rotate 45°
    axis.text.y     = element_text(size = 6),
    strip.text.y    = element_text(size = 8, face = "bold"),
    strip.text.x    = element_text(size = 8, face = "bold"),
    legend.key.height = grid::unit(4, "mm")
  )

print(p)

# 4.5*8.5