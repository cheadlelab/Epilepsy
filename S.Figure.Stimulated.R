rm(list = ls(all = TRUE))

## Only ggplot2 (avoid dplyr conflicts)
library(ggplot2)

## ------------------------------------
## Load data
## ------------------------------------
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
load(paste0(file_path, "/step2.1.Comparisons.sc.RData"))  # provides final_results

## ------------------------------------
## Axes and contrasts
## ------------------------------------
celltypes <- c("Exc_L6CT", "Exc_L6IT", "Exc_L6b")
contrasts <- c("Focal_vs_Distal", "Stimulated_vs_Distal")

## ------------------------------------
## Functional blocks (facet strips)
## Note: use line breaks to make multi-line facet titles
## ------------------------------------
gene_blocks <- list(
  "Immediate early\n&\nStress" = c(
    "EGR1","JUND","RHOB","ID2","FKBP5","TMSB10",
    "HSPA8","HSP90AB1","HSP90AA1","GPX4"
  ),
  "Synaptic\n&\nVesicle" = c(
    "SNAP25","SYT4","SYT11","VAMP2","RAB3A",
    "MARCKS","NPTX1","CDK5R1","CHGB","RGS4"
  ),
  "Mito\n&\nOxPhos" = c(
    "COX6A1","COX5B","COX7C","UQCR10","NDUFA4",
    "ATP5F1D","ATP6AP1","ATP6AP2","MT-ND3","CHCHD2","CHCHD10","SOD1"
  ),
  "Calcium\n&\nActivity modulators" = c(
    "CALM1","CALM2","CALM3","CAMK2N1","SLC12A5",
    "NRGN","VSNL1","CALY","HPCAL1"
  ),
  "Cytoskeleton\n&\nMicrotubule" = c(
    "TUBA1A","TUBA1B","TUBA4A","TUBB2A","TUBB3","TUBB4A",
    "STMN1","CFL1","ARPC5"
  ),
  "Endo–Lysosome\n&\nTrafficking" = c(
    "AP2M1","ATP6V1F","ATP6V1G2","ATP6V0C","ATP6V0B",
    "TRAPPC5","TMED7","TMED2","STX17","GABARAP","MAP1LC3B"
  ),
  "Proteostasis\n&\nUbiquitin" = c(
    "UBB","UBC","UBQLN2","UBL5","SUMO3","STUB1","PSMB4"
  ),
  "Adhesion\n&\nAxonGuidance" = c(
    "DCC","EPHA6","SLITRK1","SLITRK4","ITGAV","THY1","OGFRL1"
  )
)

## Flatten to a character vector of genes to plot
genes_y <- unique(unlist(gene_blocks, use.names = FALSE))

## ------------------------------------
## Style params (MATCH your earlier heatmap)
## ------------------------------------
expression_colors <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
)
clip_lim      <- 5
na_tile_color <- "grey80"
padj_cutoff   <- 0.05
l2fc_cutoff   <- log2(1.5)

## ------------------------------------
## Collect avg_log2FC AND p_val_adj per (celltype, contrast)
## (base R only)
## ------------------------------------
pull_deg <- function(ct, contrast, genes) {
  x <- final_results[[ct]][[contrast]][["deg"]]
  if (is.null(x) || nrow(x) == 0) return(data.frame())
  if (!all(c("gene","avg_log2FC","p_val_adj") %in% colnames(x))) return(data.frame())
  x <- x[x$gene %in% genes, , drop = FALSE]
  if (nrow(x) == 0) return(data.frame())
  data.frame(
    gene       = x$gene,
    avg_log2FC = x$avg_log2FC,
    p_val_adj  = x$p_val_adj,
    celltype   = ct,
    contrast   = contrast,
    stringsAsFactors = FALSE
  )
}

df_list <- list()
idx <- 1
for (ct in celltypes) {
  for (cc in contrasts) {
    df_list[[idx]] <- pull_deg(ct, cc, genes_y)
    idx <- idx + 1
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
## Map functional blocks
## ------------------------------------
block_vec <- setNames(rep(NA_character_, length(genes_y)), genes_y)
for (blk in names(gene_blocks)) {
  gset <- intersect(genes_y, gene_blocks[[blk]])
  if (length(gset) > 0) block_vec[gset] <- blk
}
gene2block <- data.frame(
  gene  = names(block_vec),
  block = as.vector(block_vec),
  stringsAsFactors = FALSE
)
df_full <- merge(df_full, gene2block, by = "gene", all.x = TRUE)

keep_blocks <- names(gene_blocks)
df_full <- df_full[!is.na(df_full$block) & df_full$block %in% keep_blocks, , drop = FALSE]

## Remove genes that are entirely NA across all tiles
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
## Contrast labels (left/right headers)
## ------------------------------------
df_full$contrast_label <- ifelse(
  df_full$contrast == "Focal_vs_Distal",
  "Focal vs Distal",
  "Stimulated vs Distal"
)
df_full$contrast_label <- factor(df_full$contrast_label,
                                 levels = c("Focal vs Distal", "Stimulated vs Distal"))

## Control x-axis ordering
df_full$celltype <- factor(df_full$celltype, levels = celltypes)

## ------------------------------------
## Clip values to match earlier heatmap style
## ------------------------------------
df_full$avg_log2FC_clip <- df_full$avg_log2FC
ok <- !is.na(df_full$avg_log2FC_clip)
df_full$avg_log2FC_clip[ok] <- pmax(pmin(df_full$avg_log2FC_clip[ok], clip_lim), -clip_lim)

## ------------------------------------
## Significance stars (MATCH earlier: padj + abs(log2FC) threshold)
## ------------------------------------
df_full$star <- ""
df_full$star[!is.na(df_full$p_val_adj) &
               df_full$p_val_adj < padj_cutoff &
               !is.na(df_full$avg_log2FC) &
               abs(df_full$avg_log2FC) > l2fc_cutoff] <- "*"

## ------------------------------------
## Plot (color style matches your earlier heatmap)
## ------------------------------------
p <- ggplot(df_full, aes(x = celltype, y = gene, fill = avg_log2FC_clip)) +
  geom_tile() +
  geom_text(aes(label = star), size = 2.2, color = "black") +
  facet_grid(block ~ contrast_label, scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = expression_colors,
    limits  = c(-clip_lim, clip_lim),
    name    = "log2FC",
    na.value = na_tile_color
  ) +
  labs(x = NULL, y = NULL, title = "L6 Neurons") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 6, face = "italic", color = "black"),
    strip.text.y      = element_text(size = 8, face = "bold"),
    strip.text.x      = element_text(size = 8, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.title      = element_text(size = 7),
    plot.margin       = margin(3, 3, 3, 3, "pt")
  )

print(p)

## ------------------------------------
## Save (4.5 x 8.5 inches)
## ------------------------------------
pdf(
  file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.Stimulated.L6_modules.pdf",
  width = 4.5, height = 8.5
)
print(p)
dev.off()
