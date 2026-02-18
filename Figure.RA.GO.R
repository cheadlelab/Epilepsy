rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
os_name <- Sys.info()[["sysname"]]   # "Linux", "Windows", "Darwin" …
file_path <- switch(
  os_name,
  "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output", 
  "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",
  stop(sprintf("Unsupported OS: %s", os_name))
)

# load(paste0(file_path, "/step1.3.annotation.AC.seurat4.RData"))
load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))

tmp <- subset(test, subset = Annotation.fine %in% c("Astrocyte"))
tmp <- JoinLayers(tmp)

# tmp[["RNA"]] <- as(object = tmp[["RNA"]], Class = "Assay5")
tmp[["RNA"]] <- split(tmp[["RNA"]], f = tmp$orig.ident)

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

tmp <- IntegrateLayers(
  object = tmp, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

ElbowPlot(tmp)

dims_parameter <- 20
tmp <- RunUMAP(tmp, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
tmp <- FindNeighbors(tmp, reduction = "harmony", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.3)
tmp@reductions$umap_harmony <- tmp@reductions$umap
tmp@meta.data$clusters_harmony <- tmp@meta.data$seurat_clusters



DefaultAssay(tmp) <- "RNA"
tmp <- JoinLayers(tmp)
Idents(tmp) <- 'seurat_clusters'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

Idents(tmp) <- "clusters_harmony"

new_identities <- c("0" = "Other Astrocyte",
                    "1" = "Other Astrocyte",
                    "2" = "Other Astrocyte",
                    "3" = "Reactive Astrocyte",
                    "4" = "Lipid-Accumulated Reactive Astrocyte",
                    "5" = "Other Astrocyte",
                    "6" = "Other Astrocyte"
)

tmp@meta.data$Annotation <- new_identities[as.character(tmp@meta.data$clusters_harmony)]
tmp$Annotation <- factor(tmp$Annotation, levels = c("Reactive Astrocyte",
                                                    "Lipid-Accumulated Reactive Astrocyte",
                                                    "Other Astrocyte"))


DimPlot(tmp, group.by = "Annotation", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE)

cell_subset <- subset(tmp, subset = Annotation == "Reactive Astrocyte")
Idents(cell_subset) <- cell_subset$condition

deg <- FindMarkers(cell_subset, 
                   ident.1 = "Focal",
                   ident.2 = "Distal",
                   logfc.threshold = log2(1.2),
                   min.pct = 0.1,
                   test.use = "MAST",
                   latent.vars = "patient")


################################################################################
# Params (align to your previous style)
library(openxlsx)

min_lfc     <- log2(1.5)
padj_cutoff <- 0.05

# Ensure gene symbol column
deg_df <- deg %>%
  as.data.frame() %>%
  tibble::rownames_to_column("symbol")

# If FindMarkers outputs "p_val_adj" as NA sometimes, keep it but avoid -Inf
deg_df <- deg_df %>%
  mutate(
    p_val_adj     = ifelse(is.na(p_val_adj), 1, p_val_adj),
    negLog10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin))
  )

sig_up <- deg_df %>%
  filter(avg_log2FC >  min_lfc, p_val_adj < padj_cutoff) %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(
    Direction   = "Upregulated",
    Rank_log2FC = row_number()
  )

sig_down <- deg_df %>%
  filter(avg_log2FC < -min_lfc, p_val_adj < padj_cutoff) %>%
  arrange(avg_log2FC) %>%
  mutate(
    Direction   = "Downregulated",
    Rank_log2FC = row_number()
  )

# Tidy column order
front_cols <- c("symbol", "avg_log2FC", "Rank_log2FC", "p_val_adj", "Direction", "negLog10_padj")
ord_cols <- function(d) c(intersect(front_cols, colnames(d)), setdiff(colnames(d), front_cols))

sig_up   <- sig_up[,   ord_cols(sig_up),   drop = FALSE]
sig_down <- sig_down[, ord_cols(sig_down), drop = FALSE]

# Output filename
contrast_name <- "RA_Focal_vs_Distal"
out_xlsx <- file.path(file_path, paste0("DEG_", contrast_name, "_UpDown.xlsx"))

wb <- createWorkbook()

addWorksheet(wb, "Upregulated")
writeData(wb, "Upregulated", sig_up)

addWorksheet(wb, "Downregulated")
writeData(wb, "Downregulated", sig_down)

# Optional: ALL sheet (not filtered by direction; already contains negLog10_padj)
addWorksheet(wb, "All_DEG")
writeData(wb, "All_DEG", deg_df %>% arrange(p_val_adj, desc(abs(avg_log2FC))))

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Saved: ", out_xlsx)
message("Up: ", nrow(sig_up), "  Down: ", nrow(sig_down))












################################################################################

celltype <- "RA"


ont_choice  <- "BP"                    # "BP", "CC", or "MF"
min_lfc     <- log2(1.5)
padj_cutoff <- 0.05


test.markers1 <- deg
test.markers1$symbol <- rownames(test.markers1)


###################################################################################################################################################
##################################################################### GO plot #####################################################################
###################################################################################################################################################

Go_up   <- list()       # list of per-celltype sig GO tables (Up)
Go_down <- list()       # list of per-celltype sig GO tables (Down)
Up_all   <- data.frame()  # combined Up across CTs
Down_all <- data.frame()  # combined Down across CTs




## stricter ORA with background/universe + size/overlap filters
min_deg_n   <- 10    # Minimum number of DEGs required to run ORA; skip ORA if fewer than this
min_overlap <- 3     # Minimum overlap (Count) of DEGs per GO term (require Count >= min_overlap)
minGSSize   <- 5    # Minimum GO term size (number of genes in the term) to keep
maxGSSize   <- 2000  # Maximum GO term size to keep (filter out very broad terms)

ego_up <- NULL
ego_down <- NULL


bg_universe <- unique(test.markers1$symbol[!is.na(test.markers1$symbol)])

genes.up <- unique(test.markers1$symbol[
  test.markers1$avg_log2FC >  min_lfc & test.markers1$p_val_adj < padj_cutoff
])
genes.down <- unique(test.markers1$symbol[
  test.markers1$avg_log2FC < -min_lfc & test.markers1$p_val_adj < padj_cutoff
])

# genes.up <- genes.up[!grepl("^HSP", genes.up)]
# genes.down <- genes.down[!grepl("^HSP", genes.down)]

if (length(genes.up) >= min_deg_n) {
  ego_up <- tryCatch(
    enrichGO(
      gene          = genes.up,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont_choice,     # BP/CC/MF
      universe      = bg_universe,    
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    ),
    error = function(e) NULL
  )
}

if (length(genes.down) >= min_deg_n) {
  ego_down <- tryCatch(
    enrichGO(
      gene          = genes.down,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont_choice,
      universe      = bg_universe,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    ),
    error = function(e) NULL
  )
}

if (!is.null(ego_up) && nrow(ego_up@result) > 0) {
  ego_up@result$Direction <- "Upregulated"
  ego_up_sub <- subset(ego_up@result, p.adjust < padj_cutoff & Count >= min_overlap)
  if (nrow(ego_up_sub) > 0) {
    Go_up[[celltype]] <- ego_up_sub
  } else {
    warning(sprintf("[%s] Upregulated GO: filtered to empty (p.adjust/Count).", celltype))
  }
} else {
  warning(sprintf("[%s] Upregulated GO results are missing or empty.", celltype))
}

if (!is.null(ego_down) && nrow(ego_down@result) > 0) {
  ego_down@result$Direction <- "Downregulated"
  ego_down_sub <- subset(ego_down@result, p.adjust < padj_cutoff & Count >= min_overlap)
  if (nrow(ego_down_sub) > 0) {
    Go_down[[celltype]] <- ego_down_sub
  } else {
    warning(sprintf("[%s] Downregulated GO: filtered to empty (p.adjust/Count).", celltype))
  }
} else {
  warning(sprintf("[%s] Downregulated GO results are missing or empty.", celltype))
}







combined_go <- data.frame()

if (exists("ego_up") && 
    !is.null(ego_up) && 
    nrow(ego_up@result) > 0) {
  ego_up_sub <- ego_up@result[ego_up@result$p.adjust < 0.05, ]
} else {
  ego_up_sub <- data.frame()
  warning("Upregulated GO results are missing or empty.")
}

if (exists("ego_down") && 
    !is.null(ego_down) && 
    nrow(ego_down@result) > 0) {
  ego_down_sub <- ego_down@result[ego_down@result$p.adjust < 0.05, ]
} else {
  ego_down_sub <- data.frame()
  warning("Downregulated GO results are missing or empty.")
}

combined_go <- rbind(ego_up_sub, ego_down_sub) %>%
  arrange(p.adjust) %>%
  group_by(ID) %>%
  ungroup() %>%
  head(10)


###################################################################################################################################################
############################################################### GO Result Visualization ##########################################################
###################################################################################################################################################
if (nrow(combined_go) > 0) {
  combined_go$logp <- -log10(combined_go$p.adjust)
  
  p1 <- ggplot(combined_go, 
               aes(x = logp, 
                   y = reorder(Description, logp), 
                   fill = Direction)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8")) +
    labs(
      title = celltype,
      x = "-log10(p.adjust)",
      y = "GO Term"
    ) +
    theme_classic(base_size = 6) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
      axis.text.x = element_text(size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      legend.position = "none",
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6)
    )
}

pdf(file = paste0("C:/Users/thech/Desktop/tmp.pdf"), width = 2.2, height = 1.5)
p1
dev.off() # Close the PDF device
