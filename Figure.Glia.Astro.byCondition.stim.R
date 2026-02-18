rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(openxlsx)
})

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


DefaultAssay(tmp) <- "RNA"
tmp <- JoinLayers(tmp)
Idents(tmp) <- 'Annotation'
markers <- FindAllMarkers(object = tmp, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0)
markers$symbol <- markers$gene


########################################################################################################################################################
##################################################################### Volcano plot #####################################################################
########################################################################################################################################################

celltype <- "Lipid-Accumulated Reactive Astrocyte"
celltype <- "Reactive Astrocyte"
tt.tmp <- subset(tmp, subset = Annotation == celltype)

Idents(tt.tmp) <- 'condition'
test.markers1 <- FindMarkers(tt.tmp, 
                             ident.1 = "Stimulated",
                             ident.2 = "Distal",
                             min.pct = 0.1,
                             test.use = "MAST",
                             latent.vars = "patient")

test.markers1$symbol <- rownames(test.markers1)
test.markers1$gene <- rownames(test.markers1)

fc_thr   <- log2(1.5)
eps      <- 1e-300

pattern_rm <- "^(AC|AL|AP)[0-9]+(?:\\.[0-9]+)?$|^LINC|.-AS[0-9]+$"

test.markers1 <- test.markers1 %>%
  filter(!grepl(pattern_rm, gene, ignore.case = TRUE)) %>%
  filter(!grepl("^MT-", gene, ignore.case = TRUE))

df <- test.markers1 %>%
  mutate(log10fdr = -log10(p_val_adj + eps)) %>%
  filter(!grepl(pattern_rm, gene, ignore.case = TRUE)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(sig = p_val_adj <= 0.05)

# force_gene <- c("APOE", "FTL", "AQP4", "CST3")
force_gene <- NULL

top_up   <- df %>% filter(sig, avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 10)
top_down <- df %>% filter(sig, avg_log2FC < 0) %>% arrange(avg_log2FC)       %>% slice_head(n = 10)

force_row <- df %>%
  filter(sig) %>%
  filter(toupper(gene) %in% toupper(force_gene))

top_lab <- bind_rows(top_up, top_down, force_row) %>% distinct(gene, .keep_all = TRUE)

mt_in_top <- top_lab$gene[grepl("^MT-", top_lab$gene)]

df$category <- "Other"
df$category[df$sig & df$gene %in% c(top_up$gene, force_gene)] <- "Top_Up"
df$category[df$sig & df$gene %in% top_down$gene]              <- "Top_Down"
df$category[df$sig & df$gene %in% mt_in_top]                  <- "MT_Top"

cols <- c(Other = "grey75", Top_Up = "#d62728", Top_Down = "#1f77b4", MT_Top = "#ff7f0e")

label_data <- top_lab %>%
  filter(p_val_adj <= 0.05) %>%
  filter(!grepl("^MT-", gene))


p <- ggplot(df, aes(x = avg_log2FC, y = log10fdr)) +
  geom_point(aes(color = category), size = 0.8, alpha = 0.85) +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = c(-fc_thr, fc_thr),
             linetype = "dashed",
             linewidth = 0.3,
             color = "grey50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.3,
             color = "grey50") +
  labs(title = "RA",
       x = "avg_log2FC",
       y = expression(-log[10]("adj. p"))) +
  theme_bw(base_size = 6) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold"))

pt_to_gg <- function(pt) pt / 2.845

p1 <- p + ggrepel::geom_text_repel(
  data = label_data,
  aes(label = gene),
  size = pt_to_gg(6),
  seed = 233,
  box.padding = 0.2,
  point.padding = 0.1,
  min.segment.length = 0,
  segment.size = 0.2,
  max.overlaps = Inf
)



pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.RA.Vol.stim.pdf", width = 2, height = 2)
print(p1)
dev.off() # Close the PDF device

###################################################################################################################################################

Idents(tmp) <- "clusters_harmony"

new_identities <- c("0" = "Other",
                    "1" = "Other",
                    "2" = "Other",
                    "3" = "RA",
                    "4" = "LARA",
                    "5" = "Other",
                    "6" = "Other"
)

tmp@meta.data$new <- new_identities[as.character(tmp@meta.data$clusters_harmony)]
tmp$new <- factor(tmp$new, levels = c("RA",
                                      "LARA",
                                      "Other"))

pseudo_tmp <- AggregateExpression(
  tmp,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("new")
)



pattern_rm <- "^(AC|AL|AP)[0-9]+(?:\\.[0-9]+)?$|^LINC|.-AS[0-9]+$|^MT-|^RPS|^RPL"
markers_use <- markers %>% dplyr::filter(!grepl(pattern_rm, gene, ignore.case = TRUE))

top_by_cluster <- markers_use %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

clust_levels <- levels(Idents(pseudo_tmp))
top_by_cluster$cluster <- factor(top_by_cluster$cluster, levels = clust_levels)

features_heat <- unique(top_by_cluster$gene)
features_heat <- features_heat[features_heat %in% rownames(pseudo_tmp)]

suppressMessages({
  pseudo_tmp <- ScaleData(pseudo_tmp, features = features_heat, verbose = FALSE)
})

p_heat <- DoHeatmap(
  object   = pseudo_tmp,
  features = features_heat,
  slot     = "scale.data",
  group.by = NULL,
  draw.lines = FALSE,
  raster   = TRUE,
  size     = 3,
  angle    = 0
  
) + ggplot2::theme(
  axis.text.y   = ggplot2::element_text(size = 6),
  axis.text.x   = ggplot2::element_text(size = 6, angle = 0),
  legend.title  = ggplot2::element_text(size = 6),
  legend.text   = ggplot2::element_text(size = 6)
)

rb_cols <- c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDD0A2","#EF8A62","#B2182B")
p_heat <- p_heat + scale_fill_gradientn(colors = rb_cols, limits = c(-1, 1), oob = scales::squish)

print(p_heat)

# 1.5 * 2.5

pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V3/Figure.LARA.Heatmap.pdf", width = 2.5, height = 2.8)
print(p_heat)
dev.off() # Close the PDF device


###################################################################################################################################################
##################################################################### GO plot #####################################################################
###################################################################################################################################################
celltype <- "Reactive Astrocyte" # "Lipid-Accumulated Reactive Astrocyte" "Reactive Astrocyte"
test.markers1 <- test.markers1

ont_choice  <- "BP"                    # "BP", "CC", or "MF"
min_lfc     <- log2(1.5)
padj_cutoff <- 0.05


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


pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.RA.GO.stim.pdf", width = 4, height = 1.5)
print(p1)
dev.off() # Close the PDF device


pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.LARA.GO.pdf", width = 3, height = 2.5)
print(p1)
dev.off() # Close the PDF device
