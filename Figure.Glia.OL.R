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

tmp <- subset(test, subset = Annotation.fine %in% c("Olg"))
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

p2 <- DimPlot(tmp, reduction = "umap_harmony", group.by = "clusters_harmony", raster=FALSE, label = TRUE)
p4 <- DimPlot(tmp, reduction = "umap_harmony", group.by = "condition", raster=FALSE, label = TRUE)
p2+p4

Idents(tmp) <- "clusters_harmony"

new_identities <- c("0" = "OL",
                    "1" = "OL",
                    "2" = "OxPos OL",
                    "3" = "OL",
                    "4" = "OL",
                    "5" = "OL",
                    "6" = "OL"
)

tmp@meta.data$Annotation <- new_identities[as.character(tmp@meta.data$clusters_harmony)]

tmp@meta.data$Annotation <- new_identities[as.character(tmp@meta.data$clusters_harmony)]
tmp$Annotation <- factor(tmp$Annotation, levels = c("OL",
                                                    "OxPos OL"))




DefaultAssay(tmp) <- "RNA"
tmp <- JoinLayers(tmp)
Idents(tmp) <- 'Annotation'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DimPlot(tmp, reduction = "umap_harmony", group.by = "Annotation", raster=FALSE, label = TRUE)

##################################################################################################################################
## addmodulescore
Markers.MT <- grep("^MT-", rownames(tmp), value = TRUE)
Markers.RP <- grep("^RP[SL]", rownames(tmp), value = TRUE)
Markers.CI <- c("NDUFA1", "NDUFA4", "NDUFA6", "NDUFB2", "NDUFB5", "NDUFS3", "NDUFS7")
Markers.CIV <- c("COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2", "COX7C", "COX8A")
Markers.CV <- c("ATP5F1A", "ATP5F1B", "ATP5MC1", "ATP5MC2", "ATP5PF")
Markers.CIs <- c("NDUFA1", "NDUFA4", "NDUFA6", "NDUFB2", "NDUFB5", "NDUFS3", "NDUFS7",
                 "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2", "COX7C", "COX8A",
                 "ATP5F1A", "ATP5F1B", "ATP5MC1", "ATP5MC2", "ATP5PF")

tmp <- AddModuleScore(tmp, list(Markers.MT), name="MT_Score")
tmp <- AddModuleScore(tmp, list(Markers.RP), name="RP_Score")
tmp <- AddModuleScore(tmp, features = list(Markers.CI),  name = "CI_Score")
tmp <- AddModuleScore(tmp, features = list(Markers.CIV), name = "CIV_Score")
tmp <- AddModuleScore(tmp, features = list(Markers.CV),  name = "CV_Score")
tmp <- AddModuleScore(tmp, features = list(Markers.CIs),  name = "CIs_Score")

minmax <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
tmp$MT_Score <- minmax(tmp$MT_Score1)
tmp$RP_Score <- minmax(tmp$RP_Score1)
tmp$CI_Score  <- minmax(tmp$CI_Score1)
tmp$CIV_Score <- minmax(tmp$CIV_Score1)
tmp$CV_Score  <- minmax(tmp$CV_Score1)
tmp$CIs_Score  <- minmax(tmp$CIs_Score1)

FeaturePlot(tmp, features = c("MT_Score"))
FeaturePlot(tmp, features = c("RP_Score"))
FeaturePlot(tmp, features = c("CI_Score"))
FeaturePlot(tmp, features = c("CIV_Score"))
FeaturePlot(tmp, features = c("CV_Score"))
FeaturePlot(tmp, features = c("CIs_Score"))

oxpos_minimal <- c(
  "COX4I1",
  "COX5B",
  "NDUFA4",
  "ATP5F1B"
)

FeaturePlot(tmp, features = oxpos_minimal,order = TRUE)
##################################################################################################################################
# vlnplot for QC again
# 10 3.5
VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal", "log10GenesPerUMI"), ncol = 5, group.by = "Annotation", pt.size = 0)

##################################################################################################################################
# umap for OLs
save_path <- "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5"

pdf(file = paste0(save_path, "/Figure.Glia.nolegend.pdf"), width = 5, height = 5)
DimPlot(tmp, group.by = "Annotation", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.Glia.nolegend.byCondition.pdf"), width = 15, height = 5)
DimPlot(tmp, group.by = "Annotation", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE, split.by = "condition") + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device



pdf(file = paste0(save_path, "/Figure.Glia.legend.pdf"), width = 5, height = 5)
DimPlot(tmp, group.by = "Annotation", raster = FALSE)
dev.off() # Close the PDF device


pdf(file = paste0(save_path, "/Figure.OL.MT.pdf"), width = 5, height = 5)
FeaturePlot(
  tmp, features = c("MT_Score"),
  min.cutoff = c("0.3"), 
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.5, order = TRUE,
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.OL.RP.pdf"), width = 5, height = 5)
FeaturePlot(
  tmp, features = c("RP_Score"),
  min.cutoff = c("0.3"), 
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.5, order = TRUE,
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device


pdf(file = paste0(save_path, "/Figure.OL.CIs.pdf"), width = 10, height = 10)
FeaturePlot(
  tmp, features = c("CIs_Score"),
  min.cutoff = c("0.3"), 
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.1
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device





oxpos_minimal <- c(
  "COX4I1",
  "COX5B",
  "NDUFA4",
  "ATP5F1B"
)

pdf(file = paste0(save_path, "/Figure.OL.COX4I1.pdf"), width = 5, height = 5)
FeaturePlot(
  tmp, features = c("COX4I1"),
  min.cutoff = c("0.3"), 
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.5, order = TRUE,
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.OL.COX5B.pdf"), width = 5, height = 5)
FeaturePlot(
  tmp, features = c("COX5B"),
  min.cutoff = c("0.3"), 
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.5, order = TRUE,
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.OL.NDUFA4.pdf"), width = 5, height = 5)
FeaturePlot(
  tmp, features = c("NDUFA4"),
  min.cutoff = c("0.3"), 
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.5, order = TRUE,
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.OL.ATP5F1B.pdf"), width = 5, height = 5)
FeaturePlot(
  tmp, features = c("ATP5F1B"),
  min.cutoff = c("0.3"), 
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.5, order = TRUE,
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device


#####################################################################################################################################################################
##################################################################### Volcano plot by condition #####################################################################
#####################################################################################################################################################################

celltype <- "OxPos OL"
tt.tmp <- subset(tmp, subset = Annotation == celltype)

Idents(tt.tmp) <- 'condition'
test.markers1 <- FindMarkers(tt.tmp, 
                             ident.1 = "Focal",
                             ident.2 = "Distal",
                             min.pct = 0.1,
                             test.use = "MAST",
                             latent.vars = "patient")

test.markers1$symbol <- rownames(test.markers1)
test.markers1$gene <- rownames(test.markers1)

test.markers1 <- test.markers1 %>%
  filter(!grepl("^MT-|^mt-", gene))

fc_thr   <- log2(1.5)
eps      <- 1e-300

pattern_rm <- "^(AC|AL|AP)[0-9]+(?:\\.[0-9]+)?$|^LINC|.-AS[0-9]+$"

df <- test.markers1 %>%
  mutate(log10fdr = -log10(p_val_adj + eps)) %>%
  filter(!grepl(pattern_rm, gene, ignore.case = TRUE)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(sig = p_val_adj <= 0.05)

# force_gene <- c("APOE", "FTL", "AQP4", "CST3")

top_up   <- df %>% filter(sig, avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 10) %>% filter(sig, avg_log2FC >= fc_thr)
top_down <- df %>% filter(sig, avg_log2FC < 0) %>% arrange(avg_log2FC)       %>% slice_head(n = 10) %>% filter(sig, avg_log2FC <= -fc_thr)

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
  labs(title = celltype,
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



pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.OxPosOL.FocalvsDistal.Vol.pdf", width = 2, height = 2)
print(p1)
dev.off() # Close the PDF device


####################################################################################################################################################################
##################################################################### Volcano plot by celltype #####################################################################
####################################################################################################################################################################

celltype <- "OxPos OL"
tt.tmp <- tmp
Idents(tt.tmp) <- 'Annotation'
test.markers1 <- FindMarkers(tmp, 
                             ident.1 = "OxPos OL",
                             ident.2 = "OL",
                             min.pct = 0.1)

test.markers1$symbol <- rownames(test.markers1)
test.markers1$gene <- rownames(test.markers1)

test.markers1 <- test.markers1 %>% filter(!grepl("^MT-|^mt-", gene))
test.markers1 <- test.markers1 %>% filter(!grepl("^RP[SL]", gene))

fc_thr   <- log2(1.5)
eps      <- 1e-300

pattern_rm <- "^(AC|AL|AP)[0-9]+(?:\\.[0-9]+)?$|^LINC|.-AS[0-9]+$"

df <- test.markers1 %>%
  mutate(log10fdr = -log10(p_val_adj + eps)) %>%
  filter(!grepl(pattern_rm, gene, ignore.case = TRUE)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(sig = p_val_adj <= 0.05)

# force_gene <- c("APOE", "FTL", "AQP4", "CST3")

top_up   <- df %>% filter(sig, avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 10) %>% filter(sig, avg_log2FC >= fc_thr)
top_down <- df %>% filter(sig, avg_log2FC < 0) %>% arrange(avg_log2FC)       %>% slice_head(n = 10) %>% filter(sig, avg_log2FC <= -fc_thr)

force_row <- df %>%
  filter(sig) %>%
  filter(toupper(gene) %in% toupper(force_gene))

top_lab <- bind_rows(top_up, top_down, force_row) %>% distinct(gene, .keep_all = TRUE)

# mt_in_top <- top_lab$gene[grepl("^MT-", top_lab$gene)]
# mt_in_top <- top_lab$gene[grepl("^(MT-|RP[SL])", top_lab$gene)]

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
  labs(title = celltype,
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



pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.OxPosOL.Vol.pdf", width = 2, height = 2)
print(p1)
dev.off() # Close the PDF device

###############################################################################################################################################################
##################################################################### GO plot by celltype #####################################################################
###############################################################################################################################################################
# celltype <- "OxPos OL"
# tt.tmp <- subset(tmp, subset = Annotation == celltype)

# Idents(tt.tmp) <- 'condition'
# test.markers1 <- FindMarkers(tt.tmp, 
#                              ident.1 = "Focal",
#                              ident.2 = "Distal",
#                              min.pct = 0.1,
#                              test.use = "MAST",
#                              latent.vars = "patient")
# test.markers1$gene <- rownames(test.markers1)
# no significant result by condition


DefaultAssay(tmp) <- "RNA"
tmp <- JoinLayers(tmp)
Idents(tmp) <- 'Annotation'
markers <- FindAllMarkers(object = tmp, only.pos = FALSE, min.pct = 0.1, logfc.threshold = log2(1.5))


test.markers1 <- subset(markers, subset = cluster == celltype)

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

test.markers1$symbol <- test.markers1$gene
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
  head(20)


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


pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.OxPosOL.GO.pdf", width = 3, height = 2)
print(p1)
dev.off() # Close the PDF device

