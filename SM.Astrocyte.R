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

sheet1_df <- markers %>%
  filter(cluster %in% c("Reactive Astrocyte", "Lipid-Accumulated Reactive Astrocyte", "Other Astrocyte")) %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC))

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

sheet2_df <- df %>%
  filter(sig) %>%
  mutate(symbol = gene) %>%
  arrange(desc(avg_log2FC))
########################################################################################################################################################
Idents(tt.tmp) <- 'condition'
test.markers1 <- FindMarkers(tt.tmp, 
                             ident.1 = "Focal",
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

sheet3_df <- df %>%
  filter(sig) %>%
  mutate(symbol = gene) %>%
  arrange(desc(avg_log2FC))

## ---------- Write Excel ----------
out_xlsx <- "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V7/SM.DEG.Astrocyte.xlsx"

wb <- createWorkbook()

wb <- createWorkbook()

addWorksheet(wb, "RA_LARA_Other_DEG")
writeData(wb, "RA_LARA_Other_DEG", sheet1_df)

addWorksheet(wb, "RA_Focal_vs_Distal")
writeData(wb, "RA_Focal_vs_Distal", sheet2_df)

addWorksheet(wb, "RA_Stim_vs_Distal")
writeData(wb, "RA_Stim_vs_Distal", sheet3_df)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved: ", out_xlsx)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved: ", out_xlsx)