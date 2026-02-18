rm(list = ls(all = TRUE))

#load packages
library(Seurat)
library(dplyr)
library(jsonlite)
library(dplyr)
library(arrow)
library(spacexr)
library(ggplot2)
library(scales)
options(future.globals.maxSize = 128 * 1024^3)

save.path <- "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/"

load(paste0(save.path, "step2.umap.RData"))


orig.ident <- c("1C SA", "4E", "1C SB", "5C", "2E", "3E SB", "1E", "4C", "5E", "2C", "3E SA")
patient <- c(1, 4 ,1, 5, 2, 3, 1, 4, 5, 2, 3)
condition <- c("Distal", "Focal", "Distal", "Distal", "Focal", "Focal", "Focal", "Distal", "Focal", "Distal", "Focal")

test@meta.data$patient <- patient[match(test@meta.data$orig.ident, orig.ident)]
test@meta.data$condition <- condition[match(test@meta.data$orig.ident, orig.ident)]




tmp <- subset(test, subset = predicted.celltype %in% c("Olg"))

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

dims_parameter <- 15
tmp <- RunUMAP(tmp, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
tmp <- FindNeighbors(tmp, reduction = "harmony", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.3)
tmp@reductions$umap_harmony <- tmp@reductions$umap

DimPlot(tmp, reduction = "umap_harmony", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

FeaturePlot(tmp, features = c("IDH1","IDH2","MGST1","EFHD1"))
VlnPlot(tmp, features = c("IDH1","IDH2","MGST1","EFHD1"),pt.size=0)+DimPlot(tmp, reduction = "umap_harmony", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "OL",
                    "1" = "OL",
                    "2" = "OxPos OL",
                    "3" = "OL",
                    "4" = "OxPos OL",
                    "5" = "OL",
                    "6" = "OxPos OL",
                    "7" = "OL",
                    "8" = "OxPos OL"
)

tmp@meta.data$Annotation <- new_identities[as.character(tmp@meta.data$seurat_clusters)]

Idents(tmp) <- "Annotation"

VlnPlot(tmp, features = c("IDH1","IDH2","MGST1","EFHD1"),pt.size=0)+DimPlot(tmp, reduction = "umap_harmony", group.by = "Annotation", raster=FALSE, label = TRUE)

Idents(tmp) <- "Annotation"
tmp <- JoinLayers(tmp)
markers <- FindMarkers(tmp, 
                       ident.1 = "OxPos OL",
                       ident.2 = "OL",
                       logfc.threshold = 0,
                       min.pct = 0.1)
