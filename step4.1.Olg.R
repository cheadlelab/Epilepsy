rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)

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



DefaultAssay(tmp) <- "RNA"
tmp <- JoinLayers(tmp)
Idents(tmp) <- 'seurat_clusters'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

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
save(tmp, file = paste0(file_path, "/step4.Olg.RData"))
load(paste0(file_path, "/step4.Olg.RData"))

tmp$new <- paste0(tmp$patient,"_",tmp$condition)
table(tmp$new,tmp$Annotation)






Idents(tmp) <- 'Annotation'
markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DimPlot(tmp, reduction = "umap_harmony", group.by = "Annotation", raster=FALSE, label = TRUE)


FeaturePlot(tmp, features = c("SLC38A1", "APOE","FTL","TFRC"))



tt.tmp <- subset(tmp, subset = Annotation %in% c("Reactive_Astrocytes"))
Idents(tt.tmp) <- 'condition'
markers <- FindAllMarkers(object = tt.tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
