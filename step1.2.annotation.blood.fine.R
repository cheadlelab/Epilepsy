# https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC
rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)

file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"
load(paste0(file_path, "/step1.1.umap.coarse.RData"))

test@meta.data$Annotation.fine <- "Unknown"


#####################################################################################################################################
########################################################### Label Myeloid ###########################################################
#####################################################################################################################################
tmp <- subset(test, subset = Annotation.coarse %in% c("Myeloid"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 25
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.2)


# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

# a <- FindMarkers(tmp, ident.1 = 13, ident.2 = NULL, min.pct = 0.1, only.pos = TRUE)
a <- FindAllMarkers(tmp, min.pct = 0.1, only.pos = TRUE)
top <- a %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))


markers <- c("PPP1R14A", "CLEC9A", "FCER1A", "ITM2C", "CD14", "FCGR3A", "GNLY", "CD3D", "PPBP", "MS4A1")
FeaturePlot(tmp, features = markers)


Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "CD14_Mono", # CD14
                    "1" = "CD16_Mono", # FCGR3A
                    "2" = "CD14_Mono", # CD14
                    "3" = "CD14_Mono", # CD14
                    "4" = "CD14_Mono", # CD14
                    "5" = "CD14_Mono", # CD14
                    "6" = "Myeloid/NKT", # CD14
                    "7" = "cDC2", # FCER1A # some B/cDC2
                    "8" = "Myeloid/platelets", # CD14 TUBB1
                    "9" = "cDC1" # CLEC9A
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)


# B/cDC2
tmp <- subset(test, subset = Annotation.fine %in% c("cDC2"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 25
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.2)


# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)
FeaturePlot(tmp, features = "MS4A1")

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "cDC2",
                    "1" = "cDC2", 
                    "2" = "cDC2",
                    "3" = "B/cDC2"
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)

###############################################################################################################################
########################################################### Label B ###########################################################
###############################################################################################################################
tmp <- subset(test, subset = Annotation.coarse %in% c("B"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 15
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.2)

# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

# a <- FindMarkers(tmp, ident.1 = 13, ident.2 = NULL, min.pct = 0.1, only.pos = TRUE)
a <- FindAllMarkers(tmp, min.pct = 0.1, only.pos = TRUE)
top <- a %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

markers <- c("PPP1R14A", "CLEC9A", "FCER1A", "ITM2C", "CD14", "FCGR3A", "GNLY", "CD3D", "PPBP", "MS4A1")
FeaturePlot(tmp, features = markers)

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "B", # 
                    "1" = "B", # 
                    "2" = "B", # 
                    "3" = "B/Myeloid", 
                    "4" = "B/NKT", 
                    "5" = "B/platelets", # 
                    "6" = "B"
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)

# second round
tmp <- subset(test, subset = Annotation.fine %in% c("B"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 15
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.3)

# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

# a <- FindMarkers(tmp, ident.1 = 13, ident.2 = NULL, min.pct = 0.1, only.pos = TRUE)
a <- FindAllMarkers(tmp, min.pct = 0.1, only.pos = TRUE)
top <- a %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

markers <- c("LINC01857", "TNFRSF13B", "IGHM", "COCH", "AIM2", "IGHM", "IGHD")
FeaturePlot(tmp, features = markers)

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "B_naive", # 
                    "1" = "B_intermediate", # LINC01857
                    "2" = "B_naive", # 
                    "3" = "B_naive", 
                    "4" = "B_memory", # COCH
                    "5" = "B_intermediate", # LINC01857
                    "6" = "B/other"
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)

###################################################################################################################################
########################################################### Label Other ###########################################################
###################################################################################################################################
tmp <- subset(test, subset = Annotation.coarse %in% c("Other"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 20
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.1)

# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

# a <- FindMarkers(tmp, ident.1 = 13, ident.2 = NULL, min.pct = 0.1, only.pos = TRUE)
a <- FindAllMarkers(tmp, min.pct = 0.1, only.pos = TRUE)
top <- a %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

markers <- c("MZB1")
FeaturePlot(tmp, features = markers)
FeaturePlot(tmp, features = "CPA3")

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "MKI67+_NKT", # MKI67
                    "1" = "EMP", # Erythroid Megakaryocyte Progenitor # MYCT1
                    "2" = "Plasmablast", # MZB1
                    "3" = "pro_B", # DNTT
                    "4" = "Late_Eryth", # CTSE
                    "5" = "Mast" # TPSAB1
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)

#################################################################################################################################
########################################################### Label NKT ###########################################################
#################################################################################################################################
tmp <- subset(test, subset = Annotation.coarse %in% c("NKT"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 20
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.3)

# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

# a <- FindMarkers(tmp, ident.1 = 13, ident.2 = NULL, min.pct = 0.1, only.pos = TRUE)
a <- FindAllMarkers(tmp, min.pct = 0.1, only.pos = TRUE)
top <- a %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

markers <- c("CCR7")
FeaturePlot(tmp, features = markers)
FeaturePlot(tmp, features = "KLRD1")

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "TCM_TEM_CTL",
                    "1" = "TCM_TEM_CTL",
                    "2" = "TCM_TEM_CTL",
                    "3" = "NK", # FCER1G
                    "4" = "CD4_Naive", # CCR7
                    "5" = "gdT", # TRDV2
                    "6" = "CD8_Naive", # CCR7
                    "7" = "Treg", # FOXP3
                    "8" = "NK_CD56bright", # NCAM1
                    "9" = "MAIT", # SLC4A10
                    "10" = "Platelets" # PPBP
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)

# second round
tmp <- subset(test, subset = Annotation.fine %in% c("TCM_TEM_CTL"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 20
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.3)

# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)

# a <- FindMarkers(tmp, ident.1 = 13, ident.2 = NULL, min.pct = 0.1, only.pos = TRUE)
a <- FindAllMarkers(tmp, min.pct = 0.1, only.pos = TRUE)
top <- a %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

FeaturePlot(tmp, features = c("CD4","CD8A"))

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "CD4_Memory", # TMSB10
                    "1" = "Memory",
                    "2" = "CD8_Memory",
                    "3" = "CD8_Memory", 
                    "4" = "NK",
                    "5" = "CD4_CTL" # GZMH
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)

# third round
tmp <- subset(test, subset = Annotation.fine %in% c("Memory"))

tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp)
tmp <- ScaleData(tmp)
tmp <- RunPCA(tmp)

ElbowPlot(tmp)
dims_parameter <- 15
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:dims_parameter)
tmp <- FindClusters(tmp, resolution = 0.2)

# plot cluster_label_simplify
DimPlot(tmp, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)


FeaturePlot(tmp, features = c("CD4","CD8A"))

Idents(tmp) <- "seurat_clusters"

new_identities <- c("0" = "CD8_Memory",
                    "1" = "CD8_Memory",
                    "2" = "CD4_Memory",
                    "3" = "CD4_Memory"
)

tmp@meta.data$Annotation.fine <- new_identities[tmp@meta.data$seurat_clusters]
DimPlot(tmp, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tmp))
new_annotations <- tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations
test@meta.data$Annotation.fine[test@meta.data$Annotation.fine == "Unknown"] <- test@meta.data$Annotation.coarse[test@meta.data$Annotation.fine == "Unknown"]
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)

table(test@meta.data$Annotation.fine, test@meta.data$condition)
table(test@meta.data$Annotation.fine, test@meta.data$patient)

save(test, file = paste0(file_path, "/step1.2.umap.fine.RData"))





`%!in%` <- Negate(`%in%`)
test <- subset(test, subset = Annotation.fine %!in% c("B/cDC2", "B/Myeloid", "B/NKT", "B/other", "B/platelets", "Myeloid/NKT", "Myeloid/platelets"))

DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster=FALSE, label = TRUE)







######################################################################################################################################
################################################################ UMAP ################################################################
######################################################################################################################################
save_path <- file_path
pdf(file = paste0(save_path, "/Figure1.blood.annotation.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.blood.annotation.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", raster = FALSE)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.blood.condition.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap", group.by = "condition", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.blood.condition.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap", group.by = "condition", label = TRUE, raster = FALSE)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.blood.patient.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap", group.by = "patient", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.blood.patient.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap", group.by = "patient", label = TRUE, raster = FALSE)
dev.off() # Close the PDF device