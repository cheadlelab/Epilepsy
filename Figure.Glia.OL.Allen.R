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
####################################################### load Allen Dataset #######################################################
##################################################################################################################################
os_name <- Sys.info()[["sysname"]]   # "Linux", "Windows", "Darwin" …
file_path <- switch(
  os_name,
  "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output", 
  "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",
  stop(sprintf("Unsupported OS: %s", os_name))
)

load(paste0(file_path, "/step0.AllenCortical.RData"))
test@meta.data$dataset <- "Allen_Cortical"
test@meta.data$sex <- test@meta.data$donor_sex_label
test@meta.data$donor_sex_label <- NULL
test@meta.data$condition <- "Control"
Allen1 <- test
a <- rownames(Allen1)
a <- gsub("\\.AS(\\d+)", "-AS\\1", a)
a <- gsub("\\.OT", "-OT", a)
a <- gsub("\\.DT", "-DT", a)

a <- gsub("HLA\\.", "HLA-", a)
a <- gsub("MT\\.", "MT-", a)
rownames(Allen1) <- a
# obj <- subset(obj, features = genes)
annotation_mapping <- c(
  "Unknown"         = "Unknown",
  "VIP"             = "Inh_VIP",
  "LAMP5"           = "Inh_LAMP5",
  "IT"              = "Exc_L23IT",
  "PAX6"            = "Inh_CXCL14",
  "Oligodendrocyte" = "Olg",
  "Astrocyte"       = "Astrocyte",
  "L5/6 IT Car3"    = "Exc_L56IT",
  "L5/6 NP"         = "Exc_L56NP",
  "SST"             = "Inh_SST",
  "L6 CT"           = "Exc_L6CT",
  "OPC"             = "OPC",
  "PVALB"           = "Inh_PVALB",
  "L6b"             = "Exc_L6b",
  "Microglia"       = "Microglia",
  "L5 ET"           = "Exc_L5ET",
  "Pericyte"        = "Vascular",
  "Endothelial"     = "Vascular",
  "L4 IT"           = "Exc_L4IT",
  "VLMC"            = "Vascular"
)
Allen1@meta.data$Annotation <- annotation_mapping[Allen1@meta.data$subclass_label]
Allen1 <- subset(Allen1, Annotation == "Olg")


load(paste0(file_path, "/step0.AllenMTG.RData"))
test@meta.data$dataset <- "Allen_MTG"
test@meta.data$sex <- test@meta.data$donor_sex_label
test@meta.data$donor_sex_label <- NULL
test@meta.data$cluster_label_simplify <- test@meta.data$cluster_label
test@meta.data$condition <- "Control"
Allen2 <- test
a <- rownames(Allen2)
a <- gsub("\\.AS(\\d+)", "-AS\\1", a)
a <- gsub("\\.OT", "-OT", a)
a <- gsub("\\.DT", "-DT", a)

a <- gsub("HLA\\.", "HLA-", a)
a <- gsub("MT\\.", "MT-", a) 
rownames(Allen2) <- a
# obj <- subset(obj, features = genes)
annotation_mapping_2 <- c(
  "Pax6"            = "Inh_CXCL14",
  "L5/6 NP"         = "Exc_L56NP",
  "L5 IT"           = "Exc_L5IT",
  "L6 CT"           = "Exc_L6CT",
  "L4 IT"           = "Exc_L4IT",
  "Astrocyte"       = "Astrocyte",
  "L2/3 IT"         = "Exc_L23IT",
  "Vip"             = "Inh_VIP",
  "Sst Chodl"       = "Inh_SST", 
  "L6 IT Car3"      = "Exc_L56IT",
  "L6 IT"           = "Exc_L6IT",
  "Sncg"            = "Inh_CXCL14",
  "Pvalb"           = "Inh_PVALB",
  "Oligodendrocyte" = "Olg",
  "Sst"             = "Inh_SST",
  "VLMC"            = "Vascular", # Endo/Pericytes
  "Lamp5"           = "Inh_LAMP5",
  "Microglia-PVM"   = "Microglia",
  "OPC"             = "OPC",
  "L6b"             = "Exc_L6b",
  "Lamp5 Lhx6"      = "Inh_LAMP5",
  "Endothelial"     = "Vascular",
  "L5 ET"           = "Exc_L5ET",
  "Chandelier"      = "Chandelier"
)

Allen2@meta.data$Annotation <- annotation_mapping_2[Allen2@meta.data$subclass_label]
Allen2 <- subset(Allen2, Annotation == "Olg")

###################################################################################################################################
############################################################## merge ##############################################################
###################################################################################################################################
test <- merge(Allen1, y = c(Allen2))

rm(list = setdiff(ls(), "test"))
gc()

test <- NormalizeData(test)
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test)

test <- IntegrateLayers(
  object = test, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

ElbowPlot(test)

#####################################################################################################################################
############################################################## Label 1 ##############################################################
#####################################################################################################################################

dims_parameter <- 20
test <- RunUMAP(test, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
test <- FindNeighbors(test, reduction = "harmony", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.3)
test@reductions$umap_harmony <- test@reductions$umap
test@meta.data$clusters_harmony <- test@meta.data$seurat_clusters




p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "clusters_harmony", raster=FALSE, label = TRUE)
p4 <- DimPlot(test, reduction = "umap_harmony", group.by = "orig.ident", raster=FALSE, label = TRUE)
p2+p4

##################################################################################################################################
## addmodulescore
Markers.CIs <- c("NDUFA1", "NDUFA4", "NDUFA6", "NDUFB2", "NDUFB5", "NDUFS3", "NDUFS7",
                 "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2", "COX7C", "COX8A",
                 "ATP5F1A", "ATP5F1B", "ATP5MC1", "ATP5MC2", "ATP5PF")
Markers.CIs <- c("COX4I1", "COX5B", "NDUFA4", "ATP5F1B")
  
test <- AddModuleScore(test, features = list(Markers.CIs),  name = "CIs_Score")

minmax <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
test$CIs_Score  <- minmax(test$CIs_Score1)

FeaturePlot(test, features = c("CIs_Score"))

save_path <- "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5"

pdf(file = paste0(save_path, "/Figure.OL.Allen_CIs_Score.pdf"), width = 5, height = 5)
FeaturePlot(
  test, features = c("CIs_Score"),
  cols = c("lightgrey", "red"),
  raster = FALSE, pt.size = 0.5, order = TRUE,
) + NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device


