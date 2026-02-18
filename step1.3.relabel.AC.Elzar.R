# scvi envs
rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)

library(reticulate)
use_condaenv("/grid/cheadle/home/qianyu/.conda/envs/scvi-env", required = TRUE)

##################################################################################################################################
###################################################### load Cheadle Dataset ######################################################
##################################################################################################################################
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output"
# file_path <- "/grid/Epilepsy/output"

load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))

test <- JoinLayers(test)

test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
test <- subset(test, subset = Annotation.fine != c("Other"))
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay5")
test[["RNA"]] <- split(test[["RNA"]], f = test$orig.ident)

test@meta.data$donor <- test@meta.data$patient
test1 <- test
##################################################################################################################################
####################################################### load Allen Dataset #######################################################
##################################################################################################################################

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
  "L5/6 IT Car3"    = "Exc_L56IT_CAR3",
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
Allen1@meta.data$Annotation.coarse <- annotation_mapping[Allen1@meta.data$subclass_label]
Allen1 <- subset(Allen1, subset = Annotation.coarse != c("Unknown"))
Allen1@meta.data$Annotation.fine <- Allen1@meta.data$Annotation.coarse
Allen1@meta.data$patient <- substr(Allen1@meta.data$sample_name, 1, 11)



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
  "L6 IT Car3"      = "Exc_L56IT_CAR3",
  "L6 IT"           = "Exc_L6IT",
  "Sncg"            = "Inh_CXCL14",
  "Pvalb"           = "Inh_PVALB",
  "Oligodendrocyte" = "Olg",
  "Sst"             = "Inh_SST",
  "VLMC"            = "Vascular",
  "Lamp5"           = "Inh_LAMP5",
  "Microglia-PVM"   = "Microglia",
  "OPC"             = "OPC",
  "L6b"             = "Exc_L6b",
  "Lamp5 Lhx6"      = "Inh_LAMP5_LHX6",
  "Endothelial"     = "Vascular",
  "L5 ET"           = "Exc_L5ET",
  "Chandelier"      = "Chandelier"
)

Allen2@meta.data$Annotation.coarse <- annotation_mapping_2[Allen2@meta.data$subclass_label]
Allen2@meta.data$Annotation.fine <- Allen2@meta.data$Annotation.coarse
Allen2@meta.data$patient <- Allen2@meta.data$sample_name


###################################################################################################################################
############################################################## merge ##############################################################
###################################################################################################################################
test <- merge(test1, y = c(Allen1, Allen2))

rm(list = setdiff(ls(), "test"))
gc()
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output"

test <- NormalizeData(test)
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test)

test <- IntegrateLayers(
  object = test, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

##### For scvi method #####
# test <- IntegrateLayers(
#   object = test, method = scVIIntegration,
#   new.reduction = "integrated.scvi",
#   conda_env = "/grid/cheadle/home/qianyu/.conda/envs/scvi-env", verbose = FALSE
# )


# test <- JoinLayers(test)
pdf(file = paste0(file_path, "/elbow.pdf"), width = 10, height = 10)
ElbowPlot(test)
dev.off() # Close the PDF device

# save(test, file = paste0(file_path, "/step1.0.annotation.AC.tmp.RData"))
# load(paste0(file_path, "/step1.0.annotation.AC.tmp.RData"))


##################################################################################################################################
############################################################## umap ##############################################################
##################################################################################################################################

dims_parameter <- 20
test <- RunUMAP(test, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
test <- FindNeighbors(test, reduction = "harmony", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.2)
test@reductions$umap_harmony <- test@reductions$umap
test@meta.data$clusters_harmony <- test@meta.data$seurat_clusters

######################################################################################################################################
############################################################## pre save ##############################################################
######################################################################################################################################
pdf(file = paste0(file_path, "/umap_all.pdf"), width = 30, height = 30)
p1 <- DimPlot(test, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(test, reduction = "umap", group.by = "subclass_label", label = TRUE, raster = FALSE)
p3 <- DimPlot(test, reduction = "umap", group.by = "Annotation.fine", label = TRUE, raster = FALSE)
p4 <- DimPlot(test, reduction = "umap", group.by = "condition", label = TRUE, raster = FALSE)
p1+p2+p3+p4
dev.off() # Close the PDF device


table(test@meta.data$Annotation.fine, test@meta.data$condition)

save(test, file = paste0(file_path, "/step1.3.annotation.AC.RData"))
