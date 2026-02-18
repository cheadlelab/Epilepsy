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

load(paste0(file_path, "/step1.0.cheadle.merge.Elzar.RData"))
# genes <- rownames(test1)
samples_df <- data.frame(
  orig.ident = c("Exp01", "Exp02_Control", "Exp02_ES", "Exp02_EU",
                 "Exp03_EP", "Exp03_HU",
                 "Exp04_3E",
                 "Exp04_4C", "Exp04_4E", "Exp04_4S",
                 "Exp05_5C", "Exp05_5E", "Exp05_5S",
                 "Exp06_6C", "Exp06_6E", "Exp06_6S",
                 "Exp07_4C", "Exp07_4E", "Exp07_5E"),
  age = c("28", "28", "28", "28",
          "39", "39",
          "24",
          "34", "34", "34",
          "59", "59", "59",
          "43", "43", "43",
          "34", "34", "59"),
  condition = c("Distal", "Distal", "Stimulated", "Focal", 
                "Focal", "Distal", 
                "Focal",
                "Distal", "Focal", "Stimulated",
                "Distal", "Focal", "Stimulated",
                "Distal", "Focal", "Stimulated",
                "Distal", "Focal", "Focal"),
  sex = c("F", "F", "F", "F",
          "F", "F",
          "M",
          "M", "M", "M",
          "F", "F", "F",
          "F", "F", "F",
          "M", "M", "F"),
  patient = c("1", "1", "1", "1",
              "2", "2",
              "3",
              "4", "4", "4",
              "5", "5", "5",
              "6", "6", "6",
              "4", "4", "5")
)
test@meta.data$age <- samples_df$age[match(test@meta.data$orig.ident, samples_df$orig.ident)]
test@meta.data$condition <- samples_df$condition[match(test@meta.data$orig.ident, samples_df$orig.ident)]
test@meta.data$sex <- samples_df$sex[match(test@meta.data$orig.ident, samples_df$orig.ident)]
test@meta.data$patient <- samples_df$patient[match(test@meta.data$orig.ident, samples_df$orig.ident)]

test1 <- test
test1@meta.data$dataset <- "Cheadle"

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
Allen1 <- subset(Allen1, Annotation != "Unknown")


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

###################################################################################################################################
############################################################## merge ##############################################################
###################################################################################################################################
test <- merge(test1, y = c(Allen1, Allen2))

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

##### For scvi method #####
test <- IntegrateLayers(
  object = test, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/grid/cheadle/home/qianyu/.conda/envs/scvi-env", verbose = FALSE
)


# test <- JoinLayers(test)
pdf(file = paste0(file_path, "/elbow.pdf"), width = 10, height = 10)
ElbowPlot(test)
dev.off() # Close the PDF device

# save(test, file = paste0(file_path, "/step1.0.annotation.AC.tmp.RData"))
# load(paste0(file_path, "/step1.0.annotation.AC.tmp.RData"))

#####################################################################################################################################
############################################################## Label 1 ##############################################################
#####################################################################################################################################

dims_parameter <- 20
test <- RunUMAP(test, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
test <- FindNeighbors(test, reduction = "harmony", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.3)
test@reductions$umap_harmony <- test@reductions$umap
test@meta.data$clusters_harmony <- test@meta.data$seurat_clusters


dims_parameter <- 20
test <- RunUMAP(test, reduction = "integrated.scvi", dims = 1:dims_parameter, reduction.key = "scvi_")
test <- FindNeighbors(test, reduction = "integrated.scvi", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.2)
test@reductions$umap_scvi <- test@reductions$umap
test@meta.data$clusters_scvi <- test@meta.data$seurat_clusters

# plot cluster_label_simplify
pdf(file = paste0(file_path, "/harmony.20.pdf"), width = 30, height = 15)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "clusters_harmony", raster=FALSE, label = TRUE)
# p2 <- DimPlot(test, reduction = "umap_scvi", group.by = "clusters_scvi", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "subclass_label", raster=FALSE, label = TRUE)
p1+p2
dev.off() # Close the PDF device

file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output"
# save(test, file = paste0(file_path, "/step1.0.annotation.AC.tmp.RData"))

table(test@meta.data$seurat_clusters, test@meta.data$subclass_label)

Idents(test) <- "clusters_harmony"

new_identities <- c("0" = "Exc_L23",
                    "1" = "Olg",
                    "2" = "Exc_L5",
                    "3" = "Inh_VIP",
                    "4" = "Inh_SST",
                    "5" = "Exc_L4",
                    "6" = "Inh_PVALB",
                    "7" = "Astrocyte",
                    "8" = "Exc_L6",
                    "9" = "Exc_L23",
                    "10" = "Inh_LAMP5",
                    "11" = "OPC",
                    "12" = "Exc_L56IT_CAR3",
                    "13" = "Exc_L6b",
                    "14" = "Exc_L6CT",
                    "15" = "Exc_L5",
                    "16" = "Exc_L56NP",
                    "17" = "Microglia",
                    "18" = "Vascular",
                    "19" = "Chandelier",
                    "20" = "Exc_New",
                    "21" = "Exc_L5ET",
                    "22" = "Olg",
                    "23" = "Olg"
)

test@meta.data$Annotation <- new_identities[test@meta.data$clusters_harmony]

# plot cluster_label_simplify
pdf(file = paste0(file_path, "/harmony.20.pdf"), width = 30, height = 15)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation", raster=FALSE, label = TRUE)
# p2 <- DimPlot(test, reduction = "umap_scvi", group.by = "clusters_scvi", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "subclass_label", raster=FALSE, label = TRUE)
p1+p2
dev.off() # Close the PDF device

##################################################################################################################################
############################################################## save ##############################################################
##################################################################################################################################
pdf(file = paste0(file_path, "/umap_label.pdf"), width = 30, height = 30)
p1 <- DimPlot(test, reduction = "umap", group.by = "Annotation", label = TRUE, raster = FALSE)
p2 <- DimPlot(test, reduction = "umap", group.by = "subclass_label", label = TRUE, raster = FALSE)
p3 <- DimPlot(test, reduction = "umap", group.by = "Annotation", label = TRUE, raster = FALSE)
p4 <- DimPlot(test, reduction = "umap", group.by = "condition", label = TRUE, raster = FALSE)
p1+p2+p3+p4
dev.off() # Close the PDF device


table(test@meta.data$orig.ident, test@meta.data$Annotation)
table(test@meta.data$Annotation, test@meta.data$condition)


save(test, file = paste0(file_path, "/step1.1.annotation.AC.RData"))

cell_annotation <- test@meta.data
save(cell_annotation, file = paste0(file_path, "/step1.1.annotation.meta.RData"))



