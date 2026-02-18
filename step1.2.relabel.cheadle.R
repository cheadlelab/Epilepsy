# MOSS root+scvi envs
rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)

##################################################################################################################################
###################################################### load Cheadle Dataset ######################################################
##################################################################################################################################
os_name <- Sys.info()[["sysname"]]   # returns "Linux", "Windows", "Darwin" (macOS), etc.

file_path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output",   # Ubuntu and other Linux
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",  # Windows UNC path
                    # fallback – change or add more cases if you need them
                    stop(sprintf("Unsupported OS: %s", os_name))
)

load(paste0(file_path, "/step1.0.cheadle.merge.Elzar.RData"))
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
test@meta.data$dataset <- "Cheadle"


load(paste0(file_path, "/step1.1.annotation.meta.RData"))

cells_in_both <- which(colnames(test) %in% rownames(cell_annotation))
new_annotations <- cell_annotation$Annotation[match(colnames(test)[cells_in_both], rownames(cell_annotation))]
test@meta.data$Annotation.coarse[cells_in_both] <- new_annotations

VlnPlot(test, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal", "log10GenesPerUMI"),
        ncol = 5, group.by = "Annotation.coarse", pt.size = 0)

test <- JoinLayers(test)

test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
test <- subset(test, subset = Annotation.coarse != "Exc_New")
# Exc_New are low quality cells
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay5")
test[["RNA"]] <- split(test[["RNA"]], f = test$orig.ident)



##################################################################################################################################
######################################################## harmony and scvi ########################################################
##################################################################################################################################

test <- NormalizeData(test)
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test)

test <- IntegrateLayers(
  object = test, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# could skip
test <- IntegrateLayers(
  object = test, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/home/qianyu/anaconda3/envs/scvi-env", verbose = FALSE
)



dims_parameter <- 30
test <- RunUMAP(test, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
test <- FindNeighbors(test, reduction = "harmony", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.2)
test@reductions$umap_harmony <- test@reductions$umap
test@meta.data$clusters_harmony <- test@meta.data$seurat_clusters


dims_parameter <- 30
test <- RunUMAP(test, reduction = "integrated.scvi", dims = 1:dims_parameter, reduction.key = "scvi_")
test <- FindNeighbors(test, reduction = "integrated.scvi", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.2)
test@reductions$umap_scvi <- test@reductions$umap
test@meta.data$clusters_scvi <- test@meta.data$seurat_clusters

pdf(file = paste0(file_path, "/umap_relabel.pdf"), width = 30, height = 30)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "clusters_harmony", raster=FALSE, label = TRUE)
p3 <- DimPlot(test, reduction = "umap_harmony", group.by = "condition", raster=FALSE, label = TRUE)
p4 <- DimPlot(test, reduction = "umap_harmony", group.by = "patient", raster=FALSE, label = TRUE)
p1+p2+p3+p4
dev.off() # Close the PDF device


p1 <- DimPlot(test, reduction = "umap_scvi", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_scvi", group.by = "clusters_scvi", raster=FALSE, label = TRUE)
p3 <- DimPlot(test, reduction = "umap_scvi", group.by = "condition", raster=FALSE, label = TRUE)
p4 <- DimPlot(test, reduction = "umap_scvi", group.by = "patient", raster=FALSE, label = TRUE)
p1+p2+p3+p4

#################################################################################################################################
############################################################ label 1 ############################################################
#################################################################################################################################
test@meta.data$Annotation.fine <- "Empty"
# test <- JoinLayers(test)
DefaultAssay(test) <- "RNA"
Idents(test) <- 'clusters_harmony'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers <- FindAllMarkers(object = test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
FeaturePlot(test, features = c("GAD1"), reduction = "umap_harmony")

p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "clusters_harmony", raster=FALSE, label = TRUE)
p1+p2

Idents(test) <- "clusters_harmony"

new_identities <- c("0" = "Olg", # MBP
                    "1" = "Exc",
                    "2" = "Astrocyte", # OBI1-AS1
                    "3" = "Olg", # MBP
                    "4" = "Exc",
                    "5" = "Inh", # GAD1
                    "6" = "OPC", # PDGFRA
                    "7" = "Exc",
                    "8" = "Inh", # GAD1
                    "9" = "Microglia", # TLR2
                    "10" = "Inh", # GAD1
                    "11" = "Exc",
                    "12" = "Exc",
                    "13" = "Inh", # GAD1
                    "14" = "Exc",
                    "15" = "Exc",
                    "16" = "Exc",
                    "17" = "Vascular", # EBF1
                    "18" = "NK/T", # PTPRC/CD96
                    "19" = "Exc"
)

test@meta.data$Annotation.fine <- new_identities[test@meta.data$clusters_harmony]

pdf(file = paste0(file_path, "/umap_tmp.pdf"), width = 30, height = 15)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
p1+p2
dev.off() # Close the PDF device

###################################################################################################################################################
#################################################################### Label Mg #####################################################################
###################################################################################################################################################
tt.tmp <- test
tt.tmp <- JoinLayers(tt.tmp)

tt.tmp[["RNA"]] <- as(object = tt.tmp[["RNA"]], Class = "Assay")
tt.tmp <- subset(tt.tmp, subset = Annotation.fine %in% c("Microglia"))
tt.tmp[["RNA"]] <- as(object = tt.tmp[["RNA"]], Class = "Assay5")
tt.tmp[["RNA"]] <- split(tt.tmp[["RNA"]], f = tt.tmp$orig.ident)

tt.tmp <- NormalizeData(tt.tmp)
tt.tmp <- FindVariableFeatures(tt.tmp)
tt.tmp <- ScaleData(tt.tmp)
tt.tmp <- RunPCA(tt.tmp)


tt.tmp <- IntegrateLayers(
  object = tt.tmp, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


ElbowPlot(tt.tmp)

dims_parameter <- 5
tt.tmp <- RunUMAP(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindNeighbors(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindClusters(tt.tmp, resolution = 0.1)

pdf(file = paste0(file_path, "/umap_tmp1.pdf"), width = 30, height = 15)
p1 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(tt.tmp, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
p1+p2
dev.off() # Close the PDF device

pdf(file = paste0(file_path, "/umap_tmp2.pdf"), width = 30, height = 20)
FeaturePlot(tt.tmp, features = c("MBP", "TLR2", "CD163","TMEM119", "AIF1", "P2RY12"), reduction = "umap")
dev.off() # Close the PDF device

pdf(file = paste0(file_path, "/umap_tmp3.pdf"), width = 30, height = 30)
FeaturePlot(tt.tmp, features = c("AIF1", "TMEM119", "P2RY12", "CX3CR1", "CD68", "ITGAM", "TREM2", "SALL1", "FCRLS"), reduction = "umap")
dev.off() # Close the PDF device

DefaultAssay(tt.tmp) <- "RNA"
tt.tmp <- JoinLayers(tt.tmp)
Idents(tt.tmp) <- 'seurat_clusters'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers <- FindAllMarkers(object = tt.tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
FeaturePlot(tt.tmp, features = c("TLR2","MBP"), reduction = "umap")

new_identities <- c("0" = "Microglia",
                    "1" = "Macrophage",
                    "2" = "Macrophage",
                    "3" = "Olg",
                    "4" = "Other")

tt.tmp@meta.data$Annotation.fine <- new_identities[tt.tmp@meta.data$seurat_clusters]

pdf(file = paste0(file_path, "/umap_myeloid.pdf"), width = 30, height = 15)
p1 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.fine", label = TRUE, raster = FALSE)
p1+p2
dev.off() # Close the PDF device

### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tt.tmp))
new_annotations <- tt.tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tt.tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations

pdf(file = paste0(file_path, "/umap_tmp.pdf"), width = 30, height = 15)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
p1+p2
dev.off() # Close the PDF device

###################################################################################################################################################
#################################################################### Label Inh ####################################################################
###################################################################################################################################################
tt.tmp <- test
tt.tmp <- JoinLayers(tt.tmp)

tt.tmp[["RNA"]] <- as(object = tt.tmp[["RNA"]], Class = "Assay")
tt.tmp <- subset(tt.tmp, subset = Annotation.fine %in% c("Inh"))
tt.tmp[["RNA"]] <- as(object = tt.tmp[["RNA"]], Class = "Assay5")
tt.tmp[["RNA"]] <- split(tt.tmp[["RNA"]], f = tt.tmp$orig.ident)

tt.tmp <- NormalizeData(tt.tmp)
tt.tmp <- FindVariableFeatures(tt.tmp)
tt.tmp <- ScaleData(tt.tmp)
tt.tmp <- RunPCA(tt.tmp)


tt.tmp <- IntegrateLayers(
  object = tt.tmp, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

ElbowPlot(tt.tmp)

dims_parameter <- 20
tt.tmp <- RunUMAP(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindNeighbors(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindClusters(tt.tmp, resolution = 0.1)

pdf(file = paste0(file_path, "/umap_Inh.pdf"), width = 30, height = 15)
p1 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(tt.tmp, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
p1+p2
dev.off() # Close the PDF device


DefaultAssay(tt.tmp) <- "RNA"
tt.tmp <- JoinLayers(tt.tmp)
Idents(tt.tmp) <- 'seurat_clusters'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers <- FindAllMarkers(object = tt.tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf(file = paste0(file_path, "/umap_tmp4.pdf"), width = 10, height = 10)
FeaturePlot(tt.tmp, features = c("CXCL14", "LHX6", "VIP"), reduction = "umap")
dev.off() # Close the PDF device

new_identities <- c(
  "0" = "Inh_SST",
  "1" = "Inh_VIP",
  "2" = "Inh_PVALB",
  "3" = "Inh_CXCL14",
  "4" = "Inh_VIP",
  "5" = "Inh_LAMP5_LHX6",
  "6" = "Inh_LAMP5",
  "7" = "Chandelier",
  "8" = "Other"
)

tt.tmp@meta.data$Annotation.fine <- new_identities[tt.tmp@meta.data$seurat_clusters]

pdf(file = paste0(file_path, "/umap_tmp5.pdf"), width = 20, height = 10)
p1 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.fine", label = TRUE, raster = FALSE)
p1+p2
dev.off() # Close the PDF device

### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tt.tmp))
new_annotations <- tt.tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tt.tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations

pdf(file = paste0(file_path, "/umap_tmp.pdf"), width = 20, height = 10)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
p1+p2
dev.off() # Close the PDF device


save(test, file = paste0(file_path, "/tmp1.RData"))
load(paste0(file_path, "/tmp1.RData"))

###################################################################################################################################################
#################################################################### Label Exc ####################################################################
###################################################################################################################################################
tt.tmp <- test
tt.tmp <- JoinLayers(tt.tmp)

tt.tmp[["RNA"]] <- as(object = tt.tmp[["RNA"]], Class = "Assay")
tt.tmp <- subset(tt.tmp, subset = Annotation.fine %in% c("Exc"))
tt.tmp[["RNA"]] <- as(object = tt.tmp[["RNA"]], Class = "Assay5")
tt.tmp[["RNA"]] <- split(tt.tmp[["RNA"]], f = tt.tmp$orig.ident)

tt.tmp <- NormalizeData(tt.tmp)
tt.tmp <- FindVariableFeatures(tt.tmp)
tt.tmp <- ScaleData(tt.tmp)
tt.tmp <- RunPCA(tt.tmp)


tt.tmp <- IntegrateLayers(
  object = tt.tmp, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

ElbowPlot(tt.tmp)

dims_parameter <- 30
tt.tmp <- RunUMAP(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindNeighbors(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindClusters(tt.tmp, resolution = 0.2)

pdf(file = paste0(file_path, "/umap_Exc.pdf"), width = 20, height = 10)
p1 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(tt.tmp, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
p1+p2
dev.off() # Close the PDF device

DefaultAssay(tt.tmp) <- "RNA"
tt.tmp <- JoinLayers(tt.tmp)
Idents(tt.tmp) <- 'seurat_clusters'
#markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
markers <- FindAllMarkers(object = tt.tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
save(top10, file = paste0(file_path, "/top10.RData"))
write.csv(top10, file = paste0(file_path, "/top10.csv"), row.names = FALSE)

pdf(file = paste0(file_path, "/umap_tmp6.pdf"), width = 20, height = 20)
FeaturePlot(tt.tmp, features = c("GRIN3A", "PLCH1", "RORB", "CUX2"), reduction = "umap")
dev.off() # Close the PDF device

FeaturePlot(tt.tmp, features = c("RORB", "CUX2","THEMIS","FEZF2"), reduction = "umap")
FeaturePlot(tt.tmp, features = c("ERBB4","BCAS1","RORB","CUX2"), reduction = "umap")

new_identities <- c("0" = "Exc_L23IT",
                    "1" = "Exc_L5IT_GRIN3A",
                    "2" = "Exc_L23IT",
                    "3" = "Exc_L234IT_ERBB4",
                    "4" = "Exc_L4IT",
                    "5" = "Exc_L6IT",
                    "6" = "Exc_L6b",
                    "7" = "Exc_L6CT",
                    "8" = "Exc_L56IT_CAR3",
                    "9" = "Exc_L56NP",
                    "10" = "Exc_L5IT_GRIN3A-",
                    "11" = "Exc_L4IT_PLCH1",
                    "12" = "Exc_L5ET",
                    "13" = "Other"
)


tt.tmp@meta.data$Annotation.fine <- new_identities[tt.tmp@meta.data$seurat_clusters]

pdf(file = paste0(file_path, "/umap_tmp.pdf"), width = 20, height = 10)
p1 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(tt.tmp, reduction = "umap", group.by = "Annotation.fine", label = TRUE, raster = FALSE)
p1+p2
dev.off() # Close the PDF device

### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(tt.tmp))
new_annotations <- tt.tmp@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(tt.tmp))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations

pdf(file = paste0(file_path, "/umap_Exc.pdf"), width = 20, height = 10)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
p1+p2
dev.off() # Close the PDF device

###############################################################################################
file_path <- "/home/qianyu/Desktop"
save(test, file = paste0(file_path, "/step1.2.relabel.cheadle.RData"))
load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))

########################################################################################################################################################
#################################################################### Label Mg again ####################################################################
########################################################################################################################################################


Myeloid <- subset(test, subset = Annotation.fine %in% c("Macrophage", "Microglia"))

Myeloid <- JoinLayers(Myeloid)
Myeloid[["RNA"]] <- split(Myeloid[["RNA"]], f = Myeloid$orig.ident)

Myeloid <- NormalizeData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid)
Myeloid <- ScaleData(Myeloid)
Myeloid <- RunPCA(Myeloid)

Myeloid <- IntegrateLayers(
  object = Myeloid, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

ElbowPlot(Myeloid)

dims_parameter <- 5
Myeloid <- RunUMAP(Myeloid, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
Myeloid <- FindNeighbors(Myeloid, reduction = "harmony", dims = 1:dims_parameter)
Myeloid <- FindClusters(Myeloid, resolution = 0.2)
Myeloid@reductions$umap_harmony <- Myeloid@reductions$umap
Myeloid@meta.data$clusters_harmony <- Myeloid@meta.data$seurat_clusters


new_identities <- c("0" = "Microglia",
                    "1" = "Microglia",
                    "2" = "Macrophage",
                    "3" = "Macrophage",
                    "4" = "Macrophage"
)

Myeloid@meta.data$Annotation.fine <- new_identities[Myeloid@meta.data$seurat_clusters]

p1 <- DimPlot(Myeloid, reduction = "umap", group.by = "Annotation.coarse", label = TRUE, raster = FALSE)
p2 <- DimPlot(Myeloid, reduction = "umap", group.by = "Annotation.fine", label = TRUE, raster = FALSE)
p1+p2

### pass the labels
cells_in_both <- which(colnames(test) %in% colnames(Myeloid))
new_annotations <- Myeloid@meta.data$Annotation.fine[match(colnames(test)[cells_in_both], colnames(Myeloid))]
test@meta.data$Annotation.fine[cells_in_both] <- new_annotations


p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
p1+p2

save(Myeloid, test, file = paste0(file_path, "/step1.2.relabel.cheadle.RData"))


