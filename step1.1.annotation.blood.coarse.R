# scvi envs
rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)

##################################################################################################################################
###################################################### load Cheadle Dataset ######################################################
##################################################################################################################################
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output/Blood"

# file_path <- "/grid/Epilepsy/output"

##################################################################################################################################
####################################################### load Blood Dataset #######################################################
##################################################################################################################################

load(paste0(file_path, "/step1.0.cheadle.Blood.RData"))


test <- JoinLayers(test)
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")


test <- subset(
  test, 
  subset = patient %in% c(
    "Epileptic_101323", "Epileptic_101923", "Epileptic_22924", 
    "Epileptic_2849", "Healthy_092524_400", "Healthy_092524_607",
    "Healthy_101124", "Healthy_102324", "Healthy_61023", 
    "Healthy_71324", "Healthy_91724"
  )
)

test@meta.data$condition <- ifelse(
  grepl("Epileptic", test@meta.data$patient), "epileptic", "healthy"
)

# test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay5")
# test[["RNA"]] <- split(test[["RNA"]], f = test@meta.data$orig.ident)

table(test@meta.data$patient)

test <- NormalizeData(test)
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test)

# test <- JoinLayers(test)
pdf(file = paste0(file_path, "/blood.elbow.pdf"), width = 10, height = 10)
ElbowPlot(test)
dev.off() # Close the PDF device


#####################################################################################################################################
############################################################## Label 1 ##############################################################
#####################################################################################################################################

dims_parameter <- 25
test <- RunUMAP(test, reduction = "pca", dims = 1:dims_parameter)
test <- FindNeighbors(test, reduction = "pca", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.1)


# plot cluster_label_simplify
pdf(file = paste0(file_path, "/blood.25.pdf"), width = 45, height = 15)
p1 <- DimPlot(test, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap", group.by = "orig.ident", raster=FALSE, label = TRUE)
p3 <- DimPlot(test, reduction = "umap", group.by = "condition", raster=FALSE, label = TRUE)
p1+p2+p3
dev.off() # Close the PDF device


save(test, file = paste0(file_path, "/step1.1.umap.Elzar.tmp.RData"))
# turn to Moss
load(paste0(file_path, "/step1.1.umap.Elzar.tmp.RData"))


Idents(test) = "seurat_clusters"
a <- FindMarkers(test, ident.1 = 7, ident.2 = NULL, min.pct = 0.1, only.pos = TRUE)
a <- FindAllMarkers(test, min.pct = 0.1, only.pos = TRUE)

# Combined FeaturePlot for PBMC cell type markers

pdf(file = paste0(file_path, "/blood.tmp.pdf"), width = 40, height = 90)
FeaturePlot(test, features = c(
  "PPBP", # Platelets
  "CD3D", "CD2", # T cells
  "CD4", "IL7R", # CD4+ T cells
  "CD8A", "CD8B", # CD8+ T cells
  "FOXP3", "IL2RA", # Regulatory T cells
  "CCR7", "SELL", # Naive T cells
  "GZMB", "PRF1", # Effector T cells
  "MS4A1", "CD79A", "CD79B", # B cells
  "MZB1", "XBP1", "PRDM1", # Plasma cells
  "NKG7", "GNLY", "FCGR3A", # NK cells
  "CD14", "LYZ", "MSR1", # Monocytes
  "FCN1", # Classical monocytes
  "FCGR3A", "CD16", # Non-classical monocytes
  "ITGAX", "HLA-DRA", # Dendritic cells
  "IL3RA", "LILRA4", # Plasmacytoid dendritic cells
  # "CLC", "IL5RA", # Eosinophils
  "HBA1", "HBB", # RBCs
  "ENPP3", "FCER1A", # Basophils
  "S100A8", "S100A9" # Neutrophils
))
dev.off() # Close the PDF device




Idents(test) <- "seurat_clusters"
# https://v22.proteinatlas.org/ENSG00000151789-ZNF385D/single+cell+type/pbmc

new_identities <- c("0" = "Myeloid",
                    "1" = "NKT",
                    "2" = "NKT",
                    "3" = "B",
                    "4" = "NKT",
                    "5" = "Myeloid",
                    "6" = "Myeloid",
                    "7" = "Other",
                    "8" = "pDC",
                    "9" = "B",
                    "10" = "NKT"
)

test@meta.data$Annotation.coarse <- new_identities[test@meta.data$seurat_clusters]

pdf(file = paste0(file_path, "/blood.coarse.pdf"), width = 15, height = 15)
DimPlot(test, reduction = "umap", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
dev.off() # Close the PDF device

save(test, file = paste0(file_path, "/step1.1.umap.coarse.RData"))

load(paste0(file_path, "/step1.1.umap.coarse.RData"))
