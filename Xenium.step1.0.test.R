# https://divingintogeneticsandgenomics.com/post/neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data-in-seurat-and-bioconductor/
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager

rm(list = ls(all = TRUE))

#load packages
library(Seurat)
library(dplyr)
library(jsonlite)
library(dplyr)
library(arrow)

options(future.globals.maxSize = 512 * 1024^3)

source("/home/qianyu/Desktop/grid/Epilepsy/Xenium/ReadXenium2.R")
# source("/grid/cheadle/home/qianyu/Epilepsy/Xenium/ReadXenium2.R")

save.path <- "/home/qianyu/Desktop/grid/Epilepsy/Xenium/output/"
# save.path <- "/grid/cheadle/home/qianyu/Epilepsy/Xenium/output/"
file_path <- "/home/qianyu/Desktop/grid/Epilepsy/Xenium/20240726__223041__Cheadle_CKX01/"
# file_path <- "/grid/cheadle/home/qianyu/Epilepsy/Xenium/20240726__223041__Cheadle_CKX01/"
filenames <- c("output-XETG00234__0036747__Region_1__20240726__223131",
               "output-XETG00234__0036747__Region_2__20240726__223131",
               "output-XETG00234__0036747__Region_3__20240726__223131",
               "output-XETG00234__0036747__Region_4__20240726__223131",
               "output-XETG00234__0036747__Region_5__20240726__223131",
               "output-XETG00234__0036747__Region_6__20240726__223131",
               "output-XETG00234__0036760__Region_1__20240726__223131",
               "output-XETG00234__0036760__Region_2__20240726__223131",
               "output-XETG00234__0036760__Region_3__20240726__223131",
               "output-XETG00234__0036760__Region_4__20240726__223131",
               "output-XETG00234__0036760__Region_5__20240726__223131"
)
samples <- c("1C SA",
             "4E",
             "1C SB",
             "5C",
             "2E",
             "3E SB",
             "1E",
             "4C",
             "5E",
             "2C",
             "3E SA")
p <- list()

for (i in 1:11) {
  
  filename <- filenames[i]
  data <- ReadXenium2(
    data.dir = paste0(file_path, filename),
    type = c("centroids", "segmentations"),
  )
  
  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = "Xenium"
  )
  
  
  json_file <- "/home/qianyu/Desktop/grid/Epilepsy/Xenium/7G69XV_gene_panel-permissive20240802-shlee.json"
  permissive <- fromJSON(json_file)
  a <- permissive[["payload"]][["targets"]][["type"]]
  b <- as.matrix(data[["matrix"]][["Gene Expression"]])
  legal <- a$data$name[which(a$descriptor=="gene")]
  pool <- rownames(data[["matrix"]][["Gene Expression"]])
  useful <- intersect(legal, pool)
  data$matrix$RNA <- data[["matrix"]][["Gene Expression"]][useful,]
  
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["RNA"]], assay = "RNA")
  xenium.obj[["GeneExpression"]] <- CreateAssayObject(counts = data$matrix[["Gene Expression"]])
  xenium.obj[["DeprecatedCodeword"]] <- CreateAssayObject(counts = data$matrix[["Deprecated Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["NegativeControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[["UnassignedCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["fov"]] <- coords
  
  xenium.obj <- subset(xenium.obj, subset = nCount_RNA > 0)
  xenium.obj@meta.data$orig.ident <- samples[i]
  
  assign(paste0("sample_",i), xenium.obj)
  
  # p[i] <- xenium.obj
}

# i = 11
# ImageDimPlot(p[[i]], cols = "polychrome", size = 0.75)

test <- merge(sample_1, 
              y = c(
                sample_2, sample_3, sample_4, sample_5, sample_6,
                sample_7, sample_8, sample_9, sample_10, sample_11
              ),
              add.cell.ids = c(
                "XETG00234_R1", "XETG00234_R2", "XETG00234_R3", "XETG00234_R4", "XETG00234_R5", "XETG00234_R6", 
                "XETG00234_R1", "XETG00234_R2", "XETG00234_R3", "XETG00234_R4", "XETG00234_R5"
              ),
              project = "Xenium",
              merge.data = TRUE)

VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, group.by = "orig.ident")
save(test, file = paste0(save.path, "step0.XeniumObject.RData"))


test <- merge(sample_1, 
              y = c(
                sample_2, sample_5,
                sample_7, sample_8, sample_10
              ),
              add.cell.ids = c(
                "XETG00234_R1", "XETG00234_R2", "XETG00234_R5", 
                "XETG00234_R1", "XETG00234_R2", "XETG00234_R4"
              ),
              project = "Xenium",
              merge.data = TRUE)

VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0, group.by = "orig.ident")
save(test, file = paste0(save.path, "step0.XeniumObject.QC.RData"))




test <- NormalizeData(test)
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test)

test <- SCTransform(test, assay = "RNA")
test <- RunPCA(test, npcs = 30, features = rownames(test))
test <- RunUMAP(test, dims = 1:30)
test <- FindNeighbors(test, reduction = "pca", dims = 1:30)
test <- FindClusters(test, resolution = 0.3)


test <- IntegrateLayers(
  object = test, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


pdf(file = paste0(save.path, "/elbow.pdf"), width = 10, height = 10)
ElbowPlot(test)
dev.off() # Close the PDF device

dims_parameter <- 30
test <- RunUMAP(test, reduction = "harmony", dims = 1:dims_parameter, reduction.key = "harmony_")
test <- FindNeighbors(test, reduction = "harmony", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.2)
test@reductions$umap_harmony <- test@reductions$umap
test@meta.data$clusters_harmony <- test@meta.data$seurat_clusters

dims_parameter <- 30
test <- RunUMAP(test, reduction = "pca", dims = 1:dims_parameter, reduction.key = "normal_")
test <- FindNeighbors(test, reduction = "pca", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.2)
test@reductions$umap_normal <- test@reductions$umap
test@meta.data$clusters_normal <- test@meta.data$seurat_clusters

pdf(file = paste0(save.path, "/umap.pdf"), width = 20, height = 10)
p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "clusters_harmony", raster=FALSE, label = TRUE)+NoLegend()
p2 <- DimPlot(test, reduction = "umap_normal", group.by = "clusters_normal", raster=FALSE, label = TRUE)+NoLegend()
p1+p2
dev.off() # Close the PDF device

save(test, file = paste0(save.path, "step1.test.QC.RData"))





