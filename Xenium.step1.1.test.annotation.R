# https://divingintogeneticsandgenomics.com/post/neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data-in-seurat-and-bioconductor/
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ
rm(list = ls(all = TRUE))

#load packages
library(Seurat)
library(dplyr)
library(jsonlite)
library(dplyr)
library(arrow)
library(spacexr)
library(ggplot2)
library(dplyr)

library(future)
library(future.apply)
plan(multicore, workers = 8)

options(future.globals.maxSize = 128 * 1024^3)

source("//grid/cheadle_home/qianyu/Epilepsy/Xenium/ReadXenium2.R")
# source("/grid/cheadle/home/qianyu/Epilepsy/Xenium/ReadXenium2.R")

save.path <- "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/"
# save.path <- "/grid/cheadle/home/qianyu/Epilepsy/Xenium/output/"
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/Xenium/20240726__223041__Cheadle_CKX01/"
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
exp.ident <- c("1 Non-epileptiform",
               "4 Epileptiform",
               "1 Non-epileptiform",
               "5 Non-epileptiform",
               "2 Epileptiform",
               "3 Epileptiform",
               "1 Epileptiform",
               "4 Non-epileptiform",
               "5 Epileptiform",
               "2 Non-epileptiform",
               "3 Epileptiform")

p <- list()

load("//grid/cheadle_home/qianyu/Epilepsy/output/step1.2.relabel.cheadle.RData")
table(test$Annotation.fine)
test <- JoinLayers(test)
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
test <- subset(test, subset = Annotation.fine != "Other")
DefaultAssay(test) <- "RNA"


p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.coarse", raster=FALSE, label = TRUE)
p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", raster=FALSE, label = TRUE)
p1+p2

json_file <- "//grid/cheadle_home/qianyu/Epilepsy/Xenium/7G69XV_gene_panel-permissive20240802-shlee.json"
cheadle.epilepsy.ref <- test
test <- NULL
gc()

# for (i in c(1, 2, 5, 7, 8, 10)) {
for (i in 1:11) {
  
  n.cores <- 12
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
  
  permissive <- fromJSON(json_file)
  a <- permissive[["payload"]][["targets"]][["type"]]
  b <- as.matrix(data[["matrix"]][["Gene Expression"]])
  legal <- a$data$name[which(a$descriptor=="gene")]
  pool <- rownames(data[["matrix"]][["Gene Expression"]])
  useful <- intersect(legal, pool)
  data$matrix$RNA <- data[["matrix"]][["Gene Expression"]][useful,]
  
  options(Seurat.object.assay.version = "v3")
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["RNA"]],assay = "RNA")
  #xenium.obj[["RNA"]] <- as(object = xenium.obj[["RNA"]], Class = "Assay")
  
  
  xenium.obj[["GeneExpression"]] <- CreateAssayObject(counts = data$matrix[["Gene Expression"]])
  xenium.obj[["DeprecatedCodeword"]] <- CreateAssayObject(counts = data$matrix[["Deprecated Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["NegativeControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[["UnassignedCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["fov"]] <- coords
  
  xenium.obj <- subset(xenium.obj, subset = nCount_RNA > 0)
  xenium.obj@meta.data$orig.ident <- samples[i]
  
  
  
  ######################################################################################
  ##################################### Annotation ##################################### 
  ######################################################################################
  query.counts <- GetAssayData(xenium.obj, assay = "RNA", slot = "counts")[, Cells(xenium.obj[["RNA"]])]
  coords <- GetTissueCoordinates(xenium.obj[["fov"]], which = "centroids")
  
  rownames(coords) <- coords$cell
  coords$cell <- NULL
  query <- SpatialRNA(coords, query.counts, colSums(query.counts))
  
  # process_cell_type_info error: need a minimum of 25 cells for each cell type in the reference
  # some cell types don't fit the requirement
  # cheadle.epilepsy.ref <- subset(cheadle.ref, subset = tmp == exp.ident[i])
  # cheadle.epilepsy.ref <- UpdateSeuratObject(cheadle.epilepsy.ref)
  
  Idents(cheadle.epilepsy.ref) <- "Annotation.fine"
  # remove CR cells because there aren't enough of them for annotation
  # cheadle.epilepsy.ref <- subset(cheadle.epilepsy.ref, subset = Annotation != "CR")
  counts <- GetAssayData(cheadle.epilepsy.ref, assay = "RNA", slot = "counts")
  cluster <- as.factor(cheadle.epilepsy.ref$Annotation.fine)
  names(cluster) <- colnames(cheadle.epilepsy.ref)
  nUMI <- cheadle.epilepsy.ref$nCount_RNA
  names(nUMI) <- colnames(cheadle.epilepsy.ref)
  nUMI <- colSums(counts)
  levels(cluster) <- gsub("/", "-", levels(cluster))
  reference <- Reference(counts, cluster, nUMI)
  # run RCTD with many cores
  RCTD <- create.RCTD(query, reference, max_cores = n.cores)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  annotations.df <- RCTD@results$results_df
  annotations <- annotations.df$first_type
  names(annotations) <- rownames(annotations.df)
  xenium.obj$predicted.celltype <- annotations
  keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
  xenium.obj <- subset(xenium.obj, cells = keep.cells)
  
  xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "predicted.celltype",
                                niches.k = 5, neighbors.k = 30)
  
  # celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
  # niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") + scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
  # celltype.plot | niche.plot
  
  #######################################################################################
  DefaultAssay(xenium.obj) <- "RNA"
  assign(paste0("sample_",i), xenium.obj)
  
  save(list = paste0("sample_", i), file = paste0(save.path, "step1.sample_", i, ".RData"))
  
  xenium.obj <- NULL
  gc()
  # p[i] <- xenium.obj
}

#######################################################################################
######################################## Layer ########################################
#######################################################################################


Inh <- c("Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier")
ExcL23 <- c("Exc_L234IT_ERBB4",  # L2-3-4
            "Exc_L23IT"        # L2/3
)
ExcL456 <- c("Exc_L4IT",          # L4
             "Exc_L4IT_PLCH1",    # L4
             "Exc_L5ET",          # L5
             "Exc_L5IT_GRIN3A",    # L5
             "Exc_L5IT_GRIN3A-",   # L5
             "Exc_L56NP",         # L5-6
             "Exc_L56IT_CAR3",    # L5-6
             "Exc_L6IT",          # L6
             "Exc_L6CT",          # L6
             "Exc_L6b")           # L6
glia <- c("Astrocyte", "OPC", "Olg", "Microglia", "Vascular")
other <- c( "Macrophage", "NK/T")

for (i in 1:11) {
  load(paste0(save.path, "step1.sample_", i, ".RData"))
  sample_name <- paste0("sample_", i)
  sample_obj <- get(sample_name)
  
  sample_obj@meta.data <- sample_obj@meta.data %>%
    mutate(celltype.corse = case_when(
      predicted.celltype %in% Inh ~ "Inh",
      predicted.celltype %in% ExcL23 ~ "Exc_L23",
      predicted.celltype %in% ExcL456 ~ "Exc_L456",
      predicted.celltype %in% glia ~ predicted.celltype,
      predicted.celltype %in% other ~ predicted.celltype,
      TRUE ~ "Other"
    ))
  
  assign(sample_name, sample_obj)
}
## rename the niche
# BuildNicheAssay seems like will use random seed, try several times until get proper results
#######################################################################################
# sample_1
xenium.obj <- sample_1

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 3, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-NVU",
  "2" = "Layer-Exc23",
  "3" = "Layer-NVU"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_1 <- xenium.obj
#######################################################################################
# sample_2
xenium.obj <- sample_2

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 4, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-Exc456",
  "2" = "Layer-Exc23",
  "3" = "Layer-NVU",
  "4" = "Layer-Olg"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_2 <- xenium.obj
#######################################################################################
# sample_3
xenium.obj <- sample_3

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 4, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-Olg",
  "2" = "Layer-Olg",
  "3" = "Layer-Olg",
  "4" = "Layer-Olg"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_3 <- xenium.obj
#######################################################################################
# sample_4
xenium.obj <- sample_4

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 1, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Other"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_4 <- xenium.obj
#######################################################################################
# sample_5
xenium.obj <- sample_5

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 4, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-Exc23",
  "2" = "Layer-NVU",
  "3" = "Layer-Exc456",
  "4" = "Layer-Olg"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_5 <- xenium.obj
#######################################################################################
# sample_6
xenium.obj <- sample_6

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 5, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-NVU",
  "2" = "Layer-Exc23",
  "3" = "Layer-NVU",
  "4" = "Layer-Olg",
  "5" = "Layer-Exc456"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_6 <- xenium.obj
#######################################################################################
# sample_7
xenium.obj <- sample_7

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 5, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-NVU",
  "2" = "Layer-Olg",
  "3" = "Layer-NVU",
  "4" = "Layer-Exc456",
  "5" = "Layer-Exc23"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_7 <- xenium.obj
#######################################################################################
# sample_8
xenium.obj <- sample_8

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 4, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-Olg",
  "2" = "Layer-Exc23",
  "3" = "Layer-NVU",
  "4" = "Layer-Exc456"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_8 <- xenium.obj
#######################################################################################
# sample_9
xenium.obj <- sample_9

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 3, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-Olg",
  "2" = "Layer-Olg",
  "3" = "Layer-Olg"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_9 <- xenium.obj
#######################################################################################
# sample_10
xenium.obj <- sample_10

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 4, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-NVU",
  "2" = "Layer-Exc456",
  "3" = "Layer-Exc23",
  "4" = "Layer-NVU"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_10 <- xenium.obj
#######################################################################################
# sample_11
xenium.obj <- sample_11

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "celltype.corse",
                              niches.k = 3, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "celltype.corse", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches")
celltype.plot | niche.plot
table(xenium.obj@meta.data$celltype.corse, xenium.obj@meta.data$niches)
layer_mapping <- c(
  "1" = "Layer-NVU",
  "2" = "Layer-NVU",
  "3" = "Layer-NVU",
  "4" = "Layer-NVU"
)
xenium.obj@meta.data$layer <- layer_mapping[xenium.obj@meta.data$niches]
layer.plot <- ImageDimPlot(xenium.obj, group.by = "layer", size = 1.5, dark.background = F) + ggtitle("Layers")
celltype.plot | layer.plot
sample_11 <- xenium.obj


#######################################################################################
######################################## Merge ########################################
#######################################################################################

sample_list <- list(
  "sample1" = sample_1,
  "sample2" = sample_2,
  "sample3" = sample_3,
  "sample4" = sample_4,
  "sample5" = sample_5,
  "sample6" = sample_6,
  "sample7" = sample_7,
  "sample8" = sample_8,
  "sample9" = sample_9,
  "sample10" = sample_10,
  "sample11" = sample_11
)
save(sample_list, file = paste0(save.path, "step1.list.RData"))


test <- merge(sample_1, 
              y = c(
                sample_2, sample_3, sample_4, sample_5, sample_6,
                sample_7, sample_8, sample_9, sample_10, sample_11
              ),
              project = "Xenium",
              merge.data = TRUE)
save(test, file = paste0(save.path, "step1.merged.RData"))





#######################################################################################
######################################## Image ########################################
#######################################################################################
library(scales)

identities <- unique(cheadle.epilepsy.ref@meta.data$Annotation.fine)
my_color_palette <- hue_pal()(length(identities))

xenium.obj <- sample_10
celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome", dark.background = F) + ggtitle("Cell type") +
  scale_fill_manual(values = c("#442288", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", 
                               "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", 
                               "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", 
                               "#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3", 
                               "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"
  ))
# scale_fill_manual(values = c("#F8766D", "#EC823C", "#DD8D00", "#CA9700", "#B3A000", "#97A900", "#71B000", "#2FB600", "#00BB4B", "#00BF76", "#00C098", "#00C0B7", "#00BDD1", "#00B7E8", "#00AEFA", "#3DA1FF", "#8F91FF", "#BE80FF", "#DE71F9", "#F265E7", "#FE61CF", "#FF64B3", "#FF6C92"))
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") + scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot
