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

os_name <- Sys.info()[["sysname"]]   # "Linux", "Windows", "Darwin" …
save.path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/qianyu/Epilepsy/Xenium/output/",
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/",
                    stop(sprintf("Unsupported OS: %s", os_name))
)
for (i in 1:11) {
  load(paste0(save.path, "step1.sample_", i, ".RData"))
}


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

filenames <- c("XETG00234__0036747__Region_1",
               "XETG00234__0036747__Region_2",
               "XETG00234__0036747__Region_3",
               "XETG00234__0036747__Region_4",
               "XETG00234__0036747__Region_5",
               "XETG00234__0036747__Region_6",
               "XETG00234__0036760__Region_1",
               "XETG00234__0036760__Region_2",
               "XETG00234__0036760__Region_3",
               "XETG00234__0036760__Region_4",
               "XETG00234__0036760__Region_5"
)


test <- merge(sample_1, 
              y = c(
                sample_2, sample_3, sample_4, sample_5, sample_6,
                sample_7, sample_8, sample_9, sample_10, sample_11
              ),
              add.cell.ids = filenames,
              project = "Xenium",
              merge.data = TRUE)


DefaultAssay(test) <- "RNA"
test <- NormalizeData(test)
test <- FindVariableFeatures(test)
test <- ScaleData(test)
test <- RunPCA(test)

dims_parameter <- 25
test <- RunUMAP(test, reduction = "pca", dims = 1:dims_parameter, reduction.key = "normal_")
test <- FindNeighbors(test, reduction = "pca", dims = 1:dims_parameter)
test <- FindClusters(test, resolution = 0.2)
test@reductions$umap_normal <- test@reductions$umap
test@meta.data$clusters_normal <- test@meta.data$seurat_clusters

p1 <- DimPlot(test, reduction = "umap_normal", group.by = "clusters_normal", raster=FALSE, label = TRUE)+NoLegend()
p2 <- DimPlot(test, reduction = "umap_normal", group.by = "predicted.celltype", raster=FALSE, label = TRUE)+NoLegend()
p3 <- DimPlot(test, reduction = "umap_normal", group.by = "orig.ident", raster=FALSE, label = TRUE)+NoLegend()
p2

save(test, file = paste0(save.path, "step2.umap.RData"))
load(paste0(save.path, "step2.umap.RData"))


##################################################################################################
##################################### condition and celltype #####################################
##################################################################################################
test@meta.data$Annotation.fine <- test@meta.data$predicted.celltype
table(test@meta.data[["Annotation.fine"]])

orig.ident <- c("1C SA", "4E", "1C SB", "5C", "2E", "3E SB", "1E", "4C", "5E", "2C", "3E SA")
patient <- c(1, 4 ,1, 5, 2, 3, 1, 4, 5, 2, 3)
condition <- c("Distal", "Focal", "Distal", "Distal", "Focal", "Focal", "Focal", "Distal", "Focal", "Distal", "Focal")

test@meta.data$patient <- patient[match(test@meta.data$orig.ident, orig.ident)]
test@meta.data$condition <- condition[match(test@meta.data$orig.ident, orig.ident)]

GABAergic <- c("Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier")
Glu <- c("Exc_L234IT_ERBB4",  # L2-3-4
         "Exc_L23IT",         # L2/3
         "Exc_L4IT",          # L4
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
other <- c( "Macrophage", "NK-T")

ordered_levels <- c(GABAergic, Glu, glia, other)
test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)
print(levels(test@meta.data$Annotation.fine))


default_colors <- hue_pal()(26)
names(default_colors) <- levels(test@meta.data$Annotation.fine)

custom_colors <- c(
  "Focal" = "#d62728",
  "Distal" = "#2ca02c",
  "Stimulated" = "#1f77b4"
)


#################################################################################################
############################################# image #############################################
#################################################################################################

pdf(file = paste0(save.path, "/Figure1.Xenium.condition.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap_normal", group.by = "condition", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = custom_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save.path, "/Figure1.Xenium.condition.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap_normal", group.by = "condition", raster = FALSE) + scale_color_manual(values = custom_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save.path, "/Figure1.Xenium.annotation.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap_normal", group.by = "Annotation.fine", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = default_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save.path, "/Figure1.Xenium.annotation.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap_normal", group.by = "Annotation.fine", label = TRUE, raster = FALSE) + scale_color_manual(values = default_colors)
dev.off() # Close the PDF device

























Idents(test) <- "clusters_normal"
test.markers1 <- FindAllMarkers(test, only.pos = FALSE, min.pct = 0.25, logfc.threshold = log2(1.5))

xenium.obj <- test

save(test, file = paste0(save.path, "step2.umap.RData"))

# load("/home/qianyu/Desktop/grid/Epilepsy/output/step2.Merge.relabel.RData")
# FeaturePlot(test, features = c("CHOLD", "IDO1"))
FeaturePlot(test, features = c("PTPRC"))



Idents(test) <- "condition"
markers_list <- list()

for (i in unique(test@meta.data$predicted.celltype)) {
  tmp <- subset(test, subset = predicted.celltype == i)
  markers_list[[i]] <- FindAllMarkers(tmp, only.pos = FALSE, min.pct = 0.25, logfc.threshold = log2(1.5))
}

