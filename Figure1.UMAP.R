rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(scales)


#################################################################################################################################
############################################################# snRNA #############################################################
#################################################################################################################################
os_name <- Sys.info()[["sysname"]]   # returns "Linux", "Windows", "Darwin" (macOS), etc.

file_path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output",   # Ubuntu and other Linux
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",  # Windows UNC path
                    # fallback – change or add more cases if you need them
                    stop(sprintf("Unsupported OS: %s", os_name))
)

save_path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/figures",   # Ubuntu and other Linux
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/figures",  # Windows UNC path
                    # fallback – change or add more cases if you need them
                    stop(sprintf("Unsupported OS: %s", os_name))
)

# load(paste0(file_path, "/step1.4.final.cheadle.seurat4.RData"))
load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))

# Prepare the Seurat object
test <- JoinLayers(test)
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
test <- subset(test, subset = Annotation.fine != "Other")
DefaultAssay(test) <- "RNA"

test@active.ident <- as.factor(test@meta.data$orig.ident)

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
other <- c( "Macrophage", "NK/T")

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


##################################################################################################
##################################### condition and celltype #####################################
##################################################################################################
pdf(file = paste0(save_path, "/Figure1.snRNA.condition.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap_harmony", group.by = "condition", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = custom_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.snRNA.condition.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap_harmony", group.by = "condition", raster = FALSE) + scale_color_manual(values = custom_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.snRNA.annotation.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = default_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure1.snRNA.annotation.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap_harmony", group.by = "Annotation.fine", label = TRUE, raster = FALSE) + scale_color_manual(values = default_colors)
dev.off() # Close the PDF device


###################################################################################################
############################################# patient #############################################
###################################################################################################

results_list <- list()

for (i in 1:6) {
  tt <- subset(test, subset = patient == i)
  tt[["RNA"]] <- as(object = tt[["RNA"]], Class = "Assay5")
  tt[["RNA"]] <- split(tt[["RNA"]], f = tt$orig.ident)
  
  tt <- NormalizeData(tt)
  tt <- FindVariableFeatures(tt)
  tt <- ScaleData(tt)
  tt <- RunPCA(tt)
  
  if (length(unique(tt@meta.data$orig.ident))>1) {
    tt <- IntegrateLayers(
      object = tt, method = HarmonyIntegration,
      orig.reduction = "pca", new.reduction = "harmony",
      verbose = FALSE
    )
    
    ElbowPlot(tt)
    dims_parameter <- 15
    tt <- RunUMAP(tt, reduction = "harmony", dims = 1:dims_parameter)
  } else {
    ElbowPlot(tt)
    dims_parameter <- 15
    tt <- RunUMAP(tt, reduction = "pca", dims = 1:dims_parameter)
  }
  
  results_list[[paste0("patient_", i)]] <- tt
  
}


custom_colors <- c(
  "Focal" = "#1f77b4",
  "Distal" = "#ff7f0e",
  "Stimulated" = "#2ca02c"
)
custom_colors <- c(
  "Focal" = "#d62728",
  "Distal" = "#2ca02c",
  "Stimulated" = "#1f77b4"
)

default_colors <- hue_pal()(26)
names(default_colors) <- levels(test@meta.data$Annotation.fine)

## here, pdf function doesn't work with loop
for (i in 1:6) {
  tt <- results_list[[paste0("patient_", i)]]
  
  pdf(file = paste0(save_path, "/Figure1.patient", i, ".condition.nolegend.pdf"), width = 10, height = 10)
  DimPlot(tt, group.by = "condition", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
    NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = custom_colors)
  dev.off() # Close the PDF device
  
  pdf(file = paste0(save_path, "/Figure1.patient", i, ".condition.legend.pdf"), width = 30, height = 30)
  DimPlot(tt, group.by = "condition", raster = FALSE) + scale_color_manual(values = custom_colors)
  dev.off() # Close the PDF device
  
  pdf(file = paste0(save_path, "/Figure1.patient", i, ".annotation.nolegend.pdf"), width = 10, height = 10)
  DimPlot(tt, group.by = "Annotation.fine", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
    NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = default_colors)
  dev.off() # Close the PDF device
  
  pdf(file = paste0(save_path, "/Figure1.patient", i, ".annotation.legend.pdf"), width = 30, height = 30)
  DimPlot(tt, group.by = "Annotation.fine", label = TRUE, raster = FALSE) + scale_color_manual(values = default_colors)
  dev.off() # Close the PDF device
  
}