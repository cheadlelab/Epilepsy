rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(scales)


#################################################################################################################################
############################################################# snRNA #############################################################
#################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"
save_path <- "//grid/cheadle_home/qianyu/Epilepsy/figures"


load(paste0(file_path, "/step1.2.umap.fine.RData"))

unique(test@meta.data$Annotation.fine)

# Prepare the Seurat object
`%!in%` <- Negate(`%in%`)
test <- subset(test, subset = Annotation.fine %!in% c("B/cDC2", "B/Myeloid", "B/NKT", "B/other", "B/platelets", "Myeloid/NKT", "Myeloid/platelets"))

unique(test@meta.data$Annotation.fine)
unique(test@meta.data$patient)

patient_mapping <- c(
  "Epileptic_101923" = "patient1",
  "Epileptic_101323" = "patient2",
  "Epileptic_22924" = "patient3",
  "Epileptic_2849" = "patient4"
)

# Add the new 'patient_id' column to metadata
test@meta.data$patient_id <- ifelse(
  test@meta.data$patient %in% names(patient_mapping),
  patient_mapping[test@meta.data$patient],
  as.character(test@meta.data$patient)  # Keep original ID if not in mapping
)
test@meta.data$patient_id <- factor(
  test@meta.data$patient_id,
  levels = c(
    "patient1", "patient2", "patient3", "patient4",
    "Healthy_092524_400", "Healthy_092524_607",
    "Healthy_101124", "Healthy_102324", "Healthy_61023", 
    "Healthy_71324", "Healthy_91724"
  )
)


unique(test@meta.data$patient_id)

DefaultAssay(test) <- "RNA"
test@active.ident <- as.factor(test@meta.data$orig.ident)

unique(test@meta.data$Annotation.fine)

ordered_levels <- c(
  ## Myeloid lineage
  "CD14_Mono",      # Classical monocytes
  "CD16_Mono",      # Non-classical monocytes
  "cDC1",           # Conventional DC subset 1
  "cDC2",           # Conventional DC subset 2
  "pDC",            # Plasmacytoid DCs
  "Mast",           # Mast cells
  
  ## Lymphoid lineage
  ### B-cell series (progenitor → mature → effector)
  "pro_B",          
  "B_naive",        
  "B_intermediate", 
  "B_memory",       
  "Plasmablast",    
  
  ### T-cell series
  "CD4_Naive",      # Naïve CD4+
  "CD4_Memory",     # Memory CD4+
  "Treg",           # Regulatory T cells
  "CD4_CTL",        # CD4+ cytotoxic
  "CD8_Naive",      # Naïve CD8+
  "CD8_Memory",     # Memory CD8+
  "MAIT",           # Mucosal-associated invariant T
  "gdT",            # γδ T cells
  
  ### NK-cell series
  "NK_CD56bright",  # CD56bright NK
  "NK",             # CD56dim/adaptive NK
  
  ## Proliferating populations
  "MKI67+_NKT",     # Cycling NKT
  
  ## Erythroid & platelet lineages
  "EMP",            # Erythro-myeloid progenitors
  "Late_Eryth",     # Late erythroid cells
  "Platelets"       # Platelets
)



test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)
print(levels(test@meta.data$Annotation.fine))

default_colors <- hue_pal()(25)
names(default_colors) <- levels(test@meta.data$Annotation.fine)

custom_colors <- c("epileptic" = "#E64B35",
                   "healthy" = "#357EBD")

##################################################################################################
##################################### condition and celltype #####################################
##################################################################################################
pdf(file = paste0(save_path, "/Figure.blood.condition.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap", group.by = "condition", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = custom_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.blood.condition.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap", group.by = "condition", raster = FALSE) + scale_color_manual(values = custom_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.blood.annotation.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = default_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.blood.annotation.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap", group.by = "Annotation.fine", label = TRUE, raster = FALSE) + scale_color_manual(values = default_colors)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.blood.patient.nolegend.pdf"), width = 20, height = 20)
DimPlot(test, reduction = "umap", group.by = "patient_id", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.blood.patient.legend.pdf"), width = 30, height = 30)
DimPlot(test, reduction = "umap", group.by = "patient_id", label = TRUE, raster = FALSE)
dev.off() # Close the PDF device

###################################################################################################
############################################# patient #############################################
###################################################################################################

## here, pdf function doesn't work with loop
for (i in unique(test@meta.data$patient_id)) {
  
  tt <- subset(test, subset = patient_id == i)
  
  pdf(file = paste0(save_path, "/Figure.blood.", i, ".annotation.nolegend.pdf"), width = 10, height = 10)
  DimPlot(tt, group.by = "Annotation.fine", label = FALSE, raster = FALSE, pt.size = 0.5, shuffle = TRUE) + 
    NoLegend() + ggtitle(NULL) + NoAxes() + scale_color_manual(values = default_colors)
  dev.off() # Close the PDF device
  
  pdf(file = paste0(save_path, "/Figure.blood.", i, ".annotation.legend.pdf"), width = 30, height = 30)
  DimPlot(tt, group.by = "Annotation.fine", label = TRUE, raster = FALSE) + scale_color_manual(values = default_colors)
  dev.off() # Close the PDF device
  
}