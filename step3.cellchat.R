# scvi envs
rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(patchwork)
library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)


file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output/Blood"
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
load(paste0(file_path, "/step1.2.umap.fine.RData"))
unique(test@meta.data$Annotation.fine)

options(future.globals.maxSize = 100 * 1024^3)  # Set limit to 10G

`%!in%` <- Negate(`%in%`)
test <- subset(test, subset = Annotation.fine %!in% c("Macrophage/platelets", "Myeloid/NKT", "RBC", "Platelets"))

levels_order <- c(
  "MKI67+_NKT", "NK_1", "NK_2", "NK_3",
  "B", "Plasma_B", "Pre_B",
  "Macrophage", "FCGR3A+ Mono", "DC", "cDCs", "pDC", "Mast",
  "Naive_CD4+_T", "Memory_CD4+_T", "CD8+_T_1", "CD8+_T_2", "CD8+_T_3", "regular_T", "gamma_delta_T"
)

test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = levels_order)
###############################################################################################################################
######################################################## Load Database ########################################################
###############################################################################################################################
# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

#############################################################################################################################
###################################################### loop by patient ######################################################
#############################################################################################################################

patient <- unique(test@meta.data$patient)
for (tmp in patient) {
  print(tmp)
  seurat_subset <- subset(test, subset = patient == tmp)

  Idents(seurat_subset) <- seurat_subset@meta.data$Annotation.fine
  # data.input <- GetAssayData(seurat_subset, assay = "RNA", layer = "layer", slot = "data") # normalized data matrix
  data.input <- LayerData(seurat_subset, assay = "RNA", layer = "data") # normalized data matrix
  
  labels <- Idents(seurat_subset)
  meta <- as.data.frame(seurat_subset@meta.data)
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Annotation.fine")
  cellchat@DB <- CellChatDB.use
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  # Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  
  options(future.availableCores.methods = c("mc.cores"))
  options(mc.cores = parallel::detectCores())
  future::plan("multisession", workers = 12) # do parallel
  
  #> processing ('multicore') is not supported when running R from RStudio
  #> because it is considered unstable. For more details, how to control forked
  #> processing or not, and how to silence this warning in future R sessions, see ?
  #> parallelly::supportsMulticore
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # cellchat <- projectData(cellchat, PPI.human)
  
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  # Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat)
  #> triMean is used for calculating the average gene expression per cell group. 
  #> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-1https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html1-26 08:34:38]"
  #> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 100)
  
  
  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  
  # manually save the figures
  # par(mfrow = c(1,2), xpd=TRUE)
  # netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  # netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  save(cellchat, file = paste0(file_path, "/step3.CellChat.", tmp, ".RData"))
}
