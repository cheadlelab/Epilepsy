# conda activate condaR

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
library(future)

##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
os_name <- Sys.info()[["sysname"]]   # "Linux", "Windows", "Darwin" …
file_path <- switch(
  os_name,
  "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output", 
  "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",
  stop(sprintf("Unsupported OS: %s", os_name))
)
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output"
options(future.globals.maxSize = 500 * 1024^3)  # Set limit to 50G
options(future.availableCores.methods = c("mc.cores"))
options(mc.cores = parallel::detectCores())
# future::plan("multisession", workers = 8) # rstudio
future::plan("multicore", workers = 32) # r

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
load(paste0(file_path, "/step3.1.cellchat.prepare.RData"))
table(test$Annotation.fine, test$condition)

# for (tmp in c("Focal", "Distal", "Stimulated")) {
for (tmp in c("Distal")) {
  print(tmp)
  test@meta.data$Annotation.new <- as.character(test@meta.data$Annotation.fine)
  custom_order <- c(
    "Exc_L23IT", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", 
    "Exc_L5IT_GRIN3A", "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", 
    "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", 
    "Exc_L234IT_ERBB4",
    "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", 
    "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier",
    "Other Astrocyte", "Reactive Astrocyte", "Lipid-Accumulated Reactive Astrocyte",
    "Microglia", "OL", "OxPos OL", "OPC", 
    "Macrophage", "NK/T", "Vascular"
  )
  test@meta.data$Annotation.new <- factor(test@meta.data$Annotation.new, levels = custom_order)
  unique(test@meta.data$Annotation.new)
  seurat_subset <- subset(test, subset = condition == tmp)
  gc()
  Idents(seurat_subset) <- seurat_subset@meta.data$Annotation.new
  data.input <- LayerData(seurat_subset, assay = "RNA", layer = "data") # normalized data matrix
  labels <- Idents(seurat_subset)
  meta <- as.data.frame(seurat_subset@meta.data)
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Annotation.new")
  cellchat@idents <- droplevels(cellchat@idents)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 20)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  if (tmp == "Focal") {
    Focal <- cellchat
    save(Focal, file = paste0(file_path, "/step3.CellChat.Focal.RData"))
  } else if (tmp == "Distal") {
    Distal <- cellchat
    save(Distal, file = paste0(file_path, "/step3.CellChat.Distal.RData"))
  } else if (tmp == "Stimulated") {
    Stimulated <- cellchat
    save(Stimulated, file = paste0(file_path, "/step3.CellChat.Stimulated.RData"))
  }
  seurat_subset <- NULL
  gc()
}

