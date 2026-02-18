rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
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

load("//grid/cheadle_home/qianyu/Epilepsy/output/step1.0.cheadle.merge.Elzar.RData")
Brain <- test
load("//grid/cheadle_home/qianyu/Epilepsy/output/Blood/step1.0.cheadle.Blood.RData")
Blood <- test


p1 <- VlnPlot(Brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal", "log10GenesPerUMI"), ncol = 5, pt.size = 0)
p2 <- VlnPlot(Blood, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal", "log10GenesPerUMI"), ncol = 5, pt.size = 0)

save_path <- "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V7"

pdf(file = paste0(save_path, "/Figure.snRNA.QC.pdf"), width = 5, height = 16)
VlnPlot(Brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), ncol = 1, pt.size = 0)
dev.off() # Close the PDF device

pdf(file = paste0(save_path, "/Figure.scRNA.QC.pdf"), width = 5, height = 16)
VlnPlot(Blood, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), ncol = 1, pt.size = 0)
dev.off() # Close the PDF device
