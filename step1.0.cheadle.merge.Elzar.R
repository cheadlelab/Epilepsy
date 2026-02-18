# module load EBModules
# module load Anaconda3/2021.05-SingleCell
# conda create --name velocyto python=3.7
# conda activate ElzarR
# R

rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(velocyto.R)
library(scDblFinder)

###############################################################################################################################
###################################################### Path and Metadata ######################################################
###############################################################################################################################
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/Cellranger_output"
save_path <- "/grid/cheadle/home/qianyu/Epilepsy/output"

datalist <- c(
  "Exp01",
  "Exp02_Control",
  "Exp02_ES",
  "Exp02_EU",
  "Exp03_EP",
  "Exp03_HU",
  "Exp04_3E",
  "Exp04_4C",
  "Exp04_4E",
  "Exp04_4S",
  "Exp05_5C",
  "Exp05_5E",
  "Exp05_5S",
  "Exp06_6C",
  "Exp06_6E",
  "Exp06_6S",
  "Exp07_4C",
  "Exp07_4E",
  "Exp07_5S"
)

###############################################################################################################################
######################################################### Load and QC #########################################################
###############################################################################################################################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_objects <- list()
vlnplots_list <- list()

for (sample_id in datalist) {
  
  print(paste0("Processing ", sample_id))
  # Construct file path for .h5 file
  data.file <- paste0(file_path, "/", sample_id, "/outs/output_filtered_seurat.h5")
  
  # Check if file exists before attempting to read
  if(file.exists(data.file)) {
    data.data <- Read10X_h5(filename = data.file, use.names = TRUE)
    obj <- CreateSeuratObject(counts = data.data, project = sample_id)
    obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
    
    ## Normalization and cell cycle score
    print("Normalization")
    obj <- NormalizeData(obj)
    obj <- CellCycleScoring(obj, 
                            g2m.features = g2m.genes, 
                            s.features = s.genes)
    
    ## QC rough
    print("Remove MT over expressed cells")
    obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
    obj <- PercentageFeatureSet(obj, pattern = "^RP[SL]", col.name = "percent.ribosomal")
    
    obj <- subset(obj, subset = nFeature_RNA > 200)
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- ScaleData(obj)
    
    ## Doublet identification
    print("Doublet identification")
    tmp <- as.SingleCellExperiment(obj)
    tmp <- scDblFinder(tmp)
    tmp <- as.Seurat(tmp)
    obj@meta.data <- tmp@meta.data
    rm(tmp)
    
    obj <- subset(obj, subset = scDblFinder.class == "singlet")
    pdf(file = paste0(save_path, "/", sample_id, ".pdf"), width = 20, height = 4)
    p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal", "log10GenesPerUMI"), ncol = 5)
    print(p1)
    dev.off() # Close the PDF device
    #vlnplots_list[[sample_id]] <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    ## QC fine
    obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 2.5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.8)
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- ScaleData(obj)
    assign(sample_id, obj)
    # save(obj, file = paste0(paste0(save_path, "/", sample_id, ".RData")))
    rm(obj)
  } else {
    cat("File does not exist:", data.file, "\n")
  }
}




test <- merge(Exp01, 
              y = c(Exp02_Control, Exp02_ES, Exp02_EU, Exp03_EP, Exp03_HU, Exp04_3E, Exp04_4C, Exp04_4E, Exp04_4S, 
                    Exp05_5C, Exp05_5E, Exp05_5S, Exp06_6C, Exp06_6E, Exp06_6S, Exp07_4C, Exp07_4E, Exp07_5E),
              project = "Epilepsy",
              merge.data = TRUE)

layersList <- lapply(test@assays$RNA@layers, function(x) {dim(x)})
layers_to_null <- names(layersList)[!grepl("^counts", names(layersList))]
test@assays$RNA@layers[layers_to_null] <- NULL

save(test, file = paste0(save_path, "/step1.0.cheadle.merge.Elzar.RData"))