# module load EBModules
# module load Anaconda3/2021.05-SingleCell

# conda activate ElzarR
# R

rm(list = ls(all = TRUE))
gc()
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)
library(velocyto.R)
###############################################################################################################################
###################################################### Path and Metadata ######################################################
###############################################################################################################################

file_path <- "/grid/cheadle/home/qianyu/Epilepsy/Blood/count_1lane"
save_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
demux_path <- "/grid/cheadle/home/qianyu/Epilepsy/Blood/genetic_demux_1lane"


datalist <- c(
  "Cheadle_QL02_PBMC_1", "Cheadle_QL02_PBMC_2", "Cheadle_QL02_PBMC_3", "Cheadle_QL02_PBMC_4",
  "Cheadle_QL02_PBMC_5", "Cheadle_QL02_PBMC_6", "Cheadle_QL02_PBMC_7", "Cheadle_QL02_PBMC_8"
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
  patient.file <- paste0(demux_path, "/", sample_id, "_outdir/vireo_results/donor_ids.tsv")
  patient <- read.table(patient.file, sep = "\t", header = TRUE)
  
  
  
  
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
    
    ## patient source
    a <- colnames(obj)
    matched_patient_info <- patient$donor_id[match(a, patient$cell)]
    obj@meta.data$patient <- matched_patient_info
    
    matched_indices <- match(a, patient$cell)
    non_matched_cells <- a[is.na(matched_indices)]
    
    pdf(file = paste0(save_path, "/", sample_id, ".pdf"), width = 20, height = 4)
    p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal", "log10GenesPerUMI"), ncol = 5)
    print(p1)
    dev.off() # Close the PDF device
    #vlnplots_list[[sample_id]] <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    ## QC fine
    # obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 2.5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.8)
    
    obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 5 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.ribosomal < 50 & log10GenesPerUMI > 0.8)
    # obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 10 & nCount_RNA > 500 & log10GenesPerUMI > 0.8)
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- ScaleData(obj)
    assign(sample_id, obj)
    # save(obj, file = paste0(paste0(save_path, "/", sample_id, ".RData")))
    rm(obj)
  } else {
    cat("File does not exist:", data.file, "\n")
  }
}




test <- merge(Cheadle_QL02_PBMC_1, 
              y = c(Cheadle_QL02_PBMC_2, Cheadle_QL02_PBMC_3, Cheadle_QL02_PBMC_4, 
                    Cheadle_QL02_PBMC_5, Cheadle_QL02_PBMC_6, Cheadle_QL02_PBMC_7, Cheadle_QL02_PBMC_8),
              add.cell.ids = c("Lane1", "Lane2", "Lane3", "Lane4", "Lane5", "Lane6", "Lane7", "Lane8"),
              project = "Epilepsy",
              merge.data = TRUE)

layersList <- lapply(test@assays$RNA@layers, function(x) {dim(x)})
layers_to_null <- names(layersList)[!grepl("^counts", names(layersList))]
test@assays$RNA@layers[layers_to_null] <- NULL

save(test, file = paste0(save_path, "/step1.0.cheadle.Blood.RData"))