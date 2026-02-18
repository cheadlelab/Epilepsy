rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)

# Load the Seurat object
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"

load(paste0(file_path, "/step1.2.umap.fine.RData"))
Idents(test) <- "condition"

table(test@meta.data$Annotation.fine, test@meta.data$patient)
# Create lists to store differential expression results
non_regressed <- list()  # for default test results
patient_regressed <- list()  # for MAST test results
filtered_non_regressed <- list()
filtered_patient_regressed <- list()

# Loop over each unique cell type
for (cell in unique(test@meta.data$Annotation.fine)) {
  
  if (grepl("/", cell)) {
    warning(paste("Skipping", cell, "- contains '/' in name"))
    next
  }
  
  cat("Processing cell type:", cell, "\n")
  
  cell_metadata <- test@meta.data[test@meta.data$Annotation.fine == cell, ]
  patient_counts <- table(cell_metadata$patient)
  qualified_patients <- names(patient_counts)[patient_counts > 100]

  cell_metadata <- cell_metadata[cell_metadata$patient %in% qualified_patients, ]

  if (!all(c("epileptic", "healthy") %in% unique(cell_metadata$condition))) {
    warning(paste("Skipping", cell, "- not all conditions present among qualified patients"))
    next
  }
  
  epileptic_patients <- unique(cell_metadata$patient[cell_metadata$condition == "epileptic"])
  healthy_patients   <- unique(cell_metadata$patient[cell_metadata$condition == "healthy"])
  
  if (length(epileptic_patients) <= 3 || length(healthy_patients) <= 3) {
    warning(paste("Skipping", cell, "- insufficient patients in one or both conditions (epileptic:",
                  length(epileptic_patients), ", healthy:", length(healthy_patients), ")"))
    next
  }

  # Subset the object for the current cell type
  cell_subset <- subset(test, subset = Annotation.fine == cell & patient %in% qualified_patients)
  
  pseudo_test <- AggregateExpression(object = cell_subset, 
                                     group.by = c("condition", "patient"),
                                     assays = "RNA",
                                     slot = "counts",
                                     return.seurat = T)

  Idents(pseudo_test) <- "condition"
  
  # Differential expression analysis using the default test
  deg1 <- FindMarkers(
    object = pseudo_test, 
    ident.1 = "epileptic",
    ident.2 = "healthy",
    logfc.threshold = log2(1.2),
    min.pct = 0.1
  )
  
  # Differential expression analysis using MAST with patient as a latent variable
  deg2 <- FindMarkers(
    object = pseudo_test, 
    ident.1 = "epileptic",
    ident.2 = "healthy",
    logfc.threshold = log2(1.2),
    min.pct = 0.1,
    test.use = "MAST",
    latent.vars = "patient"
  )
  
  # Store the results
  non_regressed[[cell]] <- deg1
  patient_regressed[[cell]] <- deg2
  
  # filtered_non_regressed[[cell]] <- deg1[deg1$p_val_adj < 0.05, ]
  # filtered_patient_regressed[[cell]] <- deg2[deg2$p_val_adj < 0.05, ]
  
}

# Optionally, save the results to a file
save(non_regressed, patient_regressed, file = paste0(file_path, "/step2.2.DEG.pseudobulk.RData"))
