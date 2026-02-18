rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)

# Load the Seurat object
file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output/Blood"
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"

load(paste0(file_path, "/step1.2.umap.fine.RData"))
Idents(test) <- "condition"

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
  
  # Check cell population in the metadata before subsetting
  cell_metadata <- test@meta.data[test@meta.data$Annotation.fine == cell, ]
  pop_counts <- table(cell_metadata$condition)
  
  if (!all(c("epileptic", "healthy") %in% names(pop_counts))) {
    warning(paste("Skipping", cell, "- does not contain both conditions"))
    next
  }
  
  if (pop_counts["epileptic"] <= 50 || pop_counts["healthy"] <= 50) {
    warning(paste("Skipping", cell, "- insufficient cells in one or both conditions (epileptic:",
                  pop_counts["epileptic"], ", healthy:", pop_counts["healthy"], ")"))
    next
  }
  
  # Subset the object for the current cell type
  cell_subset <- subset(test, subset = Annotation.fine == cell)
  
  # Differential expression analysis using the default test
  deg1 <- FindMarkers(
    object = cell_subset, 
    ident.1 = "epileptic",
    ident.2 = "healthy",
    logfc.threshold = log2(1.2),
    min.pct = 0.1
  )
  
  # Differential expression analysis using MAST with patient as a latent variable
  deg2 <- FindMarkers(
    object = cell_subset, 
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
  
  filtered_non_regressed[[cell]] <- deg1[deg1$p_val_adj < 0.05, ]
  filtered_patient_regressed[[cell]] <- deg2[deg2$p_val_adj < 0.05, ]
  
}

# Optionally, save the results to a file
save(non_regressed, patient_regressed, filtered_non_regressed, filtered_patient_regressed, file = paste0(file_path, "/step2.1.DEG.singlecell.RData"))
load(paste0(file_path, "/step2.1.DEG.singlecell.RData"))







# install.packages("writexl")  # if needed
library(writexl)

output_dir <- file.path(file_path, "patient_regressed_xlsx")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (cell in names(patient_regressed)) {
  # skip any cell names containing a backslash
  if (grepl("\\\\", cell)) {
    warning("Skipping ", cell, " – contains backslash '\\'")
    next
  }
  
  df <- patient_regressed[[cell]]
  if (nrow(df) == 0) {
    warning("No DEGs for ", cell, " – skipping")
    next
  }
  
  # build filepath using the original cell name
  outfile <- file.path(output_dir, paste0("patient_regressed_", cell, ".xlsx"))
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
  
  # reorder by descending avg_log2FC and add gene column
  df_ord <- df[order(-df$avg_log2FC), , drop = FALSE]
  df_ord <- cbind(gene = rownames(df_ord), df_ord)
  rownames(df_ord) <- NULL
  
  # write out
  message("Writing ", basename(outfile), " (", nrow(df_ord), " rows)…")
  write_xlsx(df_ord, path = outfile)
}
