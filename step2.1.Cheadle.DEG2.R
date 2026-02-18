rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)

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

load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))
test <- JoinLayers(test)
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
DefaultAssay(test) <- "RNA"

test <- subset(test, subset = Annotation.fine != "Other")


# Check the distribution of conditions and available cell types
table(test@meta.data$condition, test@meta.data$patient)
print(unique(test@meta.data$Annotation.fine))

test@meta.data$Annotation.new <- as.character(test@meta.data$Annotation.fine)

custom_order <- c("Exc_L23IT", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", "Exc_L5IT_GRIN3A", "Exc_L5IT_GRIN3A-",
                  "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", "Exc_L234IT_ERBB4",
                  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier",
                  "Astrocyte", "Microglia", "Olg", "OPC",
                  "Macrophage", "NK/T", "Vascular")
test@meta.data$Annotation.new <- factor(test@meta.data$Annotation.new, levels = custom_order)
unique(test@meta.data$Annotation.new)

# test <- subset(test, dataset %in% c("Cheadle", "Allen_MTG"))


##################################################################################################################################
# Analyze by Cell Type with Patient Correction using latent.vars for multiple comparisons
##################################################################################################################################
# Define the comparisons as a named list

comparisons <- list(
  "Focal_vs_Distal"       = c("Focal", "Distal"),
  # "Focal_vs_Stimulated"   = c("Focal", "Stimulated"),
  "Stimulated_vs_Distal"  = c("Stimulated", "Distal")
)


final_results <- list()

# Initialize a list to store final results for each cell type and comparison

for (cell in custom_order) {
  cat("----------------------------------------------------------\n")
  cat("Processing cell type:", cell, "\n")
  
  # Subset to current cell type and include only cells from the three conditions
  cell_subset <- subset(test, subset = Annotation.new == cell)
  
  # Skip cell types with no cells
  if (ncol(cell_subset) == 0) {
    cat("No cells found for cell type:", cell, "\n")
    next
  }
  
  DefaultAssay(cell_subset) <- "RNA"
  # Set the cell identities to the condition for differential analysis
  Idents(cell_subset) <- cell_subset$condition
  
  # Initialize a list for this cell type to store each comparison
  final_results[[cell]] <- list()
  
  # Loop over the three comparisons
  for (comp_name in names(comparisons)) {
    conds <- comparisons[[comp_name]]
    cat("Performing comparison:", comp_name, "(", conds[1], "vs", conds[2], ")\n")
    
    # Check if both conditions are present in the subset
    present_conditions <- unique(cell_subset$condition)
    if (!all(conds %in% present_conditions)) {
      cat("  Skipping", comp_name, "for cell type", cell, 
          "because one or both conditions are not present.\n")
      next
    }
    
    # Perform differential expression analysis using latent.vars = "patient"
    # deg <- FindMarkers(cell_subset, 
    #                    ident.1 = conds[1],
    #                    ident.2 = conds[2],
    #                   logfc.threshold = log2(1.2),
    #                   min.pct = 0.1,
    #                   test.use = "MAST",
    #                   latent.vars = "patient")
    
    deg <- FindMarkers(cell_subset, 
                       ident.1 = conds[1],
                       ident.2 = conds[2],
                       min.pct = 0.1,
                       test.use = "MAST",
                       latent.vars = "patient")
    
    # Add gene names as a column (rownames are gene names)
    deg$gene <- rownames(deg)
    
    # Separate DEGs into upregulated and downregulated (using adjusted p-value < 0.05)
    up_deg <- deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ]
    down_deg <- deg[deg$avg_log2FC < 0 & deg$p_val_adj < 0.05, ]
    
    # Order the upregulated genes in descending order of avg_log2FC,
    # and the downregulated genes in ascending order (most negative first)
    up_deg <- up_deg[order(up_deg$avg_log2FC, decreasing = TRUE), ]
    down_deg <- down_deg[order(down_deg$avg_log2FC, decreasing = FALSE), ]
    
    # Save the complete DEG results as well as separated upregulated and downregulated lists
    final_results[[cell]][[comp_name]] <- list(deg = deg, 
                                               upregulated = up_deg, 
                                               downregulated = down_deg)
    
    cat("Completed comparison:", comp_name, "for cell type:", cell, "\n")
  }
}

##################################################################################################################################
# Display and Save Final Results
##################################################################################################################################
print(final_results)

# Set the path to save the output (modify as needed)
save(final_results, file = paste0(file_path, "/step2.1.Comparisons.sc2.RData"))
load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))

