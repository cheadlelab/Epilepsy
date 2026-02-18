rm(list = ls(all = TRUE))

#load packages
library(Seurat)
library(dplyr)
library(jsonlite)
library(dplyr)
library(arrow)
library(spacexr)
library(ggplot2)
library(scales)
options(future.globals.maxSize = 128 * 1024^3)

save.path <- "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/"

load(paste0(save.path, "step2.umap.RData"))


orig.ident <- c("1C SA", "4E", "1C SB", "5C", "2E", "3E SB", "1E", "4C", "5E", "2C", "3E SA")
patient <- c(1, 4 ,1, 5, 2, 3, 1, 4, 5, 2, 3)
condition <- c("Distal", "Focal", "Distal", "Distal", "Focal", "Focal", "Focal", "Distal", "Focal", "Distal", "Focal")

test@meta.data$patient <- patient[match(test@meta.data$orig.ident, orig.ident)]
test@meta.data$condition <- condition[match(test@meta.data$orig.ident, orig.ident)]



Idents(test) <- "condition"
markers_list <- list()

for (i in unique(test@meta.data$predicted.celltype)) {
  tmp <- subset(test, subset = predicted.celltype == i)
  markers_list[[i]] <- FindMarkers(tmp, 
                                   ident.1 = "Focal",
                                   ident.2 = "Distal",
                                   logfc.threshold = log2(1.5),
                                   min.pct = 0.1)
}

###################################################################################
library(openxlsx)
library(dplyr)

# Output path
file_path <- save.path
condition <- "Focal_vs_Distal"  # You can change this if needed

# Create Excel: one sheet per cell type
if (length(markers_list) > 0) {
  wb <- createWorkbook()
  
  # Clean up sheet names to be Excel-safe
  sheet_names <- names(markers_list)
  if (is.null(sheet_names) || any(sheet_names == "")) {
    sheet_names <- paste0("Sheet", seq_along(markers_list))
  }
  
  invalid_pat <- "[\\\\/:*?\\[\\]]"  # forbidden characters in Excel sheet names
  max_sheet_len <- 31
  safe_names <- character(length(sheet_names))
  
  for (i in seq_along(sheet_names)) {
    nm <- gsub(invalid_pat, "_", sheet_names[i])
    nm <- substr(nm, 1, max_sheet_len)
    if (is.na(nm) || nm == "") nm <- paste0("Sheet", i)
    safe_names[i] <- nm
  }
  safe_names <- make.unique(safe_names)
  
  # Write one sheet per cell type
  for (i in seq_along(markers_list)) {
    addWorksheet(wb, safe_names[i])
    writeData(wb, safe_names[i], markers_list[[i]], rowNames = TRUE)
  }
  
  # Save workbook
  saveWorkbook(
    wb,
    file.path(file_path, paste0("Xenium.DEG_", condition, "_All_CellTypes.xlsx")),
    overwrite = TRUE
  )
  
}

# Optional: combined ALL file
if (length(markers_list) > 0) {
  all_markers <- dplyr::bind_rows(Map(function(ct, df) dplyr::mutate(df, CellType = ct),
                                      names(markers_list), markers_list))
  write.xlsx(
    all_markers,
    file = file.path(file_path, paste0("Xenium.DEG_", condition, "_ALL.xlsx")),
    overwrite = TRUE
  )
}
