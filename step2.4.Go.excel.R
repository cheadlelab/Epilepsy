rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(ggplot2)
##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"
# file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output/Blood"
# file_path <- "/grid/Epilepsy/output"
# load(paste0(file_path, "/step2.1.AllComparisons.RData"))
load(paste0(file_path, "/step2.1.DEG.singlecell.RData"))


final_results <- filtered_patient_regressed

ont_choice  <- "BP"                    # "BP", "CC", or "MF"
min_lfc     <- log2(1.5)
padj_cutoff <- 0.05

###################################################################################################################################################
##################################################################### GO plot #####################################################################
###################################################################################################################################################

Go_up   <- list()       # list of per-celltype sig GO tables (Up)
Go_down <- list()       # list of per-celltype sig GO tables (Down)
Up_all   <- data.frame()  # combined Up across CTs
Down_all <- data.frame()  # combined Down across CTs


## stricter ORA with background/universe + size/overlap filters
min_deg_n   <- 10    # Minimum number of DEGs required to run ORA; skip ORA if fewer than this
min_overlap <- 3     # Minimum overlap (Count) of DEGs per GO term (require Count >= min_overlap)
minGSSize   <- 5    # Minimum GO term size (number of genes in the term) to keep
maxGSSize   <- 2000  # Maximum GO term size to keep (filter out very broad terms)

# celltype <- "CD14_Mono"

for (celltype in names(final_results)) {
  ego_up <- NULL
  ego_down <- NULL
  
  if (is.null(final_results[[celltype]])) next
  
  test.markers1 <- final_results[[celltype]]
  test.markers1$symbol <- rownames(test.markers1)
  
  if (!"symbol" %in% colnames(test.markers1)) {
    if ("gene" %in% colnames(test.markers1)) {
      test.markers1$symbol <- test.markers1$gene
    } else if (!is.null(rownames(test.markers1))) {
      test.markers1$symbol <- rownames(test.markers1)
    }
  }
  
  bg_universe <- unique(test.markers1$symbol[!is.na(test.markers1$symbol)])
  
  genes.up <- unique(test.markers1$symbol[
    test.markers1$avg_log2FC >  min_lfc & test.markers1$p_val_adj < padj_cutoff
  ])
  genes.down <- unique(test.markers1$symbol[
    test.markers1$avg_log2FC < -min_lfc & test.markers1$p_val_adj < padj_cutoff
  ])
  
  # genes.up <- genes.up[!grepl("^HSP", genes.up)]
  # genes.down <- genes.down[!grepl("^HSP", genes.down)]
  
  if (length(genes.up) >= min_deg_n) {
    ego_up <- tryCatch(
      enrichGO(
        gene          = genes.up,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont_choice,     # BP/CC/MF
        universe      = bg_universe,    
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
      ),
      error = function(e) NULL
    )
  }
  
  if (length(genes.down) >= min_deg_n) {
    ego_down <- tryCatch(
      enrichGO(
        gene          = genes.down,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = ont_choice,
        universe      = bg_universe,
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
      ),
      error = function(e) NULL
    )
  }
  
  if (!is.null(ego_up) && nrow(ego_up@result) > 0) {
    ego_up@result$Direction <- "Upregulated"
    ego_up_sub <- subset(ego_up@result, p.adjust < padj_cutoff & Count >= min_overlap)
    if (nrow(ego_up_sub) > 0) {
      Go_up[[celltype]] <- ego_up_sub
    } else {
      warning(sprintf("[%s] Upregulated GO: filtered to empty (p.adjust/Count).", celltype))
    }
  } else {
    warning(sprintf("[%s] Upregulated GO results are missing or empty.", celltype))
  }
  
  if (!is.null(ego_down) && nrow(ego_down@result) > 0) {
    ego_down@result$Direction <- "Downregulated"
    ego_down_sub <- subset(ego_down@result, p.adjust < padj_cutoff & Count >= min_overlap)
    if (nrow(ego_down_sub) > 0) {
      Go_down[[celltype]] <- ego_down_sub
    } else {
      warning(sprintf("[%s] Downregulated GO: filtered to empty (p.adjust/Count).", celltype))
    }
  } else {
    warning(sprintf("[%s] Downregulated GO results are missing or empty.", celltype))
  }
}

#################################################################################
##################################### Excel ##################################### 
#################################################################################
# ======================= Write Excel outputs (each cell type = one sheet) =======================
# install.packages("openxlsx") # uncomment if not installed
library(openxlsx)

# helper vars for safe Excel sheet names (no functions)
.invalid_pat <- "[:\\\\/\\?\\*\\[\\]]"
.max_sheet_len <- 31

# ------------------- Upregulated workbook -------------------
if (exists("Go_up") && is.list(Go_up) && length(Go_up) > 0) {
  wb_up <- createWorkbook()
  
  up_sheet_names <- names(Go_up)
  if (is.null(up_sheet_names) || any(up_sheet_names == "")) {
    up_sheet_names <- paste0("Sheet", seq_along(Go_up))
  }
  
  # sanitize + truncate + make unique
  up_safe_names <- character(length(up_sheet_names))
  for (i in seq_along(up_sheet_names)) {
    nm <- gsub(.invalid_pat, "_", up_sheet_names[i])
    nm <- substr(nm, 1, .max_sheet_len)
    if (is.na(nm) || nm == "") nm <- paste0("Sheet", i)
    up_safe_names[i] <- nm
  }
  up_safe_names <- make.unique(up_safe_names)
  
  for (i in seq_along(Go_up)) {
    df <- Go_up[[i]]
    if (is.data.frame(df) && nrow(df) > 0) {
      addWorksheet(wb_up, up_safe_names[i])
      writeData(wb_up, up_safe_names[i], df)
    }
  }
  
  saveWorkbook(
    wb_up,
    file.path(file_path, paste0("step2.4.GO_", ont_choice, "_Upregulated.xlsx")),
    overwrite = TRUE
  )
}

# ------------------- Downregulated workbook -------------------
if (exists("Go_down") && is.list(Go_down) && length(Go_down) > 0) {
  wb_dn <- createWorkbook()
  
  dn_sheet_names <- names(Go_down)
  if (is.null(dn_sheet_names) || any(dn_sheet_names == "")) {
    dn_sheet_names <- paste0("Sheet", seq_along(Go_down))
  }
  
  # sanitize + truncate + make unique
  dn_safe_names <- character(length(dn_sheet_names))
  for (i in seq_along(dn_sheet_names)) {
    nm <- gsub(.invalid_pat, "_", dn_sheet_names[i])
    nm <- substr(nm, 1, .max_sheet_len)
    if (is.na(nm) || nm == "") nm <- paste0("Sheet", i)
    dn_safe_names[i] <- nm
  }
  dn_safe_names <- make.unique(dn_safe_names)
  
  for (i in seq_along(Go_down)) {
    df <- Go_down[[i]]
    if (is.data.frame(df) && nrow(df) > 0) {
      addWorksheet(wb_dn, dn_safe_names[i])
      writeData(wb_dn, dn_safe_names[i], df)
    }
  }
  
  saveWorkbook(
    wb_dn,
    file.path(file_path, paste0("step2.4.GO_", ont_choice, "_Downregulated.xlsx")),
    overwrite = TRUE
  )
}
