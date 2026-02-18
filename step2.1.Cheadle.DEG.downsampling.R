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

# remove "Other"
test <- subset(test, subset = Annotation.fine != "Other")

# only keep the 3 conditions we care about (safe-guard)
conds_keep <- c("Focal", "Distal", "Stimulated")
test <- subset(test, subset = condition %in% conds_keep)

# Check the distribution of conditions and available cell types
print(table(test@meta.data$condition, test@meta.data$patient))
print(unique(test@meta.data$Annotation.fine))

##################################################################################################################################
# Relevel cell types
##################################################################################################################################
test@meta.data$Annotation.new <- as.character(test@meta.data$Annotation.fine)

custom_order <- c(
  "Exc_L23IT", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", "Exc_L5IT_GRIN3A", "Exc_L5IT_GRIN3A-",
  "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", "Exc_L234IT_ERBB4",
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier",
  "Astrocyte", "Microglia", "Olg", "OPC",
  "Macrophage", "NK/T", "Vascular"
)
test@meta.data$Annotation.new <- factor(test@meta.data$Annotation.new, levels = custom_order)
print(unique(test@meta.data$Annotation.new))

##################################################################################################################################
# Downsampling: per cell type x condition, sample up to 500 cells
##################################################################################################################################
set.seed(123)

meta <- test@meta.data
meta$cell <- rownames(meta)

keep_cells <- meta %>%
  filter(!is.na(Annotation.new)) %>%
  group_by(Annotation.new, condition) %>%
  group_modify(~{
    n <- nrow(.x)
    if (n <= 500) {
      .x
    } else {
      .x[sample.int(n, 500), , drop = FALSE]
    }
  }) %>%
  pull(cell)

test <- subset(test, cells = keep_cells)

cat("\nAfter downsampling (<=500 each Annotation.new x condition):\n")
print(table(test$Annotation.new, test$condition))

##################################################################################################################################
# Analyze by Cell Type with Patient Correction (latent.vars="patient") for multiple comparisons
##################################################################################################################################
comparisons <- list(
  "Focal_vs_Distal"      = c("Focal", "Distal"),
  "Stimulated_vs_Distal" = c("Stimulated", "Distal")
)

final_results <- list()

for (cell in custom_order) {
  cat("----------------------------------------------------------\n")
  cat("Processing cell type:", cell, "\n")
  
  cell_subset <- subset(test, subset = Annotation.new == cell)
  
  if (ncol(cell_subset) == 0) {
    cat("No cells found for cell type:", cell, "\n")
    next
  }
  
  DefaultAssay(cell_subset) <- "RNA"
  Idents(cell_subset) <- cell_subset$condition
  
  final_results[[cell]] <- list()
  
  for (comp_name in names(comparisons)) {
    conds <- comparisons[[comp_name]]
    cat("Performing comparison:", comp_name, "(", conds[1], "vs", conds[2], ")\n")
    
    present_conditions <- unique(cell_subset$condition)
    if (!all(conds %in% present_conditions)) {
      cat("  Skipping", comp_name, "for cell type", cell,
          "because one or both conditions are not present.\n")
      next
    }
    
    deg <- FindMarkers(
      object      = cell_subset,
      ident.1     = conds[1],
      ident.2     = conds[2],
      min.pct     = 0.1,
      test.use    = "MAST",
      latent.vars = "patient"
    )
    
    deg$gene <- rownames(deg)
    
    up_deg   <- deg[deg$avg_log2FC > 0 & deg$p_val_adj < 0.05, ]
    down_deg <- deg[deg$avg_log2FC < 0 & deg$p_val_adj < 0.05, ]
    
    up_deg   <- up_deg[order(up_deg$avg_log2FC, decreasing = TRUE), ]
    down_deg <- down_deg[order(down_deg$avg_log2FC, decreasing = FALSE), ]
    
    final_results[[cell]][[comp_name]] <- list(
      deg          = deg,
      upregulated  = up_deg,
      downregulated= down_deg
    )
    
    cat("Completed comparison:", comp_name, "for cell type:", cell, "\n")
  }
}

##################################################################################################################################
# Save Final Results
##################################################################################################################################
print(final_results)

save(final_results, file = paste0(file_path, "/step2.1.Comparisons.sc.downsampling.RData"))
# optional reload check
# load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))




####################################################################################################################################

rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
save_path <- "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V7"

in_rdata <- file.path(file_path, "step2.1.Comparisons.sc.downsampling.RData")
load(in_rdata)  # loads: final_results

min_lfc     <- log2(1.5)
padj_cutoff <- 0.05

collect_deg_all <- function(final_results, comparison, direction = c("Up", "Down"),
                            min_lfc = log2(1.5), padj_cutoff = 0.05,
                            rm_mt = TRUE) {
  direction <- match.arg(direction)
  out_list <- list()
  
  for (celltype in names(final_results)) {
    # safety checks
    if (is.null(final_results[[celltype]][[comparison]])) next
    if (is.null(final_results[[celltype]][[comparison]][["deg"]])) next
    
    df <- final_results[[celltype]][[comparison]][["deg"]]
    
    if (is.null(df) || nrow(df) == 0) next
    
    # make sure gene symbol column exists
    if (!"symbol" %in% colnames(df)) {
      if ("gene" %in% colnames(df)) {
        df$symbol <- df$gene
      } else if (!is.null(rownames(df))) {
        df$symbol <- rownames(df)
      } else {
        next
      }
    }
    
    # remove MT genes if requested
    if (rm_mt) df <- df %>% filter(!grepl("^MT-", symbol))
    
    # require key columns
    req_cols <- c("avg_log2FC", "p_val_adj", "symbol")
    if (!all(req_cols %in% colnames(df))) next
    
    if (direction == "Up") {
      sig <- df %>%
        filter(avg_log2FC >  min_lfc, p_val_adj < padj_cutoff) %>%
        arrange(desc(avg_log2FC)) %>%
        mutate(
          CellType      = celltype,
          Direction     = "Upregulated",
          negLog10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
          Rank_log2FC   = row_number()
        )
    } else {
      sig <- df %>%
        filter(avg_log2FC < -min_lfc, p_val_adj < padj_cutoff) %>%
        arrange(avg_log2FC) %>%
        mutate(
          CellType      = celltype,
          Direction     = "Downregulated",
          negLog10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
          Rank_log2FC   = row_number()
        )
    }
    
    if (nrow(sig) > 0) out_list[[celltype]] <- sig
  }
  
  if (length(out_list) == 0) return(data.frame())
  
  all_df <- bind_rows(out_list)
  
  # reorder columns (front)
  front_cols <- c("CellType", "symbol", "avg_log2FC", "Rank_log2FC", "p_val_adj", "Direction", "negLog10_padj")
  all_df <- all_df[, c(intersect(front_cols, colnames(all_df)),
                       setdiff(colnames(all_df), front_cols)), drop = FALSE]
  
  all_df
}

# 4 outputs
focal_up   <- collect_deg_all(final_results, "Focal_vs_Distal",      "Up",   min_lfc, padj_cutoff)
focal_down <- collect_deg_all(final_results, "Focal_vs_Distal",      "Down", min_lfc, padj_cutoff)
stim_up    <- collect_deg_all(final_results, "Stimulated_vs_Distal", "Up",   min_lfc, padj_cutoff)
stim_down  <- collect_deg_all(final_results, "Stimulated_vs_Distal", "Down", min_lfc, padj_cutoff)

# write xlsx
out_xlsx <- file.path(save_path, "SM.DEG.downsampling.xlsx")
wb <- createWorkbook()

addWorksheet(wb, "Focal_Up_ALL")
writeData(wb, "Focal_Up_ALL", focal_up)

addWorksheet(wb, "Focal_Down_ALL")
writeData(wb, "Focal_Down_ALL", focal_down)

addWorksheet(wb, "Stimulated_Up_ALL")
writeData(wb, "Stimulated_Up_ALL", stim_up)

addWorksheet(wb, "Stimulated_Down_ALL")
writeData(wb, "Stimulated_Down_ALL", stim_down)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Loaded RData: ", in_rdata)
message("Saved: ", out_xlsx)
message("Rows: Focal_Up=", nrow(focal_up),
        " | Focal_Down=", nrow(focal_down),
        " | Stim_Up=", nrow(stim_up),
        " | Stim_Down=", nrow(stim_down))

