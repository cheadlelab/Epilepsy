rm(list = ls(all = TRUE))

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
save_path <- "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V7"

load(file.path(file_path, "step2.1.Comparisons.sc2.RData"))  # final_results

min_lfc     <- log2(1.5)
padj_cutoff <- 0.05

collect_deg_all <- function(final_results, condition, direction = c("Up", "Down"),
                            min_lfc = log2(1.5), padj_cutoff = 0.05) {
  direction <- match.arg(direction)
  out_list <- list()
  
  for (celltype in names(final_results)) {
    if (is.null(final_results[[celltype]][[condition]]) ||
        is.null(final_results[[celltype]][[condition]][["deg"]])) next
    
    df <- final_results[[celltype]][[condition]][["deg"]]
    
    if (!"symbol" %in% colnames(df)) {
      if ("gene" %in% colnames(df)) {
        df$symbol <- df$gene
      } else if (!is.null(rownames(df))) {
        df$symbol <- rownames(df)
      } else {
        next
      }
    }
    
    df <- df %>% filter(!grepl("^MT-", symbol))
    
    if (direction == "Up") {
      sig <- df %>%
        filter(avg_log2FC > min_lfc, p_val_adj < padj_cutoff) %>%
        arrange(desc(avg_log2FC)) %>%
        mutate(
          Direction     = "Upregulated",
          negLog10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
          Rank_log2FC   = row_number(),
          CellType      = celltype
        )
    } else {
      sig <- df %>%
        filter(avg_log2FC < -min_lfc, p_val_adj < padj_cutoff) %>%
        arrange(avg_log2FC) %>%
        mutate(
          Direction     = "Downregulated",
          negLog10_padj = -log10(pmax(p_val_adj, .Machine$double.xmin)),
          Rank_log2FC   = row_number(),
          CellType      = celltype
        )
    }
    
    if (nrow(sig) > 0) out_list[[celltype]] <- sig
  }
  
  if (length(out_list) == 0) return(data.frame())
  
  all_df <- bind_rows(out_list)
  
  front_cols <- c("CellType", "symbol", "avg_log2FC", "Rank_log2FC", "p_val_adj", "Direction", "negLog10_padj")
  all_df <- all_df[, c(intersect(front_cols, colnames(all_df)),
                       setdiff(colnames(all_df), front_cols)), drop = FALSE]
  
  all_df
}

focal_up   <- collect_deg_all(final_results, "Focal_vs_Distal", "Up",   min_lfc, padj_cutoff)
focal_down <- collect_deg_all(final_results, "Focal_vs_Distal", "Down", min_lfc, padj_cutoff)
stim_up    <- collect_deg_all(final_results, "Stimulated_vs_Distal", "Up",   min_lfc, padj_cutoff)
stim_down  <- collect_deg_all(final_results, "Stimulated_vs_Distal", "Down", min_lfc, padj_cutoff)

out_xlsx <- file.path(save_path, "SM.DEG.xlsx")

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
message("Saved: ", out_xlsx)
message("Rows: Focal_Up=", nrow(focal_up),
        " | Focal_Down=", nrow(focal_down),
        " | Stim_Up=", nrow(stim_up),
        " | Stim_Down=", nrow(stim_down))
