rm(list = ls(all = TRUE))

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
# file_path <- "/grid/Epilepsy/output"
# load(paste0(file_path, "/step2.1.AllComparisons.RData"))
load(paste0(file_path, "/step2.1.Comparisons.sc.RData"))

cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", "Exc_L5IT_GRIN3A",
  "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", 
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", 
  "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier", "Astrocyte", "Microglia", 
  "Olg", "OPC", "Macrophage", "NK/T", "Vascular"
)

cell_types <- rev(cell_types)



condition <- "Focal_vs_Distal"
minGSSize   <- 5    # Minimum GO term size (number of genes in the term) to keep
maxGSSize   <- 2000  # Maximum GO term size to keep (filter out very broad terms)

go_all <- list()

for (celltype in cell_types) {
  
  ego_up <- NULL
  ego_down <- NULL
  
  test.markers1 <- final_results[[celltype]][[condition]][["deg"]] #deg upregulated downregulated
  test.markers1$symbol <- test.markers1$gene
  
  genes.up <- test.markers1$symbol[test.markers1$avg_log2FC > log2(1.5) & test.markers1$p_val_adj < 0.05]
  genes.down <- test.markers1$symbol[test.markers1$avg_log2FC < -log2(1.5) & test.markers1$p_val_adj < 0.05]
  
  
  if (length(genes.up) > 0) {
    # Perform GO enrichment for all three ontologies: BP, CC, MF
    ego_up <- enrichGO(
      gene = genes.up,
      OrgDb = org.Hs.eg.db, 
      keyType = "SYMBOL",
      ont = "BP",  # BP, CC, MF
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize
    )
  }
  if (length(genes.down) > 0) {
    # Perform GO enrichment for all three ontologies: BP, CC, MF
    ego_down <- enrichGO(
      gene = genes.down,
      OrgDb = org.Hs.eg.db, 
      keyType = "SYMBOL",
      ont = "BP",  # BP, CC, MF
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize
    )
  }
  
  
  if (exists("ego_up") && 
      !is.null(ego_up) && 
      nrow(ego_up@result) > 0) {
    ego_up@result$Direction <- "Upregulated"
    ego_up_sub <- ego_up@result[ego_up@result$p.adjust < 0.05, ]
    
  } else {
    ego_up_sub <- data.frame()
    warning("Upregulated GO results are missing or empty.")
  }
  
  if (exists("ego_down") && 
      !is.null(ego_down) && 
      nrow(ego_down@result) > 0) {
    ego_down@result$Direction <- "Downregulated"
    ego_down_sub <- ego_down@result[ego_down@result$p.adjust < 0.05, ]
    
  } else {
    ego_down_sub <- data.frame()
    warning("Downregulated GO results are missing or empty.")
  }
  
  combined_go <- NULL
  
  if (!is.null(ego_up_sub) && nrow(ego_up_sub) > 0 ||
      !is.null(ego_down_sub) && nrow(ego_down_sub) > 0) {
    
    combined_go <- dplyr::bind_rows(ego_up_sub, ego_down_sub) %>%
      dplyr::arrange(p.adjust)
    
    go_all[[celltype]] <- combined_go
    
  } else {
    message("Skip ", celltype, ": no GO terms found.")
    go_all[[celltype]] <- data.frame()
  }
}

go_all_compact <- purrr::compact(go_all)
go_all_nonempty <- go_all_compact[ vapply(go_all_compact, nrow, integer(1)) > 0 ]
go_merged <- dplyr::bind_rows(go_all_nonempty, .id = "CellType")


#####################################
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
})

stopifnot(all(c("ID","Description","CellType","Direction","p.adjust") %in% colnames(go_merged)))

go_summary <- go_merged %>%
  group_by(ID, Description) %>%
  summarise(
    total_occurrences = n(),                            
    n_celltypes       = n_distinct(CellType),             
    n_up              = sum(Direction == "Upregulated", na.rm = TRUE),
    n_down            = sum(Direction == "Downregulated", na.rm = TRUE),
    celltypes         = paste(sort(unique(CellType)), collapse = "; "),
    min_p             = min(p.adjust, na.rm = TRUE),        
    median_p          = median(p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_occurrences), min_p)

go_summary_unique_ct <- go_merged %>%
  distinct(ID, Description, CellType) %>% 
  count(ID, Description, name = "n_celltypes_unique") %>%
  arrange(desc(n_celltypes_unique))

go_summary_full <- go_summary %>%
  left_join(go_summary_unique_ct, by = c("ID","Description")) %>%
  relocate(n_celltypes_unique, .after = n_celltypes)




setwd("C:/Users/thech/Desktop")
write.csv(go_merged, "C:/Users/thech/Desktop/tt.csv", row.names = FALSE)

write_csv(go_summary_full, file.path("C:/Users/thech/Desktop", paste0("GO_summary.csv")))



A_sorted <- sort(table(go_merged$ID), decreasing = TRUE)



suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
block_ids <- list(
  "Metabolism" = c(
    "GO:0033993", # response to lipid
    "GO:0051248", # negative regulation of protein metabolic process
    "GO:0046034", # ATP metabolic process
    "GO:0010563", # negative regulation of phosphorus metabolic process
    "GO:0006163"  # purine nucleotide metabolic process
  ),
  "Stress & Death" = c(
    "GO:0006915", # apoptotic process
    "GO:0043067", # regulation of programmed cell death
    "GO:0043069", # negative regulation of programmed cell death
    "GO:0061684", # chaperone-mediated autophagy
    "GO:0043066", # negative regulation of apoptotic process
    "GO:0006914", # autophagy
    "GO:0016236", # macroautophagy
    "GO:0006979"  # response to oxidative stress
  ),
  "Development & Remodeling" = c(
    "GO:0048699", # generation of neurons
    "GO:0031175", # neuron projection development
    "GO:0030030", # cell projection organization
    "GO:0030182", # neuron differentiation
    "GO:0048858", # cell projection morphogenesis
    "GO:0140238", # presynaptic endocytosis
    "GO:0006897", # endocytosis
    "GO:0099536", # synaptic signaling
    "GO:0007409", # axonogenesis
    "GO:0042552"  # myelination
  ),
  "Signaling & Transport" = c(
    "GO:0051592", # response to calcium ion
    "GO:0009725", # response to hormone
    "GO:0032870", # cellular response to hormone stimulus
    "GO:1901379"	# regulation of potassium ion transmembrane transport
  ),
  "Homeostasis & Regulation" = c(
    "GO:0042592", # homeostatic process
    "GO:0061077", # chaperone-mediated protein folding
    "GO:0034605"  # cellular response to heat
  ),
  "Immunity" = c(
    "GO:0034097", # response to cytokine
    "GO:0042110", # T cell activation
    "GO:0046649", # lymphocyte activation
    "GO:0002682"  # regulation of immune system process
  )
)

### New version
block_ids <- list(
  "Metabolism" = c(
    "GO:0033993", # response to lipid
    "GO:0051248", # negative regulation of protein metabolic process
    "GO:0046034", # ATP metabolic process
    "GO:0010563" # negative regulation of phosphorus metabolic process
  ),
  "Stress & Death" = c(
    "GO:0006915", # apoptotic process
    "GO:0043067", # regulation of programmed cell death
    "GO:0043069", # negative regulation of programmed cell death
    "GO:0061684", # chaperone-mediated autophagy
    "GO:0043066", # negative regulation of apoptotic process
    "GO:0006914", # autophagy
    "GO:0016236", # macroautophagy
    "GO:0006979"  # response to oxidative stress
  ),
  "Development & Remodeling" = c(
    "GO:0048699", # generation of neurons
    "GO:0031175", # neuron projection development
    "GO:0030030", # cell projection organization
    "GO:0030182", # neuron differentiation
    "GO:0048858", # cell projection morphogenesis
    "GO:0140238", # presynaptic endocytosis
    "GO:0006897", # endocytosis
    "GO:0099536", # synaptic signaling
    "GO:0007409", # axonogenesis
    "GO:0042552"  # myelination
  )
)

two_line <- function(s) sub("\\s+(\\S+)$", "\n\\1", s)
block_labels <- setNames(vapply(names(block_ids), two_line, character(1)),
                         names(block_ids))

block_map <- stack(block_ids) %>%
  dplyr::rename(ID = values, Block = ind)

df_min <- go_merged %>%
  dplyr::filter(ID %in% block_map$ID) %>%
  dplyr::group_by(CellType, ID) %>%
  dplyr::summarise(p = suppressWarnings(min(p.adjust, na.rm = TRUE)), .groups = "drop")

df_min$p[!is.finite(df_min$p)] <- NA_real_

all_ct <- factor(cell_types, levels = cell_types)                
all_id <- factor(block_map$ID, levels = unlist(block_ids))       

df_full <- tidyr::expand_grid(CellType = all_ct, ID = all_id) %>%
  dplyr::left_join(df_min, by = c("CellType","ID")) %>%
  dplyr::mutate(score = ifelse(!is.na(p) & p < 0.05, -log10(p), NA_real_)) %>%
  dplyr::left_join(block_map, by = "ID")

df_full$Block    <- factor(df_full$Block, levels = names(block_ids))
df_full$CellType <- factor(df_full$CellType, levels = cell_types)
df_full$ID       <- factor(df_full$ID, levels = unlist(block_ids))


# --- drop CellType rows where all scores are NA ---
keep_ct <- df_full %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(any_value = any(!is.na(score)), .groups = "drop") %>%
  dplyr::filter(any_value) %>%
  dplyr::pull(CellType)

df_full <- df_full %>% dplyr::filter(CellType %in% keep_ct)

# reset y-axis order after filtering
df_full$CellType <- factor(df_full$CellType,
                           levels = cell_types[cell_types %in% keep_ct])




p <- ggplot(df_full, aes(x = ID, y = CellType, fill = score)) +
  geom_tile(color = "grey90", size = 0.2) +
  facet_grid(
    . ~ Block,
    scales = "free_x",
    space  = "free_x",
    labeller = labeller(Block = as_labeller(block_labels))
  ) +
  scale_fill_gradientn(
    colours  = c("#FFFFFF", "#fee8c8", "#fdbb84", "#e34a33"),
    na.value = "grey85",
    name     = expression(-log[10](p.adjust))
  ) +
  labs(x = "GO term IDs (grouped by module)", y = "Cell types",
       title = "GO enrichment heatmap (Focal vs Distal)") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid        = element_blank(),
    strip.placement   = "outside",
    strip.text.x      = element_text(size = 6, face = "bold", lineheight = 0.9),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    plot.title        = element_text(face = "bold")
  )

print(p)


go_terms_annot <- go_merged %>%
  filter(ID %in% unlist(block_ids)) %>%
  distinct(ID, Description)

id_order <- unlist(block_ids)
go_terms_annot$ID <- factor(go_terms_annot$ID, levels = id_order)
go_terms_annot <- go_terms_annot[order(go_terms_annot$ID), ]

p_annot <- ggplot(go_terms_annot, aes(y = rev(ID), x = 1)) +
  geom_text(aes(label = paste0(ID, "  ", Description)), hjust = 0, size = 2.5) +
  scale_y_discrete(position = "right") +
  theme_void(base_size = 6) +
  theme(
    plot.margin = margin(5, 20, 5, 5), 
    axis.text.y = element_blank()
  )


################################################################################

final_plot <- p + p_annot + plot_layout(widths = c(3, 1))
final_plot

p
# ===================== 5) 可选：导出图 =====================

pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.GO.Module.pdf", width = 5, height = 3)
print(p)
dev.off() # Close the PDF device
