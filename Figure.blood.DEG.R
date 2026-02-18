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
#################################################################################################################################
############################################################# scRNA #############################################################
#################################################################################################################################
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
save_path <- "/grid/cheadle/home/qianyu/Epilepsy/figures"


load(paste0(file_path, "/step1.2.umap.fine.RData"))

# Prepare the Seurat object
`%!in%` <- Negate(`%in%`)
test <- subset(test, subset = Annotation.fine %!in% c("Macrophage/platelets", "Myeloid/NKT", "RBC", "Platelets"))


test@active.ident <- as.factor(test@meta.data$orig.ident)
ordered_levels <- c(
  # Myeloid lineage (ordered by differentiation hierarchy)
  "FCGR3A+ Mono",        # Classical monocytes (CD16+)
  "Macrophage",          # Tissue-resident macrophages
  # "Macrophage/platelets",# Macrophage-platelet interacting population
  "DC",                  # Conventional dendritic cells
  "cDCs",                # Classical dendritic cells
  "pDC",                 # Plasmacytoid dendritic cells
  "Mast",                # Mast cells
  # "Myeloid/NKT",         # Myeloid-NKT hybrid population
  
  # Lymphoid lineage
  ## T-cell series (ordered by differentiation state)
  "Naive_CD4+_T",        # Naïve CD4+ T cells
  "Memory_CD4+_T",       # Memory CD4+ T cells
  "regular_T",           # Conventional T cells (unspecified subset)
  "CD8+_T_1",            # CD8+ T cell subset 1 (likely effector)
  "CD8+_T_2",            # CD8+ T cell subset 2 (memory)
  "CD8+_T_3",            # CD8+ T cell subset 3 (exhausted)
  "gamma_delta_T",       # γδ T cells
  
  ## NK-cell series 
  "NK_1",                # NK subset 1 (CD56bright)
  "NK_2",                # NK subset 2 (CD56dim)
  "NK_3",                # NK subset 3 (adaptive)
  
  ## B-cell series (developmental order)
  "Pre_B",               # Pre-B cells (immature)
  "B",                   # Mature B cells
  "Plasma_B",            # Plasma cells (antibody-secreting)
  
  # Proliferating populations
  "MKI67+_NKT"           # Proliferating NKT cells (cell cycle active)
  
  # Terminal functional elements
  # "Platelets"          # Platelets (thrombocytes)
  # "RBC"                # Red blood cells (erythrocytes)
)
test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)
print(levels(test@meta.data$Annotation.fine))

default_colors <- hue_pal()(20)
names(default_colors) <- levels(test@meta.data$Annotation.fine)

##################################################################################################################################
########################################################## DEG analysis ##########################################################
##################################################################################################################################
DEG_results <- list()
for (celltype in unique(test@meta.data$Annotation.fine)) {
  
  tt.tmp <- subset(test, Annotation.fine == celltype)
  
  tt.tmp <- AggregateExpression(
    tt.tmp,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = c("Annotation.fine", "patient", "condition")
  )
  
  Idents(tt.tmp) <- "condition"
  markers <- FindMarkers(
    object = tt.tmp,
    ident.1 = "epileptic",
    ident.2 = "healthy",
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "DESeq2" # for pseudobulk
  )
  
  markers <- markers %>%
    mutate(
      log2FoldChange = avg_log2FC,
      padj = p_val_adj,
      symbol = rownames(.)
    )
  
  
  DEG_results[[celltype]] <- markers
}

save(DEG_results, file = paste0(file_path, "/blood.DEGs.RData"))

##################################################################################################################################
########################################################## volcano plot ##########################################################
##################################################################################################################################
source("VolcanoPlot.R")
library(ggplot2)
library(ggrepel)


for (celltype in names(DEG_results)) {
  filtered_data <- DEG_results[[celltype]] %>%
    dplyr::filter(
      !grepl("^(AC|AL|AP)[0-9]+\\.[0-9]+", rownames(.)),
      !grepl("^LINC|^MIR|^SNORD|^RNU", rownames(.)),
      !is.na(padj)
    )
  
  safe_name <- gsub("[^[:alnum:]]", "_", celltype)
  pdf_file <- paste0("Volcano_", safe_name, ".pdf")
  
  pdf(file = pdf_file, width = 6, height = 6)
  
  p <- VolcanoPlot(
    dif = filtered_data,
    log2FC = log2(1.5),
    padj = 0.05,
    title = celltype,
    cols = c("#377EB8", "#E41A1C"),
    label.max = 30
  )
  print(p)
  dev.off()
}

###################################################################################################################################
########################################################## GO enrichment ##########################################################
###################################################################################################################################

up_genes <- DEG_results[[celltype]] %>%
  filter(padj < 0.05, log2FoldChange > log2(1.5)) %>%
  pull(symbol)
down_genes <- DEG_results[[celltype]] %>%
  filter(padj < 0.05, log2FoldChange < -log2(1.5)) %>%
  pull(symbol)

safe_enrichGO <- function(genes) {
  if (length(genes) == 0) {
    return(NULL)
  } else {
    tryCatch({
      enrichGO(
        gene          = genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = c("ALL"), # "BP", "MF", "CC"
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
      )
    }, error = function(e) {
      message("Enrichment failed: ", e$message)
      return(NULL)
    })
  }
}

ego.up <- safe_enrichGO(up_genes)
ego.down <- safe_enrichGO(down_genes)

extract_terms <- function(ego, direction) {
  if (is.null(ego)) return(data.frame())
  if (nrow(ego@result) == 0) return(data.frame())
  
  ego@result %>%
    filter(p.adjust < 0.05) %>%
    head(10) %>%
    mutate(
      Direction = direction,
      logP = -log10(p.adjust)
    )
}

up_terms <- extract_terms(ego.up, "Up")
down_terms <- extract_terms(ego.down, "Down")

plot_data <- rbind(up_terms, down_terms)

if (nrow(plot_data) == 0) {
  message("No significant GO terms to display.")
} else {
  max_value <- if (nrow(plot_data) > 0) max(plot_data$logP) * 1.1 else 10
  
  p <- ggplot(plot_data, aes(x = logP, 
                             y = reorder(Description, logP), 
                             fill = Direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
    scale_x_continuous(
      breaks = seq(0, ceiling(max_value), by = 2),
      labels = abs,
      limits = c(0, max_value)
    ) +
    labs(x = "-log10(Adjusted p-value)", 
         y = "GO Biological Processes") +
    theme_minimal(base_size = 12) +
    theme(
      text = element_text(color = "black", size = 12),  # 修复字号设置
      panel.grid.major.y = element_blank(),
      axis.line.x = element_line(color = "black"),
      legend.position = "none"
    )
  
  print(p)
}

ggsave(paste0("GO_",celltype,".pdf"), plot = p, width = 2, height = 1, units = "in")





