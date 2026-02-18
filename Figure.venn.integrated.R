rm(list = ls(all = TRUE))
library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)

##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
# file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output"
load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))

cell_types <- names(final_results)

cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", "Exc_L5IT_GRIN3A",
  "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", 
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", 
  "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier", "Astrocyte", "Microglia", 
  "Olg", "OPC", "Macrophage", "NK/T", "Vascular"
)

# deg <- final_results[[cell_type]][["Focal_vs_Distal"]][["deg"]]
# upregulated <- final_results[[cell_type]][["Focal_vs_Distal"]][["upregulated"]]
# downregulated <- final_results[[cell_type]][["Focal_vs_Distal"]][["downregulated"]]


get_unique_genes <- function(condition, direction) {
  all_genes <- c()
  
  for (ct in cell_types) {
    deg <- final_results[[ct]][[condition]][[direction]]
    
    deg <- deg[
      deg$p_val_adj < 0.05 &
        abs(deg$avg_log2FC) >= log2(1.5),
    ]
    
    
    genes <- deg$gene
    # genes <- genes[!grepl("^MT-", genes)]
    all_genes <- c(all_genes, genes)
  }
  
  unique(all_genes)
}

# ---- Focal ----
focal_up   <- get_unique_genes("Focal_vs_Distal", "upregulated")
focal_down <- get_unique_genes("Focal_vs_Distal", "downregulated")

# ---- Stim ----
stim_up    <- get_unique_genes("Stimulated_vs_Distal", "upregulated")
stim_down  <- get_unique_genes("Stimulated_vs_Distal", "downregulated")


gene_lists_up <- list(
  Focal = focal_up,
  Stim  = stim_up
)

gene_lists_down <- list(
  Focal = focal_down,
  Stim  = stim_down
)

# common_genes <- Reduce(intersect, gene_lists)
# common_genes


venn_plot_up <- ggVennDiagram(
  gene_lists_up,
  label_alpha = 0,
  edge_size = 0.5,
  label = "count",
  set_size = 4.333,
  category.names = names(gene_lists_up)
) +
  scale_fill_gradient(low = "#F7FBFF", high = "#E78C8C") + # 4676B5
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 7),
    plot.title = element_text(size = 7),
    plot.caption = element_text(size = 7)
  ) +
  theme(  # for ggVennDiagram-specific labels
    set.text = element_text(size = 7),
    region.text = element_text(size = 7)
  )



venn_plot_down <- ggVennDiagram(
  gene_lists_down,
  label_alpha = 0,
  edge_size = 0.5,
  label = "count",
  set_size = 4.333,
  category.names = names(gene_lists_down)
) +
  scale_fill_gradient(low = "#F7FBFF", high = "#4676B5") + # 4676B5
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 7),
    plot.title = element_text(size = 7),
    plot.caption = element_text(size = 7)
  ) +
  theme(  # for ggVennDiagram-specific labels
    set.text = element_text(size = 7),
    region.text = element_text(size = 7)
  )


pdf(file = "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V6/Figure.Focal.Stim.venn.up.pdf", width = 2.5, height = 2.5)
print(venn_plot_up)
dev.off() # Close the PDF device

pdf(file = "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V6/Figure.Focal.Stim.venn.down.pdf", width = 2.5, height = 2.5)
print(venn_plot_down)
dev.off() # Close the PDF device

##################################################################################################################################
# Barplot: per cell type stacked counts (Focal / Shared / Stim)
##################################################################################################################################
cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", "Exc_L5IT_GRIN3A",
  "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", 
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", 
  "Inh_LAMP5", "Inh_LAMP5_LHX6", "Astrocyte", "Microglia", 
  "Olg", "OPC", "Macrophage", "NK/T", "Vascular"
)
cell_types <- rev(cell_types)
library(dplyr)
library(ggplot2)

get_ct_counts <- function(condition_focal = "Focal_vs_Distal",
                          condition_stim  = "Stimulated_vs_Distal",
                          direction = c("upregulated", "downregulated"),
                          padj_cutoff = 0.05,
                          lfc_cutoff  = log2(1.5)) {
  direction <- match.arg(direction)
  
  df_list <- lapply(cell_types, function(ct) {
    
    deg_focal <- final_results[[ct]][[condition_focal]][[direction]]
    deg_stim  <- final_results[[ct]][[condition_stim ]][[direction]]
    
    deg_focal <- deg_focal %>%
      filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > lfc_cutoff)
    deg_stim  <- deg_stim %>%
      filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > lfc_cutoff)
    
    g_focal <- unique(deg_focal$gene)
    g_stim  <- unique(deg_stim$gene)
    g_focal <- g_focal[!grepl("^MT-", g_focal)]
    g_stim  <- g_stim[!grepl("^MT-", g_stim)]
    
    shared    <- intersect(g_focal, g_stim)
    focal_only <- setdiff(g_focal, g_stim)
    stim_only  <- setdiff(g_stim,  g_focal)
    
    data.frame(
      cell_type = ct,
      group = c("Focal", "Shared", "Stim"),
      count = c(length(focal_only), length(shared), length(stim_only)),
      stringsAsFactors = FALSE
    )
  })
  
  df <- bind_rows(df_list)
  
  df$cell_type <- factor(df$cell_type, levels = cell_types)
  
  df$group <- factor(df$group, levels = c("Stim", "Shared", "Focal"))
  
  df
}

# ---- make data for UP and DOWN ----
df_up   <- get_ct_counts(direction = "upregulated")
df_down <- get_ct_counts(direction = "downregulated")

plot_stacked_bar <- function(df, title = NULL) {
  ggplot(df, aes(x = cell_type, y = count, fill = group)) +
    geom_col(width = 0.85) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Stim"   = "#4676B5",
      "Shared" = "#BDBDBD",
      "Focal"  = "#E78C8C"
    )) +
    labs(x = NULL, y = "DEG count", title = title) +
    theme_classic(base_size = 8) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      plot.title = element_text(size = 8)
    )
}

p_bar_up   <- plot_stacked_bar(df_up,   title = "Upregulated: Focal / Shared / Stim")
p_bar_down <- plot_stacked_bar(df_down, title = "Downregulated: Focal / Shared / Stim")

p_bar_up+p_bar_down
# ---- save ----
pdf("//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V6/Figure.Focal.Stim.barplot.up.pdf",
    width = 4.2, height = 4.8)
print(p_bar_up)
dev.off()

pdf("//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V6/Figure.Focal.Stim.barplot.down.pdf",
    width = 4.2, height = 4.8)
print(p_bar_down)
dev.off()


##################################################################################################################################
focal_only_up <- setdiff(focal_up, stim_up)
shared_up     <- intersect(focal_up, stim_up)
stim_only_up  <- setdiff(stim_up, focal_up)


get_universe_from_final_results <- function(condition) {
  all <- c()
  for (ct in cell_types) {
    obj <- final_results[[ct]][[condition]]
    
    all <- c(all, obj[["deg"]]$gene)
    
  }
  unique(all[!is.na(all)])
}

universe_focal  <- get_universe_from_final_results("Focal_vs_Distal")
universe_stim   <- get_universe_from_final_results("Stimulated_vs_Distal")
universe_global <- union(universe_focal, universe_stim)


ont_choice  <- "BP"                    # "BP", "CC", or "MF"
min_lfc     <- log2(1.5)
padj_cutoff <- 0.05


## stricter ORA with background/universe + size/overlap filters
min_deg_n   <- 10    # Minimum number of DEGs required to run ORA; skip ORA if fewer than this
min_overlap <- 3     # Minimum overlap (Count) of DEGs per GO term (require Count >= min_overlap)
minGSSize   <- 5    # Minimum GO term size (number of genes in the term) to keep
maxGSSize   <- 2000  # Maximum GO term size to keep (filter out very broad terms)

runGO <- function(genelist,
                  universe,
                  ont = "BP",
                  padj_cutoff = 0.05,
                  min_deg_n = 10,
                  min_overlap = 3,
                  minGSSize = 5,
                  maxGSSize = 2000,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  pAdjustMethod = "BH") {
  
  genes <- unique(genelist)
  genes <- genes[!is.na(genes)]
  genes <- intersect(genes, universe)
  
  if (length(genes) < min_deg_n) return(NULL)
  
  ego <- tryCatch(
    enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      universe      = universe,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff  = pvalueCutoff,
      qvalueCutoff  = qvalueCutoff
    ),
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  df <- ego@result %>%
    filter(p.adjust < padj_cutoff, Count >= min_overlap) %>%
    mutate(log10padj = -log10(p.adjust))
  
  if (nrow(df) == 0) return(NULL)
  df
}


plotGO <- function(ego_df,
                   top_n = 20,
                   title = NULL,
                   fill = "#E41A1C") {
  
  if (is.null(ego_df) || nrow(ego_df) == 0) {
    return(ggplot() + theme_void() + ggtitle(ifelse(is.null(title), "No GO terms", title)))
  }
  
  df <- ego_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n)
  
  ggplot(df, aes(x = log10padj, y = reorder(Description, log10padj))) +
    geom_col(width = 0.8, fill = fill) +
    labs(title = title, x = "-log10(p.adjust)", y = "GO term") +
    theme_classic(base_size = 7) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
      axis.text.x = element_text(size = 7, color = "black"),
      axis.text.y = element_text(size = 6, color = "black")
    )
}



ego_focal_only_up <- runGO(focal_only_up, universe = universe_global, ont = ont_choice,
                           padj_cutoff = padj_cutoff, min_deg_n = min_deg_n,
                           min_overlap = min_overlap, minGSSize = minGSSize, maxGSSize = maxGSSize)

ego_shared_up <- runGO(shared_up, universe = universe_global, ont = ont_choice,
                       padj_cutoff = padj_cutoff, min_deg_n = min_deg_n,
                       min_overlap = min_overlap, minGSSize = minGSSize, maxGSSize = maxGSSize)

ego_stim_only_up <- runGO(stim_only_up, universe = universe_global, ont = ont_choice,
                          padj_cutoff = padj_cutoff, min_deg_n = min_deg_n,
                          min_overlap = min_overlap, minGSSize = minGSSize, maxGSSize = maxGSSize)

p_focal_only_up <- plotGO(ego_focal_only_up, top_n = 10, title = "Focal-only Up", fill = "#E41A1C")
p_shared_up     <- plotGO(ego_shared_up,     top_n = 10, title = "Shared Up",     fill = "#E41A1C")
p_stim_only_up  <- plotGO(ego_stim_only_up,  top_n = 10, title = "Stim-only Up",  fill = "#E41A1C")

(p_focal_only_up | p_shared_up | p_stim_only_up)

######################################################################################################
# One-axis GO plot (stack 3 picked lists top-to-bottom, keep idx order)
######################################################################################################

make_go_pick_df <- function(ego, idx, group_name) {
  if (is.null(ego) || !is.data.frame(ego) || nrow(ego) == 0) return(NULL)
  
  df0 <- ego %>% dplyr::arrange(p.adjust)  
  
  idx <- idx[idx >= 1 & idx <= nrow(df0)]
  if (length(idx) == 0) return(NULL)
  
  df <- df0[idx, , drop = FALSE]          
  
  if (!("log10padj" %in% colnames(df))) {
    df <- df %>% dplyr::mutate(log10padj = -log10(p.adjust))
  }
  
  df$group <- group_name
  df$y_lab <- paste0(group_name, " | ", df$Description)
  df
}

# ---- indices you picked (non-redundant top10) ----
idx_focal_only_up <- c(1, 5, 6, 10, 13, 17, 21, 23, 24, 25)
idx_shared_up     <- c(1, 2, 3, 8, 14, 23, 26, 27, 31, 76)
idx_stim_only_up  <- c(1, 2, 14, 15, 16, 35, 36, 38, 53, 74)

df_focal  <- make_go_pick_df(ego_focal_only_up, idx_focal_only_up, "Focal-only Up")
df_shared <- make_go_pick_df(ego_shared_up,     idx_shared_up,     "Shared Up")
df_stim   <- make_go_pick_df(ego_stim_only_up,  idx_stim_only_up,  "Stim-only Up")

df_all <- dplyr::bind_rows(df_focal, df_shared, df_stim)

y_levels <- c(df_focal$y_lab, df_shared$y_lab, df_stim$y_lab)
df_all$y_lab <- factor(df_all$y_lab, levels = rev(y_levels))  

p_go_oneaxis_up <- ggplot(df_all, aes(x = log10padj, y = y_lab)) +
  geom_col(width = 0.8, fill = "#E41A1C") +
  labs(title = "GO enrichment (Upregulated)", x = "-log10(p.adjust)", y = NULL) +
  theme_classic(base_size = 7) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 8),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 6, color = "black")
  )

p_go_oneaxis_up

pdf("//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V6/Figure.Integrated.DEG.GO.pdf",
    width = 4, height = 4)
p_go_oneaxis_up
dev.off()







##################################################################################################################################
# Barplot: per cell type stacked counts (Focal / Shared / Stim) %
##################################################################################################################################
cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L5ET", "Exc_L5IT_GRIN3A",
  "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", 
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", 
  "Inh_LAMP5", "Inh_LAMP5_LHX6", "Astrocyte", "Microglia", 
  "Olg", "OPC", "Macrophage", "NK/T", "Vascular"
)
cell_types <- rev(cell_types)

get_ct_counts <- function(condition_focal = "Focal_vs_Distal",
                          condition_stim  = "Stimulated_vs_Distal",
                          direction = c("upregulated", "downregulated"),
                          padj_cutoff = 0.05,
                          lfc_cutoff  = log2(1.5)) {
  direction <- match.arg(direction)
  
  df_list <- lapply(cell_types, function(ct) {
    deg_focal <- final_results[[ct]][[condition_focal]][[direction]]
    deg_stim  <- final_results[[ct]][[condition_stim ]][[direction]]
    
    # threshold filter
    deg_focal <- deg_focal %>%
      dplyr::filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > lfc_cutoff)
    deg_stim  <- deg_stim %>%
      dplyr::filter(p_val_adj < padj_cutoff, abs(avg_log2FC) > lfc_cutoff)
    
    g_focal <- deg_focal$gene
    g_stim  <- deg_stim$gene
    
    # remove MT-
    g_focal <- g_focal[!grepl("^MT-", g_focal)]
    g_stim  <- g_stim[!grepl("^MT-", g_stim)]
    
    # unique just in case
    g_focal <- unique(g_focal)
    g_stim  <- unique(g_stim)
    
    shared <- intersect(g_focal, g_stim)
    focal_only <- setdiff(g_focal, g_stim)
    stim_only  <- setdiff(g_stim,  g_focal)
    
    counts <- c(length(focal_only), length(shared), length(stim_only))
    total <- sum(counts)
    percent <- if (total == 0) rep(0, 3) else counts / total * 100
    
    data.frame(
      cell_type = ct,
      group = c("Focal", "Shared", "Stim"),
      count = counts,
      percent = percent,
      stringsAsFactors = FALSE
    )
  })
  
  df <- dplyr::bind_rows(df_list)
  
  # keep your stacked order
  df$group <- factor(df$group, levels = c("Stim", "Focal", "Shared"))
  
  df$cell_type <- factor(df$cell_type, levels = cell_types)
  
  df
}


# ---- make data for UP and DOWN ----
# ---- make data for UP and DOWN ----
df_up   <- get_ct_counts(direction = "upregulated")
df_down <- get_ct_counts(direction = "downregulated")

plot_stacked_bar <- function(df, title = NULL) {
  ggplot(df, aes(x = cell_type, y = percent, fill = group)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = c(
      "Stim"   = "#4676B5",
      "Shared" = "#BDBDBD",
      "Focal"  = "#E78C8C"
    )) +
    scale_y_continuous(
      expand = c(0, 0),
      labels = function(x) paste0(x, "%")
    ) +
    coord_cartesian(ylim = c(0, 100)) + 
    coord_flip() +                    
    labs(x = NULL, y = "Percent of DEGs (%)", title = title) +
    theme_classic(base_size = 8) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      plot.title  = element_text(size = 8)
    )
}

p_bar_up   <- plot_stacked_bar(df_up,   title = "Upregulated: Focal / Shared / Stim")
p_bar_down <- plot_stacked_bar(df_down, title = "Downregulated: Focal / Shared / Stim")

p_bar_up + p_bar_down


pdf("//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V6/Figure.Focal.Stim.barplot.pdf",
    width = 10, height = 6)
p_bar_up + p_bar_down
dev.off()
