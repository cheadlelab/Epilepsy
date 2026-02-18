rm(list = ls(all = TRUE))

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
# file_path <- "/grid/Epilepsy/output"
# load(paste0(file_path, "/step2.1.AllComparisons.RData"))
load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))

cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", "Exc_L5IT_GRIN3A",
  "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", 
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", 
  "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier", "Astrocyte", "Microglia", 
  "Olg", "OPC", "Macrophage", "NK/T", "Vascular"
)


celltype <- "Inh_PVALB"
# ===============================
# Stimulated_vs_Distal
# ===============================
stim.markers <- final_results[[celltype]][["Stimulated_vs_Distal"]][["deg"]] 
genes.up.stim   <- stim.markers$gene[stim.markers$avg_log2FC >  log2(1.5) & stim.markers$p_val_adj < 0.05 & !grepl("^MT-", stim.markers$gene)]
genes.down.stim <- stim.markers$gene[stim.markers$avg_log2FC < -log2(1.5) & stim.markers$p_val_adj < 0.05 & !grepl("^MT-", stim.markers$gene)]
# ===============================
# Focal_vs_Distal
# ===============================
focal.markers <- final_results[[celltype]][["Focal_vs_Distal"]][["deg"]] 
genes.up.focal  <- focal.markers$gene[focal.markers$avg_log2FC >  log2(1.5) & focal.markers$p_val_adj < 0.05 & !grepl("^MT-", focal.markers$gene)]
genes.down.focal<- focal.markers$gene[focal.markers$avg_log2FC < -log2(1.5) & focal.markers$p_val_adj < 0.05 & !grepl("^MT-", focal.markers$gene)]


df1 <- stim.markers %>% 
  dplyr::select(gene, avg_log2FC, p_val_adj) %>%
  dplyr::rename(log2FC_stim = avg_log2FC, p_adj_stim = p_val_adj)

df2 <- focal.markers %>% 
  dplyr::select(gene, avg_log2FC, p_val_adj) %>%
  dplyr::rename(log2FC_focal = avg_log2FC, p_adj_focal = p_val_adj)

df.merge <- inner_join(df1, df2, by = "gene")


# Keep only genes significant in both contrasts (no log2FC filtering)
df.plot <- df.merge %>% dplyr::filter(p_adj_stim < 0.05, p_adj_focal < 0.05)

# Classification by thresholds for coloring
th <- log2(1.5)
df.plot <- df.plot %>%
  dplyr::mutate(
    color_class = dplyr::case_when(
      log2FC_stim >=  th & log2FC_focal >=  th ~ "both_up_thr",
      log2FC_stim <= -th & log2FC_focal <= -th ~ "both_down_thr",
      TRUE ~ "other"
    ),
    rank_metric = abs(log2FC_stim) * abs(log2FC_focal)
  )

# Top10 labels within each group
df.top_up <- df.plot %>%
  dplyr::filter(color_class == "both_up_thr") %>%
  dplyr::arrange(dplyr::desc(rank_metric)) %>%
  dplyr::slice(1:min(10, n()))

df.top_down <- df.plot %>%
  dplyr::filter(color_class == "both_down_thr") %>%
  dplyr::arrange(dplyr::desc(rank_metric)) %>%
  dplyr::slice(1:min(10, n()))

suppressPackageStartupMessages(require(ggrepel))

p <- ggplot(df.plot, aes(x = log2FC_stim, y = log2FC_focal)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_hline(yintercept =  th, linetype = "dotted", linewidth = 0.3, color = "grey60") +
  geom_hline(yintercept = -th, linetype = "dotted", linewidth = 0.3, color = "grey60") +
  geom_vline(xintercept =  th, linetype = "dotted", linewidth = 0.3, color = "grey60") +
  geom_vline(xintercept = -th, linetype = "dotted", linewidth = 0.3, color = "grey60") +
  geom_point(aes(color = color_class), alpha = 0.75, size = 1.1) +
  ggrepel::geom_text_repel(
    data = df.top_up,
    aes(label = gene),
    size = 3, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0
  ) +
  ggrepel::geom_text_repel(
    data = df.top_down,
    aes(label = gene),
    size = 3, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      "both_up_thr" = "#D62728",   # red
      "both_down_thr" = "#1F77B4", # blue
      "other" = "grey65"           # gray
    ),
    breaks = c("both_up_thr", "both_down_thr", "other"),
    labels = c("double positive", "double negative", "other")
  ) +
  labs(
    title = paste0(celltype, " — log2FC: Stimulated_vs_Distal (x) vs Focal_vs_Distal (y)"),
    x = "Stimulated_vs_Distal (log2FC)",
    y = "Focal_vs_Distal (log2FC)",
    color = "class",
    caption = "Significance: both p_adj < 0.05. Dotted lines show ±log2(1.5). Gray = below threshold in at least one axis."
  ) +
  theme_bw(base_size = 12)

print(p)

# Optional: save
# plot_dir <- file.path(file_path, "plots_stim_vs_focal_scatter")
# if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
# ggsave(file.path(plot_dir, paste0("scatter_", gsub("[^A-Za-z0-9_]+","_",celltype), ".pdf")),
#        p, width = 4.8, height = 4.8, units = "in")



# Panel plot for a selected subset of cell types

selected_cells <- c(
  "Exc_L23IT", "Exc_L4IT", "Exc_L5ET", "Exc_L5IT_GRIN3A",
  "Exc_L6CT", "Exc_L6IT", "Exc_L6b",
  "Inh_SST", "Inh_PVALB", "Inh_VIP"
)

selected_cells <- c("Exc_L23IT", "Exc_L5IT_GRIN3A", "Exc_L6IT", "Inh_SST")

th <- log2(1.5)

build_cell_df <- function(ct) {
  stim <- final_results[[ct]][["Stimulated_vs_Distal"]][["deg"]]
  focal <- final_results[[ct]][["Focal_vs_Distal"]][["deg"]]
  if (is.null(stim) || is.null(focal)) return(NULL)
  d1 <- dplyr::select(stim, gene, avg_log2FC, p_val_adj) %>%
    dplyr::rename(log2FC_stim = avg_log2FC, p_adj_stim = p_val_adj) %>%
    dplyr::filter(!grepl("^MT-", gene))
  
  d2 <- dplyr::select(focal, gene, avg_log2FC, p_val_adj) %>%
    dplyr::rename(log2FC_focal = avg_log2FC, p_adj_focal = p_val_adj) %>%
    dplyr::filter(!grepl("^MT-", gene))
  
  m <- dplyr::inner_join(d1, d2, by = "gene")
  m <- dplyr::filter(m, p_adj_stim < 0.05, p_adj_focal < 0.05)
  if (!nrow(m)) return(NULL)
  m %>%
    dplyr::mutate(
      CellType = ct,
      color_class = dplyr::case_when(
        log2FC_stim >=  th & log2FC_focal >=  th ~ "both_up",
        log2FC_stim <= -th & log2FC_focal <= -th ~ "both_down",
        TRUE ~ "other"
      ),
      rank_metric = abs(log2FC_stim) * abs(log2FC_focal)
    )
}

panel_df_list <- lapply(selected_cells, build_cell_df)
panel_df <- dplyr::bind_rows(panel_df_list)

top_up_df <- panel_df %>%
  dplyr::filter(color_class == "both_up") %>%
  dplyr::group_by(CellType) %>%
  dplyr::slice_max(order_by = rank_metric, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

top_down_df <- panel_df %>%
  dplyr::filter(color_class == "both_down") %>%
  dplyr::group_by(CellType) %>%
  dplyr::slice_max(order_by = rank_metric, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

# --- make gene labels italic (for ggrepel) ---
top_up_df$gene_italic   <- paste0("italic('", top_up_df$gene, "')")
top_down_df$gene_italic <- paste0("italic('", top_down_df$gene, "')")

p_panel <- ggplot(panel_df, aes(x = log2FC_stim, y = log2FC_focal)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, color = "grey60") +
  geom_hline(yintercept =  th, linetype = "dotted", linewidth = 0.25, color = "grey70") +
  geom_hline(yintercept = -th, linetype = "dotted", linewidth = 0.25, color = "grey70") +
  geom_vline(xintercept =  th, linetype = "dotted", linewidth = 0.25, color = "grey70") +
  geom_vline(xintercept = -th, linetype = "dotted", linewidth = 0.25, color = "grey70") +
  geom_point(data = subset(panel_df, color_class == "other"),
             color = "grey70", alpha = 0.6, size = 0.8) +
  geom_point(data = subset(panel_df, color_class != "other"),
             aes(color = color_class), alpha = 0.8, size = 1.0) +
  ggrepel::geom_text_repel(
    data = top_up_df, aes(label = gene_italic),
    size = 2.6, color = "#D62728", max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.2, min.segment.length = 0,
    parse = TRUE
  ) +
  ggrepel::geom_text_repel(
    data = top_down_df, aes(label = gene_italic),
    size = 2.6, color = "#1F77B4", max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.2, min.segment.length = 0,
    parse = TRUE
  ) +
  scale_color_manual(
    values = c("both_up" = "#D62728", "both_down" = "#1F77B4"),
    breaks = c("both_up", "both_down"),
    labels = c("double positive", "double negative"),
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  facet_wrap(~ CellType, scales = "free", nrow = 1) +
  labs(
    title = "Stimulated_vs_Distal (x) vs Focal_vs_Distal (y)",
    x = "Stimulated_vs_Distal (log2FC)",
    y = "Focal_vs_Distal (log2FC)",
    color = "class"
  ) +
  theme_bw(base_size = 10)


pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/Figure.overlap.pdf", width = 9, height = 2.8)
print(p_panel)
dev.off() # Close the PDF device
