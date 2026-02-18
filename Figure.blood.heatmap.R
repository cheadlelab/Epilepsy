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
  "CD8+_T_1",            # CD8+ T cell subset 1 Effector
  "CD8+_T_2",            # CD8+ T cell subset 2 Memory
  "CD8+_T_3",            # CD8+ T cell subset 3 Exhausted
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







Idents(test) <- test@meta.data$Annotation.fine
all.markers <- FindAllMarkers(object = test, assay = "RNA",
                              logfc.threshold = log2(1.5),
                              min.pct = 0.1
)
save(all.markers, file = paste0(file_path, "/celltype.markers.RData"))

top10_genes <- all.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

default_colors <- hue_pal()(20)
names(default_colors) <- levels(test@meta.data$Annotation.fine)


FeaturePlot(test, features = c("CD79A"))

DEGs <- c(
  "FCGR3A",   # FCGR3A+ Mono
  "CD163",    # Macrophage
  "FCER1A",     # DC
  "CLEC9A",   # cDCs
  "LILRA4",   # pDC
  "CPA3",     # Mast
  
  "PTPRC",    # immune
  "CD3D",     # T
  "CCR7",     # Naive_CD4+_T
  "IL7R",     # Memory_CD4+_T
  "FOXP3",     # regular_T
  "GZMH",     # CD8+_T_1
  "LEF1",     # CD8+_T_2
  "ME1",   # CD8+_T_3
  "TRDC",     # gamma_delta_T
  
  "KLRC3",    # NK_1
  "CD160",   # NK_2
  "XCL1",    # NK_3
  
  "DNTT",   # Pre_B
  "MS4A1",     # B
  "CD79A",     # Plasma_B
  
  "MKI67"    # MKI67+_NKT
)

dotplot_data <- DotPlot(
  object = test,
  features = DEGs,
  assay = "RNA",
  cluster.idents = FALSE
)$data

dotplot_data$id <- factor(dotplot_data$id, levels = ordered_levels)
dotplot_data$features.plot <- factor(dotplot_data$features.plot, levels = rev(DEGs))

# Calculate expression limits from data (rounded to nearest integer)
exp_min <- floor(min(dotplot_data$avg.exp.scaled))
exp_max <- ceiling(max(dotplot_data$avg.exp.scaled))
exp_limits <- c(exp_min, exp_max)

# Generate bubble plot
p <- ggplot(dotplot_data, 
            aes(x = features.plot, y = id)) +       # X: Gene symbols, Y: Cell clusters
  
  geom_point(
    aes(size = pct.exp,                    # Point size = detection percentage
        color = avg.exp.scaled),           # Point color = Z-scored expression
    shape = 19,                            # Solid circles
    stroke = 0.5                           # Point border thickness
  ) +
  
  scale_color_gradientn(
    colors = c("#f0f0f0", "#6B8CAC", "#0D2B52"),  # White → Medium blue → Navy (3-step gradient)
    values = scales::rescale(seq(exp_limits[1], exp_limits[2], length = 3)), # Linear scaling
    limits = exp_limits,                            # Data-driven value boundaries
    guide = guide_colorbar(
      direction = "horizontal",                     # Horizontal legend orientation
      title.position = "top",                       # Title above color bar
      barwidth = unit(20, "mm"),                    # Color bar length (match figure e/f)
      barheight = unit(2, "mm"),                    # Color bar thickness
      frame.colour = "black"                        # Legend border
    )
  ) +
  scale_size(
    range = c(0.5, 7),                                 # Minimum and maximum point diameters
    breaks = c(20, 40, 60, 80),                     # Standardized percentage intervals
    guide = guide_legend(
      direction = "horizontal",                     # Horizontal legend layout
      title.position = "top",                       # Title above size legend
      label.position = "bottom"                     # Percentage labels below points
    )
  ) +
  theme_classic(base_size = 11) +                    # Clean theme with 11pt base font
  labs(x = NULL, y = NULL) +                        # Remove default axis labels
  theme(
    # Gene label styling (X-axis)
    axis.text.x = element_text(
      angle = 50,                                   # 50-degree text rotation
      hjust = 1,                                    # Right alignment
      vjust = 1,                                    # Top vertical positioning
      face = "plain",                              # Italic gene symbols
      color = "black",                              # High-contrast text
      size = rel(1.2)                              # 5% larger than base
    ),
    # Cluster label styling (Y-axis)
    axis.text.y = element_text(
      face = "plain",                                # Bold cluster names
      color = "black",                           # Dark gray for readability
      size = rel(1.2)                              # 5% larger than base
    ),
    # Panel border customization
    panel.border = element_rect(                    # Solid black border
      colour = "black",
      fill = NA,
      linewidth = 0.8                               # 0.8pt line thickness
    ),
    # Legend panel configuration
    legend.position = "bottom",                     # Unified legend at bottom
    legend.box = "horizontal",                      # Horizontal legend stacking
    legend.direction = "horizontal",               # Left-to-right layout
    legend.spacing.x = unit(0.5, "cm")              # 0.5cm spacing between legends
    
  ) +
  coord_flip()                                      # Swap axes for vertical gene display
p


# 8 inches x 8.5 inches
ggsave(filename = paste0(file_path, "/celltype_bubbleplot.pdf"),
       plot = p,
       width = 12,
       height = 8,
       dpi = 300)







