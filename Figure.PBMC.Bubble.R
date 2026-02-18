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
save_path <- "//grid/cheadle_home/qianyu/Epilepsy/figures/manuscript"

load(paste0(file_path, "/step1.2.umap.fine.RData"))

# Prepare the Seurat object
`%!in%` <- Negate(`%in%`)
test <- subset(test, subset = Annotation.fine %!in% c("Myeloid/platelets", "Myeloid/NKT", "RBC", "B/Myeloid", "B/NKT", "B/cDC2", "B/platelets", "B/other"))
# test$Annotation.fine <- droplevels(test$Annotation.fine)

unique(test@meta.data$Annotation.fine)
unique(test@meta.data$patient)

patient_mapping <- c(
  "Epileptic_101923" = "patient1",
  "Epileptic_101323" = "patient2",
  "Epileptic_22924" = "patient3",
  "Epileptic_2849" = "patient4"
)

# Add the new 'patient_id' column to metadata
test@meta.data$patient_id <- ifelse(
  test@meta.data$patient %in% names(patient_mapping),
  patient_mapping[test@meta.data$patient],
  as.character(test@meta.data$patient)  # Keep original ID if not in mapping
)
unique(test@meta.data$patient_id)

DefaultAssay(test) <- "RNA"
test@active.ident <- as.factor(test@meta.data$orig.ident)

ordered_levels <- c(
  ## Myeloid lineage
  "CD14_Mono",      # Classical monocytes
  "CD16_Mono",      # Non-classical monocytes
  "cDC1",           # Conventional DC subset 1
  "cDC2",           # Conventional DC subset 2
  "pDC",            # Plasmacytoid DCs
  "Mast",           # Mast cells
  
  ## Lymphoid lineage
  ### B-cell series (progenitor → mature → effector)
  "pro_B",          
  "B_naive",        
  "B_intermediate", 
  "B_memory",       
  "Plasmablast",    
  
  ### T-cell series
  "CD4_Naive",      # Naïve CD4+
  "CD4_Memory",     # Memory CD4+
  "Treg",           # Regulatory T cells
  "CD4_CTL",        # CD4+ cytotoxic
  "CD8_Naive",      # Naïve CD8+
  "CD8_Memory",     # Memory CD8+
  "MAIT",           # Mucosal-associated invariant T
  "gdT",            # γδ T cells
  
  ### NK-cell series
  "NK_CD56bright",  # CD56bright NK
  "NK",             # CD56dim/adaptive NK
  
  ## Proliferating populations
  "MKI67+_NKT",     # Cycling NKT
  
  ## Erythroid & platelet lineages
  "EMP",            # Erythro-myeloid progenitors
  "Late_Eryth",     # Late erythroid cells
  "Platelets"       # Platelets
)

test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)
print(levels(test@meta.data$Annotation.fine))



Idents(test) <- test@meta.data$Annotation.fine
all.markers <- FindAllMarkers(object = test, assay = "RNA",
                              logfc.threshold = log2(1.5),
                              min.pct = 0.1
)




top10_genes <- all.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

default_colors <- hue_pal()(25)
names(default_colors) <- levels(test@meta.data$Annotation.fine)



DEGs <- c(
  "CD14",     # CD14_Mono
  "FCGR3A",   # CD16_Mono
  "CLEC9A",   # cDC1
  "FCER1A",   # cDC2
  "CLEC4C",   # pDC
  "HDC",      # Mast
  
  "CD79A", "MS4A1", # B
  "DNTT",     # pro_B
  "TCL1A",    # B_naive
  "SOX5",     # B_intermediate
  "SSPN",     # B_memory
  "IGHA1",     # Plasmablast
  
  "CD3D",     # T
  "LEF1",     # CD4_Naive
  "IL7R",     # CD4_Memory
  "FOXP3",    # Treg
  "GZMH",     # CD4_CTL
  "CD8A",     # CD8
  "CCR7",     # CD8_Naive
  "TRGC2",     # CD8_Memory
  "SLC4A10",  # MAIT
  "KLRG1",    # γδ T cells
  
  "NCAM1",    # NK_CD56bright
  "FCER1G",   # NK
  
  "MKI67",    # MKI67+_NKT
  
  "MYCT1",    # EMP
  "CTSE",     # Late_Eryth
  "PPBP"      # Platelets
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
            aes(x = features.plot, y = id)) +
  
  geom_point(
    aes(size = pct.exp,
        color = avg.exp.scaled),
    shape = 19,
    stroke = 0.5
  ) +
  
  scale_color_gradientn(
    colors = c("#2166AC", "#F7F7F7", "#B2182B"),
    values = scales::rescale(seq(exp_limits[1], exp_limits[2], length = 3)),
    limits = exp_limits,
    guide = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      barwidth = unit(2, "mm"),
      barheight = unit(20, "mm"),
      frame.colour = "black"
    )
  ) +
  scale_size(
    range = c(0, 2),
    breaks = c(20, 40, 60, 80)
  ) +
  theme_classic(base_size = 7) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      face = "plain",
      color = "black",
      size = 6
    ),
    axis.text.y = element_text(
      face = "plain",
      color = "black",
      size = 6
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.position = "right",
    # legend.box = "vertical", #horizontal
    # legend.direction = "vertical", #horizontal
    legend.spacing.x = unit(0.5, "cm")
  ) +
  coord_flip()

p



pdf(file = paste0(save_path, "/S.Figure2.PMBC.Bubble.pdf"), width = 3.6, height = 3.2)
print(p)
dev.off()
