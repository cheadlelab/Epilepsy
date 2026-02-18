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
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
save_path <- "C:\\Users\\thech\\Desktop"

load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))

# Prepare the Seurat object
test <- JoinLayers(test)
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
test <- subset(test, subset = Annotation.fine != "Other")
DefaultAssay(test) <- "RNA"

GABAergic <- c("Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier")
Glu <- c("Exc_L234IT_ERBB4",  # L2-3-4
         "Exc_L23IT",         # L2/3
         "Exc_L4IT",          # L4
         "Exc_L4IT_PLCH1",    # L4
         "Exc_L5ET",          # L5
         "Exc_L5IT_GRIN3A",    # L5
         "Exc_L5IT_GRIN3A-",   # L5
         "Exc_L56NP",         # L5-6
         "Exc_L56IT_CAR3",    # L5-6
         "Exc_L6IT",          # L6
         "Exc_L6CT",          # L6
         "Exc_L6b")           # L6
glia <- c("Astrocyte", "OPC", "Olg", "Microglia", "Vascular")
other <- c( "Macrophage", "NK/T")

ordered_levels <- c(GABAergic, Glu, glia, other)
test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)

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


save(all.markers, file = paste0(file_path, "/celltype.markers.RData"))

default_colors <- hue_pal()(26)
names(default_colors) <- levels(test@meta.data$Annotation.fine)


DEGs <- c(
  "SST",      # Inh_SST
  "KCNS3",    # Inh_PVALB
  "VIP",      # Inh_VIP
  "CXCL14",   # Inh_CXCL14
  "LAMP5",    # Inh_LAMP5
  "LHX6",     # Inh_LAMP5_LHX6
  "UNC5B",    # Chandelier
  
  "SLC17A7",  # Exc
  "ERBB4",    # Exc_L234IT_ERBB4
  "CUX2",     # Exc_L2/3
  "RORB",     # Exc_L4
  "PLCH1",    # Exc_L4IT_PLCH1
  "BCL11B",   # Exc_L5ET
  "FEZF2",    # Exc_L5
  "GRIN3A",   # Exc_L5IT_GRIN3A
  "TLE4",     # Exc_L56NP
  "TSHZ2",    # Exc_L5/6
  "CA3",      # Exc_L56IT_CAR3
  "FOXP2",    # Exc_L6
  "SULF1",    # Exc_L6CT
  "NPFFR2",   # Exc_L6b
  
  "OBI1-AS1", # Astrocyte
  "PDGFRA",   # OPC
  "MOG",      # Olg
  "CX3CR1",   # Microglia
  "EBF1",     # Vascular
  "CD163",    # Macrophage
  "PTPRC"   # NK/T
  
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
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(20, "mm"),
      barheight = unit(2, "mm"),
      frame.colour = "black"
    )
  ) +
  scale_size(
    range = c(0, 2.5),
    breaks = c(20, 40, 60, 80),
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom"
    )
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
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.spacing.x = unit(0.5, "cm")
  ) +
  coord_flip()

p



pdf(file = paste0(save_path, "/Figure1G.pdf"), width = 4.2, height = 3.8)
print(p)
dev.off()
