# MOSS root+scvi envs
rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)

library(reticulate)
use_condaenv("/home/qianyu/anaconda3/envs/scvi-env", required = TRUE)

##################################################################################################################################
###################################################### load Cheadle Dataset ######################################################
##################################################################################################################################
file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output"

load(paste0(file_path, "/step2.Merge.relabel.RData"))


# p1 <- DimPlot(test, reduction = "umap_harmony", group.by = "Annotation2", raster=FALSE, label = TRUE)
# p2 <- DimPlot(test, reduction = "umap_harmony", group.by = "clusters_harmony", raster=FALSE, label = TRUE)
# p3 <- DimPlot(test, reduction = "umap_harmony", group.by = "condition", raster=FALSE, label = TRUE)
# p4 <- DimPlot(test, reduction = "umap_harmony", group.by = "patient", raster=FALSE, label = TRUE)
# p1+p2+p3+p4


test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
tt.tmp <- subset(test, subset = Annotation2 == "Microglia")

orig_ident_table <- table(tt.tmp@meta.data[["orig.ident"]])
selected_orig_ident <- names(orig_ident_table[orig_ident_table > 10])
tt.tmp <- subset(tt.tmp, subset = orig.ident %in% selected_orig_ident)

tt.tmp[["RNA"]] <- as(object = tt.tmp[["RNA"]], Class = "Assay5")
tt.tmp[["RNA"]] <- split(tt.tmp[["RNA"]], f = tt.tmp$orig.ident)
cheadle <- tt.tmp
##################################################################################################################################
####################################################### load Allen Dataset #######################################################
##################################################################################################################################

load(paste0(file_path, "/step0.AllenCortical.RData"))
test@meta.data$dataset <- "Allen_Cortical"
test@meta.data$sex <- test@meta.data$donor_sex_label
test@meta.data$donor_sex_label <- NULL
test@meta.data$condition <- "Control"
Allen1 <- test
# unique(Allen1@meta.data$cluster_label_simplify)
Allen1 <- subset(Allen1, subset = cluster_label_simplify %in% c("Micro L1-6"))
a <- rownames(Allen1)
a <- gsub("\\.AS", "-AS", a)
a <- gsub("\\.OT", "-OT", a)
a <- gsub("HLA\\.", "HLA-", a)
rownames(Allen1) <- a
# obj <- subset(obj, features = genes)

load(paste0(file_path, "/step0.AllenMTG.RData"))
test@meta.data$dataset <- "Allen_MTG"
test@meta.data$sex <- test@meta.data$donor_sex_label
test@meta.data$donor_sex_label <- NULL
test@meta.data$cluster_label_simplify <- test@meta.data$cluster_label
test@meta.data$condition <- "Control"
Allen2 <- test
# unique(Allen2@meta.data$subclass_label)
Allen2 <- subset(Allen2, subset = subclass_label %in% c("Microglia-PVM"))
a <- rownames(Allen2)
a <- gsub("\\.AS", "-AS", a)
a <- gsub("\\.OT", "-OT", a)
a <- gsub("HLA\\.", "HLA-", a)
rownames(Allen2) <- a
# obj <- subset(obj, features = genes)

###################################################################################################################################
############################################################## merge ##############################################################
###################################################################################################################################
tt.tmp <- merge(cheadle, y = c(Allen1, Allen2))


tt.tmp <- NormalizeData(tt.tmp)
tt.tmp <- FindVariableFeatures(tt.tmp)
tt.tmp <- ScaleData(tt.tmp)
tt.tmp <- RunPCA(tt.tmp)


tt.tmp <- IntegrateLayers(
  object = tt.tmp, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


ElbowPlot(tt.tmp)

dims_parameter <- 15
tt.tmp <- RunUMAP(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindNeighbors(tt.tmp, reduction = "harmony", dims = 1:dims_parameter) #pca
tt.tmp <- FindClusters(tt.tmp, resolution = 0.05)

p1 <- DimPlot(tt.tmp, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
p2 <- DimPlot(tt.tmp, reduction = "umap", group.by = "condition", label = TRUE, raster = FALSE)
p1+p2


tt.tmp <- JoinLayers(tt.tmp)
DefaultAssay(tt.tmp) <- "RNA"
Idents(tt.tmp) <- 'seurat_clusters'
allmarkers <- FindAllMarkers(object = tt.tmp, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
seurat <- allmarkers$gene

####################################################################################################################################
############################################################## Xenium ##############################################################
####################################################################################################################################

save.path <- "/home/qianyu/Desktop/grid/Epilepsy/Xenium/output/"
load(paste0(save.path, "step2.umap.RData"))

category_mapping <- c(
  "1C SA" = "Distal",
  "1E" = "Focal",
  "2C" = "Distal",
  "2E" = "Focal",
  "4C" = "Distal",
  "4E" = "Focal"
)

a <- rownames(test@assays$RNA)


# Apply the mapping to create the new broad category column
test@meta.data$condition <- category_mapping[test@meta.data$orig.ident]

DimPlot(test, reduction = "umap_normal", group.by = "condition", raster=FALSE, label = TRUE)

t.marker <- subset(allmarkers, subset = gene %in% a)
overlap <- intersect(a, seurat)


xenium.microglia <- subset(test, subset = predicted.celltype == "Microglia")



xenium.microglia <- NormalizeData(xenium.microglia)
xenium.microglia <- FindVariableFeatures(xenium.microglia)
xenium.microglia <- ScaleData(xenium.microglia)
xenium.microglia <- RunPCA(xenium.microglia)

ElbowPlot(xenium.microglia)


dims_parameter <- 10
xenium.microglia <- RunUMAP(xenium.microglia, reduction = "pca", dims = 1:dims_parameter, reduction.key = "normal_")
xenium.microglia <- FindNeighbors(xenium.microglia, reduction = "pca", dims = 1:dims_parameter)
xenium.microglia <- FindClusters(xenium.microglia, resolution = 0.1)
xenium.microglia@reductions$umap_normal <- xenium.microglia@reductions$umap
xenium.microglia@meta.data$clusters_normal <- xenium.microglia@meta.data$seurat_clusters

DimPlot(xenium.microglia, reduction = "umap_normal", group.by = "clusters_normal", raster=FALSE, label = TRUE)+NoLegend()

VlnPlot(xenium.microglia, features = c("CD163"), group.by = "clusters_normal", pt.size=0)
VlnPlot(xenium.microglia, features = c("TREM2", "CX3CR1", "P2RY13", "CD68", "CD14", "ITGAM"), group.by = "clusters_normal", pt.size=0)


# in CD163+ microglia, CD163+ expression realated to the layer distribution
tt.tmp <- subset(xenium.microglia, subset = clusters_normal == 2)

Idents(tt.tmp) <- "layer"
deg_results <- FindMarkers(tt.tmp, 
                           ident.1 = "Layer-Exc23", 
                           ident.2 = "Layer-Exc456", 
                           features = "CD163", # or a list of genes you're interested in
                           group.by = "layer", 
                           logfc.threshold = 0, # Set log fold change threshold to 0 to get all genes
                           test.use = "wilcox") # Wilcoxon rank sum test is default
Idents(tt.tmp) <- "condition"
deg_results <- FindMarkers(tt.tmp, ident.1 = "Focal", ident.2 = "Distal", 
                           only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0)






table(xenium.microglia@meta.data$clusters_normal, xenium.microglia@meta.data$layer)

genes <- t.marker$gene[t.marker$avg_log2FC > 0 & t.marker$cluster == 0]
genes <- list(genes)
xenium.microglia <- AddModuleScore(
  object = xenium.microglia,
  features = genes,
  ctrl = 5,
  name = 'mg'
)

DimPlot(xenium.microglia, reduction = "umap_normal", group.by = "clusters_normal", raster=FALSE, label = TRUE)+NoLegend()
DimPlot(xenium.microglia, reduction = "umap_normal", group.by = "condition", raster=FALSE, label = TRUE)

FeaturePlot(xenium.microglia, features = c("CD163"), raster=FALSE, label = TRUE)


VlnPlot(xenium.microglia, features = c("mg1"), group.by = "clusters_normal")
VlnPlot(xenium.microglia, features = c("CD163"), group.by = "clusters_normal")




library(ggplot2)
ggplot(xenium.microglia@meta.data, aes(x=seurat_clusters, fill=condition)) + geom_bar()


#####################################################################################################################################
############################################################## Figures ##############################################################
#####################################################################################################################################


pdf(file = paste0("/home/qianyu/Desktop", "/Mg.umap.snRNA.pdf"), width = 10, height = 10)
FeaturePlot(tt.tmp, features = "CD163", label = FALSE, raster = FALSE, pt.size = 0.5) + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device



pdf(file = paste0("/home/qianyu/Desktop", "/Mg.umap.Xenium.pdf"), width = 10, height = 10)
FeaturePlot(xenium.microglia, features = "CD163", label = FALSE, raster = FALSE, pt.size = 0.5) + 
  NoLegend() + ggtitle(NULL) + NoAxes()
dev.off() # Close the PDF device


# 5 inch x 5 inch
library(ggplot2)
library(dplyr)

# Create a data frame for cluster 2 (row 2)
cluster2_data <- data.frame(
  Layer = factor(c("Layer-NVU", "Layer-Exc23", "Layer-Exc456", "Layer-Olg/OPC"), 
                 levels = c("Layer-NVU", "Layer-Exc23", "Layer-Exc456", "Layer-Olg/OPC")),
  Count = c(1032, 477, 286, 268)
)

# Calculate proportions and add them to the data frame
cluster2_data <- cluster2_data %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1))

# Create a pie chart for cluster 2 with proportion labels
ggplot(cluster2_data, aes(x = "", y = Count, fill = Layer)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5)) +
  labs(title = "Layer Distribution in Cluster 2",
       fill = "Layer")


#####################################################################################################################################
############################################################### CD163 ###############################################################
#####################################################################################################################################

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scales)
library(ggsignif)

# Assuming tt.tmp is a Seurat object and 'CD163' is the gene of interest
# Extract data for plotting
gene_data <- FetchData(tt.tmp, vars = c("CD163", "layer"))

# Check the structure of gene_data to make sure it has the correct columns
str(gene_data)

# Order the 'layer' variable according to the specified order
gene_data$layer <- factor(gene_data$layer, levels = c("Layer-NVU", "Layer-Exc23", "Layer-Exc456", "Layer-Olg/OPC"))

# Get the default ggplot2 color palette for the first four categories
default_colors <- hue_pal()(4)

# Creating the boxplot with ggplot2 with the specified order and applying the default colors
p <- ggplot(gene_data, aes(x = layer, y = CD163, fill = layer)) +
  geom_boxplot() +
  scale_fill_manual(values = default_colors) +
  labs(title = "CD163 Expression in CD163+ Microglia Across Layers", x = "Layer", y = "Expression Level") +
  theme_minimal()  # Adds a minimalistic theme to the plot

# Add pairwise comparisons and annotations for significance levels with adjusted heights
# p <- p + geom_signif(
#   comparisons = list(c("Layer-NVU", "Layer-Exc23"), c("Layer-NVU", "Layer-Exc456"), 
#                      c("Layer-NVU", "Layer-Olg/OPC"), c("Layer-Exc23", "Layer-Exc456"), 
#                      c("Layer-Exc23", "Layer-Olg/OPC"), c("Layer-Exc456", "Layer-Olg/OPC")),
#   map_signif_level = TRUE,
#   test = "wilcox.test",
#   test.args = list(p.adjust.method = "BH"),
#   y_position = c(9, 9.5, 10, 7.5, 8.5, 8)
# )
# Define the pairs of layers you want to compare
comparisons_list <- list(
  c("Layer-NVU", "Layer-Exc23"),
  c("Layer-NVU", "Layer-Exc456"),
  c("Layer-NVU", "Layer-Olg/OPC"),
  c("Layer-Exc23", "Layer-Exc456"),
  c("Layer-Exc23", "Layer-Olg/OPC"),
  c("Layer-Exc456", "Layer-Olg/OPC")
)

# Initialize a vector to store the significance levels
custom_significance_levels <- vector("character", length(comparisons_list))

# Loop through each pair and compute the adjusted p-value
for (i in seq_along(comparisons_list)) {
  ident1 <- comparisons_list[[i]][1]
  ident2 <- comparisons_list[[i]][2]
  
  # Set the identities to the 'layer' metadata column
  Idents(tt.tmp) <- "layer"
  
  # Perform differential expression analysis for 'CD163' between the two groups
  deg_results <- FindMarkers(tt.tmp, 
                             ident.1 = ident1, 
                             ident.2 = ident2, 
                             features = "CD163", 
                             logfc.threshold = 0, 
                             test.use = "wilcox")
  
  # Extract the adjusted p-value
  p_val_adj <- deg_results["CD163", "p_val_adj"]
  
  # Map p-values to significance levels
  if (is.na(p_val_adj)) {
    significance <- "ns"
  } else if (p_val_adj <= 0.001) {
    significance <- "***"
  } else if (p_val_adj <= 0.01) {
    significance <- "**"
  } else if (p_val_adj <= 0.05) {
    significance <- "*"
  } else {
    significance <- "ns"
  }
  
  # Store in the custom_significance_levels vector
  custom_significance_levels[i] <- significance
}

# Now, integrate these significance levels into your plot
p <- p +
  geom_signif(
    comparisons = comparisons_list,
    annotations = custom_significance_levels,
    y_position = c(9, 9.5, 10, 7.5, 8.5, 8),
    tip_length = 0.01,
    textsize = 5,   # Adjust text size as needed
    vjust = 0.5
  )

# Display the plot
print(p)


# Customize the theme to remove background grid and display only axis lines
p <- p + theme(
  panel.background = element_blank(),   # Set the panel background to be blank
  panel.grid.major = element_blank(),   # Remove major grid lines
  panel.grid.minor = element_blank(),   # Remove minor grid lines
  axis.line = element_line(color = "black"),  # Show axis lines
  axis.ticks = element_line(color = "black"),  # Show axis ticks
  plot.background = element_rect(fill = "white", color = NA)  # Set the plot background to white
)

# save as another file only for legend in 10x10 inches
# Add custom legend using geom_text
legend_data <- data.frame(
  x = 1.5,   # Adjust x position as needed
  y = max(gene_data$CD163) * 1.1,  # Place legend above the plot
  label = "* p < 0.05\n** p < 0.01\n*** p < 0.001\n**** p < 0.0001"
)

p <- p + geom_text(data = legend_data, aes(x, y, label = label, color = "Significance"),
                   inherit.aes = FALSE)

# Add custom legend explanation
p + scale_color_manual(name = "Legend", values = c("Significance" = "black"),
                       labels = c("Significance" = "* p < 0.05\n** p < 0.01\n*** p < 0.001\n**** p < 0.0001")) +
  guides(color = guide_legend(title = "Significance Levels"))
