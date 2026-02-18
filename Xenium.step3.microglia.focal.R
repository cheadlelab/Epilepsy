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

# in CD163+ microglia, CD163+ expression realated to the layer distribution
tt.tmp <- subset(xenium.microglia, subset = clusters_normal == 2)

#######################################################################################################
tt.tmp <- subset(tt.tmp, subset = condition == "Focal")

#####################################################################################################################################
############################################################## Figures ##############################################################
#####################################################################################################################################
# 5 inch x 5 inch
library(ggplot2)
library(dplyr)

table(tt.tmp@meta.data$layer)
# Create a data frame for cluster 2 (row 2)
cluster2_data <- data.frame(
  Layer = factor(c("Layer-NVU", "Layer-Exc23", "Layer-Exc456", "Layer-Olg/OPC"), 
                 levels = c("Layer-NVU", "Layer-Exc23", "Layer-Exc456", "Layer-Olg/OPC")),
  Count = c(383, 131, 139, 225)
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
p <- p + geom_signif(
  comparisons = list(c("Layer-NVU", "Layer-Exc23"), c("Layer-NVU", "Layer-Exc456"), 
                     c("Layer-NVU", "Layer-Olg/OPC"), c("Layer-Exc23", "Layer-Exc456"), 
                     c("Layer-Exc23", "Layer-Olg/OPC"), c("Layer-Exc456", "Layer-Olg/OPC")),
  map_signif_level = TRUE,
  test = "wilcox.test",
  y_position = c(9, 9.5, 10, 7.5, 8.5, 8)  # Adjusted heights for clarity
)

# Customize the theme to remove background grid and display only axis lines
p <- p + theme(
  panel.background = element_blank(),   # Set the panel background to be blank
  panel.grid.major = element_blank(),   # Remove major grid lines
  panel.grid.minor = element_blank(),   # Remove minor grid lines
  axis.line = element_line(color = "black"),  # Show axis lines
  axis.ticks = element_line(color = "black"),  # Show axis ticks
  plot.background = element_rect(fill = "white", color = NA)  # Set the plot background to white
)
