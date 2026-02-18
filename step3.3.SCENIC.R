#####################################################################################################
########################################## 0. Setup #################################################
#####################################################################################################
rm(list = ls(all = TRUE))

# Load required packages
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(SCopeLoomR)
library(dplyr)
library(Matrix)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(ggplot2)

#####################################################################################################
########################################## 1. Load Data #############################################
#####################################################################################################
# Detect system and define file path
os_name <- Sys.info()[["sysname"]]
file_path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output",
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",
                    stop(sprintf("Unsupported OS: %s", os_name))
)

# Load Seurat object
load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))
test <- JoinLayers(test)
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
DefaultAssay(test) <- "RNA"

#####################################################################################################
##################################### 2. Prepare Input for SCENIC ###################################
#####################################################################################################
# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scenic.html

exprMat <- as.matrix(GetAssayData(test, assay = "RNA", slot = "data"))
cellInfo <- test@meta.data

loci1 <- which(rowSums(exprMat) > 1*.01*ncol(exprMat))
exprMat_filter <- exprMat[loci1, ]


add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

# loom <- build_loom(paste0(file_path, "/step3.3.SCENIC.exprMat_filtered.loom"), dgem=exprMat_filter)
loom <- build_loom("/home/qianyu/Desktop/step3.3.SCENIC.exprMat_filtered.loom", dgem=exprMat_filter)

loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

#################################################################################################
##################################### 3. in conda environment ###################################
#################################################################################################

# https://github.com/aertslab/pySCENIC/blob/master/resources/hs_hgnc_tfs.txt

# 1) GRN
pyscenic grn \
/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.exprMat_filtered.loom \
/home/qianyu/Desktop/grid/resource/scenic/hs_hgnc_tfs.txt \
-o /home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.adj.csv \
--num_workers 32 --method grnboost2 --seed 777

# 2) cisTarget
pyscenic ctx \
/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.adj.csv \
/home/qianyu/Desktop/grid/resource/scenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
/home/qianyu/Desktop/grid/resource/scenic/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /home/qianyu/Desktop/grid/resource/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.exprMat_filtered.loom \
--output /home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.regulons.csv \
--mask_dropouts \
--num_workers 32

# 3) AUCell
pyscenic aucell \
/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.exprMat_filtered.loom \
/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.regulons.csv \
--output /home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.output.loom \
--num_workers 32

#################################################################################################
# grid connection slow
# 1) GRN
pyscenic grn \
/home/qianyu/Desktop/step3.3.SCENIC.exprMat_filtered.loom \
/home/qianyu/Desktop/grid/resource/scenic/hs_hgnc_tfs.txt \
-o /home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.adj.csv \
--num_workers 16 --method grnboost2 --seed 777

# 2) cisTarget
pyscenic ctx \
/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.adj.csv \
/home/qianyu/Desktop/grid/resource/scenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
/home/qianyu/Desktop/grid/resource/scenic/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname /home/qianyu/Desktop/grid/resource/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /home/qianyu/Desktop/step3.3.SCENIC.exprMat_filtered.loom \
--output /home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.regulons.csv \
--mask_dropouts \
--num_workers 4

# 3) AUCell
pyscenic aucell \
/home/qianyu/Desktop/step3.3.SCENIC.exprMat_filtered.loom \
/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.regulons.csv \
--output /home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.output.loom \
--num_workers 8

#########################################################################################################
##################################### 4. Retrieve AUC scores per cell ###################################
#########################################################################################################

loom_file <- "/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.output.loom"
loom <- open_loom(loom_file)

regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

saveRDS(
  list(
    'regulons' = regulons,
    'regulonAUC' = regulonAUC
  ),
  file = "/home/qianyu/Desktop/grid/Epilepsy/output/step3.3.SCENIC.output.AUC.rds"
)


#########################################################################################################
##################################### 5. Integrate with Seurat #########################################
#########################################################################################################
rm(list = ls(all = TRUE))

# step 0: Load required libraries
library(Seurat)
library(dplyr)
library(AUCell)
library(reshape2)
library(ggplot2)

# step 1: Load Seurat object and SCENIC AUC matrix
os_name <- Sys.info()[["sysname"]]
file_path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output",
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",
                    stop(sprintf("Unsupported OS: %s", os_name)))

# Load Seurat object containing microglia clusters
load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))  # This loads the object named "test"
Idents(test) <- "seurat_clusters"  # Set cluster identity

# Load SCENIC AUC results (regulon activity)
scenic_data <- readRDS(paste0(file_path, "/step3.3.SCENIC.output.AUC.rds"))

regulonAUC <- scenic_data$regulonAUC
auc_mat <- as.matrix(getAUC(regulonAUC))  # Regulon × Cell matrix

# Clean regulon names (remove brackets/parentheses)
rownames(auc_mat) <- gsub("[()+]", "", rownames(auc_mat))

# Add AUC matrix to Seurat as a new assay
test[["AUC"]] <- CreateAssayObject(data = auc_mat)
DefaultAssay(test) <- "AUC"
test <- ScaleData(test, assay = 'AUC', features = rownames(auc_mat))








tt.tmp <- subset(test, subset = Annotation.fine == "Microglia")

# step 2: Prepare long-format data frame with cluster info
auc_data <- GetAssayData(tt.tmp, assay = "AUC", slot = "data")  # Regulon × Cell
auc_df <- as.data.frame(t(as.matrix(auc_data)))               # Cell × Regulon

Idents(tt.tmp) <- "condition"
auc_df$cluster <- Idents(tt.tmp)                                # Add cluster info

# Make a long-format table: cell × regulon with AUC and cluster
auc_df$cell <- rownames(auc_df)
auc_long <- reshape2::melt(
  auc_df,
  id.vars = c("cell", "cluster"),
  variable.name = "regulon",
  value.name = "AUC"
)

# (Optional) set factor order for clusters
auc_long$cluster <- factor(
  auc_long$cluster,
  levels = c("Focal","Distal","Stimulated")
)
clusters <- levels(auc_long$cluster)

# Define comparison groups: Cluster 6 vs all others
auc_long$group <- ifelse(auc_long$cluster == "Distal", "Distal", "Others")

# step 3: Statistical test — Wilcoxon test for each regulon
auc_stat <- auc_long %>%
  group_by(regulon) %>%
  summarize(
    p_value = wilcox.test(AUC ~ group)$p.value,
    COI_mean = mean(AUC[group == "Distal"]),
    Others_mean = mean(AUC[group == "Others"]),
    .groups = "drop"
  ) %>%
  arrange(p_value)

# step 4: Identify significantly downregulated regulons in Cluster 6
downregulated_regulons <- auc_stat %>%
  filter(p_value < 0.05 & COI_mean < Others_mean)
upregulated_regulons <- auc_stat %>%
  filter(p_value < 0.05 & COI_mean > Others_mean)

# View top results
head(downregulated_regulons, 10)
head(upregulated_regulons, 10)

# step 5: Visualize top downregulated regulons using violin plots
top_down <- head(downregulated_regulons$regulon, 15)  # Top 6 regulons
top_up <- head(upregulated_regulons$regulon, 15)  # Top 6 regulons



# UMAP projection of downregulated regulon activity
DefaultAssay(test) <- "AUC"
FeaturePlot(tt.tmp, features = as.character(top_down), reduction = "umap_harmony")
VlnPlot(tt.tmp, features = as.character(top_down), pt.size = 0)


# Heatmap for top downregulated regulons
DoHeatmap(
  tt.tmp, 
  features = c(top_up, top_down),
  label = TRUE,
  slot = 'scale.data', raster = F) + scale_fill_gradient2(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
    midpoint = 0,
    guide = "colourbar",
    aesthetics = "fill"
  )

