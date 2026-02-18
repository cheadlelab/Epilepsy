# https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq
rm(list = ls(all = TRUE))

library(Seurat)
library(dplyr)
library(Matrix)

# raw_counts <- read.csv(file="/media/qianyu/MOSS1/Epilepsy/AllenData/matrix.csv", header = TRUE)
raw_counts <- read.csv(file="/grid/cheadle/home/qianyu/Epilepsy/OtherDataset/AllenData/matrix.csv", header = TRUE)
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[, -1]
tt <- t(raw_counts)
test <- CreateSeuratObject(count = tt, min.cells = 3, min.features = 200, project = "Allen_Cortical")
test@meta.data$cellid <- colnames(test)

# metadata <- read.csv(file="/media/qianyu/MOSS1/Epilepsy/AllenData/metadata.csv", header = TRUE)
metadata <- read.csv(file="/grid/cheadle/home/qianyu/Epilepsy/OtherDataset/AllenData/metadata.csv", header = TRUE)
metadata[] <- lapply(metadata, function(x) ifelse(is.na(x) | x == "", "Unknown", x))

test@meta.data$class_label <- metadata$class_label[match(test@meta.data$cellid, metadata$sample_name)]
test@meta.data$subclass_label <- metadata$subclass_label[match(test@meta.data$cellid, metadata$sample_name)]
test@meta.data$donor_sex_label <- metadata$donor_sex_label[match(test@meta.data$cellid, metadata$sample_name)]
test@meta.data$region_label <- metadata$region_label[match(test@meta.data$cellid, metadata$sample_name)]
test@meta.data$cluster_label <- metadata$cluster_label[match(test@meta.data$cellid, metadata$sample_name)]
test@meta.data$donor <- metadata$external_donor_name_label[match(test@meta.data$cellid, metadata$sample_name)]


test@meta.data$cluster_label_simplify <- sapply(strsplit(test@meta.data$cluster_label, "\\s+"), function(x) paste(x[1:min(length(x), 2)], collapse=" "))
test@meta.data$sample_name <- test@meta.data$cellid


test <- subset(test, subset = region_label %in% c("MTG", "A1C"))
save_path <- "/grid/cheadle/home/qianyu/Epilepsy/output"
save(test, file = paste0(save_path, "/step0.AllenCortical.RData"))
