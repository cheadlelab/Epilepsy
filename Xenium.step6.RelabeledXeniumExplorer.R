# https://divingintogeneticsandgenomics.com/post/neighborhood-cellular-niches-analysis-with-spatial-transcriptome-data-in-seurat-and-bioconductor/
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager
# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ
rm(list = ls(all = TRUE))

#load packages
library(Seurat)
library(dplyr)
library(jsonlite)
library(dplyr)
library(arrow)
library(spacexr)
library(ggplot2)
library(dplyr)
library(scales)

save.path <- "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/"

load(paste0(save.path, "step1.list.RData"))

xenium.obj <- sample_list$sample8

xenium.obj$plot_group <- "Other"
xenium.obj$plot_group[xenium.obj$layer == "Layer-Exc23"] <- "Layer-Exc23"
xenium.obj$plot_group[
  xenium.obj$layer == "Layer-Exc23" &
    xenium.obj$predicted.celltype == "Microglia"
] <- "Layer-Exc23_Microglia"



# --- build Xenium CSV (cell_id, group, color) ---
df_xenium <- data.frame(
  cell_id = colnames(xenium.obj),
  group   = as.character(xenium.obj$plot_group),
  stringsAsFactors = FALSE
)

# drop cells with NA group (optional but usually helpful)
df_xenium <- df_xenium[!is.na(df_xenium$group) & df_xenium$group != "", ]

# map group -> hex color
default_colors <- c(
  "Other"                 = "#D3D3D3",
  "Layer-Exc23"           = "#FFD700",  # yellow
  "Layer-Exc23_Microglia" = "#FF0000"   # red
)

df_xenium$color <- unname(default_colors[df_xenium$group])

# final column order exactly as Xenium expects
df_xenium <- df_xenium[, c("cell_id", "group", "color")]

write.csv(df_xenium, file = "//grid/cheadle_home/qianyu/Epilepsy/manuscript/V5/xenium_sample8.csv", row.names = FALSE, quote = FALSE)
