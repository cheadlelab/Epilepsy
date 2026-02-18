rm(list = ls(all = TRUE))

#load packages
library(Seurat)
library(dplyr)
library(jsonlite)
library(dplyr)
library(arrow)
library(spacexr)
library(ggplot2)
library(scales)
options(future.globals.maxSize = 128 * 1024^3)

os_name <- Sys.info()[["sysname"]]   # "Linux", "Windows", "Darwin" …
save.path <- switch(os_name,
                    "Linux"   = "/home/qianyu/Desktop/grid/qianyu/Epilepsy/Xenium/output/",
                    "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/Xenium/output/",
                    stop(sprintf("Unsupported OS: %s", os_name))
)
load(paste0(save.path, "step1.list.RData"))


orig.ident <- c("1C SA", "4E", "1C SB", "5C", "2E", "3E SB", "1E", "4C", "5E", "2C", "3E SA")
patient <- c(1, 4 ,1, 5, 2, 3, 1, 4, 5, 2, 3)
condition <- c("Distal", "Focal", "Distal", "Distal", "Focal", "Focal", "Focal", "Distal", "Focal", "Distal", "Focal")




################################################################################
xenium.obj <- sample_list[["sample2"]]

xenium.obj$plot_group <- "Other"
xenium.obj$plot_group[xenium.obj$layer == "Layer-Exc23"] <- "Layer-Exc23"
xenium.obj$plot_group[
  xenium.obj$layer == "Layer-Exc23" &
    xenium.obj$predicted.celltype == "Microglia"
] <- "Layer-Exc23_Microglia"


p1 <- ImageDimPlot(
  xenium.obj,
  group.by = "plot_group",
  size = 0.3,
  dark.background = FALSE
) +
  scale_fill_manual(values = c(
    "Other"                 = "grey80",
    "Layer-Exc23"           = "#FFD700",   # yellow
    "Layer-Exc23_Microglia" = "red"
  ))

################################################################################
xenium.obj <- sample_list[["sample8"]]

xenium.obj$plot_group <- "Other"
xenium.obj$plot_group[xenium.obj$layer == "Layer-Exc23"] <- "Layer-Exc23"
xenium.obj$plot_group[
  xenium.obj$layer == "Layer-Exc23" &
    xenium.obj$predicted.celltype == "Microglia"
] <- "Layer-Exc23_Microglia"


p2 <- ImageDimPlot(
  xenium.obj,
  group.by = "plot_group",
  size = 0.19,
  dark.background = FALSE
) +
  scale_fill_manual(values = c(
    "Other"                 = "grey80",
    "Layer-Exc23"           = "#FFD700",   # yellow
    "Layer-Exc23_Microglia" = "red"
  ))

p1+p2

# 12*6