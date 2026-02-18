################################################################################
# Step 0: Load dependencies
library(CellChat)
library(patchwork)
library(ggplot2)
library(reticulate)
library(ComplexHeatmap)
library(circlize)
library(grid)
source("//grid/cheadle_home/qianyu/Epilepsy/manuscript/V3/netVisual_chord_gene_plus.R")

# Step 1: Load CellChat objects and merge
os_name <- Sys.info()[["sysname"]]
file_path <- switch(
  os_name,
  "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output", 
  "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",
  stop(sprintf("Unsupported OS: %s", os_name))
)

load(paste0(file_path, "/step3.CellChat.Focal.RData"))  # loads Focal, Control, etc.
load(paste0(file_path, "/step3.CellChat.Distal.RData"))  # loads Focal, Control, etc.
load(paste0(file_path, "/step3.CellChat.Stimulated.RData"))  # loads Focal, Control, etc.


cellchat.F <- updateCellChat(Focal)
cellchat.D <- updateCellChat(Distal)
cellchat.S <- updateCellChat(Stimulated)

cellchat.F <- netAnalysis_computeCentrality(cellchat.F)
cellchat.D <- netAnalysis_computeCentrality(cellchat.D)
cellchat.S <- netAnalysis_computeCentrality(cellchat.S)



cellchat.F <- updateCellChat(cellchat.F)
cellchat.D <- updateCellChat(cellchat.D)
cellchat.S <- updateCellChat(cellchat.S)

object.FD <- list(Distal = cellchat.D, Focal = cellchat.F)
object.SD <- list(Distal = cellchat.D, Stimu = cellchat.S)
object.FS <- list(Stimu = cellchat.S, Focal = cellchat.F)
object.All <- list(Focal = cellchat.F, Stimu = cellchat.S, Distal = cellchat.D)

cellchat.FD <- mergeCellChat(object.FD, add.names = names(object.FD))
cellchat.SD <- mergeCellChat(object.SD, add.names = names(object.SD))
cellchat.FS <- mergeCellChat(object.FS, add.names = names(object.FD))
cellchat.All <- mergeCellChat(object.All, add.names = names(object.All))

###############################################################################
# Compare the total number of interactions and interaction strength
###############################################################################
gg1 <- compareInteractions(cellchat.All, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat.All, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2


###############################################################################
# Compare the number of interactions and interaction strength among different cell populations
###############################################################################
par(mfrow = c(1,3), xpd=TRUE)
# netVisual_diffInteraction(cellchat.FD, weight.scale = T)
netVisual_diffInteraction(cellchat.FD, weight.scale = T, measure = "weight")
netVisual_diffInteraction(cellchat.SD, weight.scale = T, measure = "weight")
netVisual_diffInteraction(cellchat.FS, weight.scale = T, measure = "weight")

netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

###############################################################################
# Figure A

gg1 <- netVisual_heatmap(
  cellchat.FD,
  measure = "weight",
  sources.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP",
                  "Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte", "Other Astrocyte", "Microglia","OPC","OxPos OL", "OL"),
  targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP"),
  row.show    = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP",
                  "Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte", "Other Astrocyte", "Microglia","OPC","OxPos OL", "OL"),
  col.show    = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP"),
  remove.isolate = TRUE,
  title.name = "Focal vs. Distal")
gg2 <- netVisual_heatmap(
  cellchat.SD,
  measure = "weight",
  sources.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP",
                  "Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte", "Other Astrocyte", "Microglia","OPC","OxPos OL", "OL"),
  targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP"),
  row.show    = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP",
                  "Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte", "Other Astrocyte", "Microglia","OPC","OxPos OL", "OL"),
  col.show    = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP"),
  remove.isolate = TRUE,
  title.name = "Stimulated vs. Distal")
gg3 <- netVisual_heatmap(
  cellchat.FS,
  measure = "weight",
  sources.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP",
                  "Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte", "Other Astrocyte", "Microglia","OPC","OxPos OL", "OL"),
  targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP"),
  row.show    = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP",
                  "Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte", "Other Astrocyte", "Microglia","OPC","OxPos OL", "OL"),
  col.show    = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP"),
  remove.isolate = TRUE,
  title.name = "Focal vs. Stimulated")
gg1 + gg2 + gg3

# 8.5 4.2

cellchat.FS@meta$Annotation.new <- as.numeric(cellchat.FS@meta$Annotation.new)
netVisual_heatmap1(cellchat.FS, measure = "weight")
# 6.2 6


# Chord diagram
c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT", "Inh_SST","Inh_PVALB","Inh_VIP",
  "Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte", "Other Astrocyte", "Microglia","OPC","OxPos OL", "OL")


#######################################################################################################
# Focal vs Distal, Glia -> Inh
# 15*10
pos.dataset = "Focal"
# pos.dataset = "Stimu"
features.name = pos.dataset
cellchat <- cellchat.FD
object.list <- object.FD

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Focal",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Distal",ligand.logFC = -0.1, receptor.logFC = -0.1)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[["Focal"]],
                     sources.use = c("Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte","OPC"),
                     targets.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.up,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))
netVisual_chord_gene(object.list[["Distal"]],
                     sources.use = c("Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte","OPC"),
                     targets.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.down,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#######################################################################################################
# Focal vs Distal, Inh -> Exc
# 15*10
pos.dataset = "Focal"
# pos.dataset = "Stimu"
features.name = pos.dataset
cellchat <- cellchat.FD
object.list <- object.FD

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Focal",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Distal",ligand.logFC = -0.1, receptor.logFC = -0.1)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat.FD)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat.FD)


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[["Focal"]],
                          targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT"),
                          sources.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                          slot.name = 'net',
                          net = net.up,
                          lab.cex = 0.8,
                          small.gap = 3.5,
                          title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))
netVisual_chord_gene(object.list[["Distal"]],
                          targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT"),
                          sources.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                          slot.name = 'net',
                          net = net.down,
                          lab.cex = 0.8,
                          small.gap = 3.5,
                          title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


#######################################################################################################
# Stimu vs Distal, Inh -> Exc
# 15*10
pos.dataset = "Stimu"
# pos.dataset = "Stimu"
features.name = pos.dataset
cellchat <- cellchat.SD
object.list <- object.SD

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Stimu",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Distal",ligand.logFC = -0.1, receptor.logFC = -0.1)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat.FD)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat.FD)


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[["Stimu"]],
                     targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT"),
                     sources.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.up,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))
netVisual_chord_gene(object.list[["Distal"]],
                     targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT"),
                     sources.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.down,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))




#######################################################################################################
# Stimu vs Distal, Inh -> Exc
# 15*10
pos.dataset = "Focal"
# pos.dataset = "Stimu"
features.name = pos.dataset
cellchat <- cellchat.FS
object.list <- object.FS

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Focal",ligand.logFC = 0.2, receptor.logFC = NULL)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat.FD)


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[["Focal"]],
                     targets.use = c("Exc_L23IT","Exc_L4IT","Exc_L5IT_GRIN3A","Exc_L6IT"),
                     sources.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.up,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))


#######################################################################################################
# Focal vs Distal, Glia -> Inh
# 15*10
pos.dataset = "Stimu"
# pos.dataset = "Stimu"
features.name = pos.dataset
cellchat <- cellchat.SD
object.list <- object.SD

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Stimu",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Distal",ligand.logFC = -0.1, receptor.logFC = -0.1)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[["Stimu"]],
                     sources.use = c("Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte","OPC"),
                     targets.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.up,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))
netVisual_chord_gene(object.list[["Distal"]],
                     sources.use = c("Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte","OPC"),
                     targets.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.down,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#######################################################################################################
# Focal vs Distal, Glia -> Inh
# 15*10
pos.dataset = "Focal"
# pos.dataset = "Stimu"
features.name = pos.dataset
cellchat <- cellchat.FS
object.list <- object.FS

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Focal",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "Stimu",ligand.logFC = -0.1, receptor.logFC = -0.1)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[["Focal"]],
                     sources.use = c("Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte","OPC"),
                     targets.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.up,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(object.list)[2]))
netVisual_chord_gene(object.list[["Stimu"]],
                     sources.use = c("Lipid-Accumulated Reactive Astrocyte", "Reactive Astrocyte","OPC"),
                     targets.use = c("Inh_SST","Inh_PVALB","Inh_VIP"),
                     slot.name = 'net',
                     net = net.down,
                     lab.cex = 0.8,
                     small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))







