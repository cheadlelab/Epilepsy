# https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ
library(spacexr)
library(ggplot2)

xenium.obj <- sample_1


query.counts <- GetAssayData(xenium.obj, assay = "RNA", slot = "counts")[, Cells(xenium.obj[["RNA"]])]
coords <- GetTissueCoordinates(xenium.obj[["fov"]], which = "centroids")

rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# allen.corted.ref can be downloaded here:
# https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
load("/run/user/1000/gvfs/smb-share:server=grid-hs,share=cheadle_home,user=qianyu/qianyu/Epilepsy/output/step2.Merge.relabel.RData")
cheadle.epilepsy.ref <- test
test <- NULL
cheadle.epilepsy.ref <- UpdateSeuratObject(cheadle.epilepsy.ref)

Idents(cheadle.epilepsy.ref) <- "Annotation"
# remove CR cells because there aren't enough of them for annotation
# cheadle.epilepsy.ref <- subset(cheadle.epilepsy.ref, subset = Annotation != "CR")
counts <- GetAssayData(cheadle.epilepsy.ref, assay = "RNA", slot = "counts")
cluster <- as.factor(cheadle.epilepsy.ref$Annotation)
names(cluster) <- colnames(cheadle.epilepsy.ref)
nUMI <- cheadle.epilepsy.ref$nCount_RNA
names(nUMI) <- colnames(cheadle.epilepsy.ref)
nUMI <- colSums(counts)
levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)
# run RCTD with many cores
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
xenium.obj <- subset(xenium.obj, cells = keep.cells)

xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov", group.by = "predicted.celltype",
                              niches.k = 5, neighbors.k = 30)

celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype", size = 1.5, cols = "polychrome",
                              dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot