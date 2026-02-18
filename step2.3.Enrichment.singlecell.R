rm(list = ls(all = TRUE))

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)

file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"

load(paste0(file_path, "/step2.1.DEG.singlecell.RData"))

ontologies <- c("BP", "CC", "MF")

GO.up <- list()
GO.down <- list()

KEGG.up <- list()
KEGG.down <- list()

for (celltype in names(patient_regressed)) {
  tmp <- patient_regressed[[celltype]]
  
  pos.tmp <- tmp %>% filter(p_val_adj < 0.05, avg_log2FC > log2(1.2))
  neg.tmp <- tmp %>% filter(p_val_adj < 0.05, avg_log2FC < -log2(1.2))
  geneListUp <- rownames(pos.tmp)
  geneListDown <- rownames(neg.tmp)
  
  for (ont in ontologies) {
    up.ego <- enrichGO(
      gene          = geneListUp,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    
    down.ego <- enrichGO(
      gene          = geneListDown,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    
    GO.up[[celltype]][[ont]] <- up.ego
    GO.down[[celltype]][[ont]] <- down.ego
  }
  
  gene.df.up <- bitr(geneListUp, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  gene.df.down <- bitr(geneListDown, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  up.kegg <- enrichKEGG(
    gene         = gene.df.up$ENTREZID,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  down.kegg <- enrichKEGG(
    gene         = gene.df.down$ENTREZID,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  KEGG.up[[celltype]] <- up.kegg
  KEGG.down[[celltype]] <- down.kegg
}

save(GO.up, GO.down, KEGG.up, KEGG.down, file = paste0(file_path, "/step2.3.Enrichment.singlecell.RData"))





barplot(KEGG.down[["regular_T"]], showCategory = 20)

a <- barplot(GO.up[["Macrophage"]][["BP"]])
b <- barplot(GO.up[["Macrophage"]][["CC"]])
c <- barplot(GO.up[["Macrophage"]][["MF"]])
a+b+c

a <- barplot(KEGG.up[["Macrophage"]])
b <- barplot(KEGG.up[["Macrophage"]])
c <- barplot(KEGG.up[["Macrophage"]])
a+b+c
a <- barplot(GO.down[["B"]][["BP"]])
b <- barplot(GO.down[["B"]][["CC"]])
c <- barplot(GO.down[["B"]][["MF"]])
a+b+c