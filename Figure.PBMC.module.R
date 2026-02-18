rm(list = ls(all = TRUE))
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratWrappers)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
})

################################################################################
# Load Seurat object
################################################################################
file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output/Blood"
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"

load(paste0(file_path, "/step1.2.umap.fine.RData"))  # creates 'test'
Idents(test) <- "condition"
load(paste0(file_path, "/step2.1.DEG.singlecell.RData"))

patient_mapping <- c(
  "Epileptic_101923" = "Patient 1",
  "Epileptic_101323" = "Patient 2",
  "Epileptic_22924"  = "Patient 3",
  "Epileptic_2849"   = "Patient 4",
  "Healthy_092524_400" = "Healthy 1",
  "Healthy_092524_607" = "Healthy 2",
  "Healthy_101124"     = "Healthy 3",
  "Healthy_102324"     = "Healthy 4",
  "Healthy_61023"      = "Healthy 5",
  "Healthy_71324"      = "Healthy 6",
  "Healthy_91724"      = "Healthy 7"
)

# Add the new 'patient_id' column to metadata
test@meta.data$patient_id <- patient_mapping[test@meta.data$patient]
test@meta.data$patient_id <- factor(
  test@meta.data$patient_id,
  levels = c(
    "Patient 1", "Patient 2", "Patient 3", "Patient 4",
    "Healthy 1", "Healthy 2", "Healthy 3", "Healthy 4", "Healthy 5", "Healthy 6", "Healthy 7"
  )
)
test@meta.data$patient <- test@meta.data$patient_id

################################################################################
# Gene Sets
################################################################################
gene_sets <- list(
  Myeloid_Inflammation = c(
    "IL1B","TNF","IL6","CCL2","CCL3","CCL4","CCL5","CXCL8","PTGS2",
    "S100A8","S100A9","NLRP3","ICAM1","TLR2","TLR4"
  ),
  Myeloid_Antigen = c(
    "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
    "CD74","CD86","CD80","CIITA","B2M","TAP1","TAP2","PSMB8","PSMB9","CD40"
  ),
  Myeloid_IFN_Antiviral = c(
    "ISG15","MX1","IFIT1","IFIT3","OAS1","OAS2","RSAD2",
    "IRF7","STAT1","STAT2","BST2","DDX58","IFIH1"
  ),
  Myeloid_M1 = c(
    "IL1B","TNF","IL6","IL12A","IL12B",
    "CXCL9","CXCL10","CCL2","CCL3","CCL4",
    "STAT1","IRF5","SOCS3",
    "CCR7","CD80","CD86","CD40",
    "HLA-DRA","HLA-DRB1","B2M","ICAM1","TLR2","TLR4"
  ),
  Myeloid_M2 = c(
    "MRC1","CD163","MSR1","MERTK","STAB1","TGM2",
    "CCL17","CCL22","CCL18",
    "IL10","TGFB1",
    "PPARG","KLF4","LILRB1","LILRB2","C1QA","C1QB","C1QC"
  ),
  
  # B cells
  Bcell_Activation = c(
    "MS4A1","CD79A","CD79B","PAX5","BLK","BANK1","SYK","LYN","CD19","CD22","CD72","SPIB","IRF8"
  ),
  Bcell_Plasma = c(
    "PRDM1","XBP1","IRF4","MZB1","SDC1","IGKC","IGHM","IGHG1","IGHA1","CD38"
  ),
  
  # T cells
  Tcell_TCR_Activation = c(
    "CD3D","CD3E","CD3G","CD247","ZAP70","LCK","LAT","VAV1","IL7R","CCR7","SELL","TCF7","LEF1","CD28","ICOS","PDCD1"
  ),
  Tcell_Cytotoxic = c(
    "GZMB","GZMH","GZMA","PRF1","GNLY","NKG7","IFNG","FASLG","KLRK1","CD160","CCL5","XCL1"
  ),
  Tcell_Treg = c(
    "FOXP3","IL2RA","CTLA4","TIGIT","IKZF2","ENTPD1","TNFRSF18","LAG3","HAVCR2","CCR8"
  ),
  
  # NK
  NK_Cytotoxic = c(
    "GZMB","GZMH","GZMA","PRF1","GNLY","NKG7","KLRD1","KLRC1","KLRC2","NCR1","FCGR3A","TYROBP"
  ),
  
  # Platelet
  Platelet_Coagulation = c(
    "ITGA2B","ITGB3","GP1BA","GP1BB","GP9","PF4","PPBP","SELP","VWF","F13A1","THBS1","PROS1","F2R"
  )
)

################################################################################
# AddModuleScore
################################################################################
added_modules <- c()
for (set_name in names(gene_sets)) {
  test <- AddModuleScore(
    object   = test,
    features = list(gene_sets[[set_name]]),
    name     = set_name,
    nbin     = 24,
    ctrl     = 100,
    assay    = DefaultAssay(test)
  )
  new_cols <- grep(paste0("^", set_name), colnames(test@meta.data), value = TRUE)
  last_col <- tail(new_cols, 1)
  colnames(test@meta.data)[match(last_col, colnames(test@meta.data))] <- set_name
  added_modules <- c(added_modules, set_name)
}
message(sprintf("[DONE] Added modules: %s", paste(added_modules, collapse = ", ")))

################################################################################
# Pseudobulk
################################################################################

module_cols <- c(
  "Myeloid_Inflammation","Myeloid_Antigen","Myeloid_IFN_Antiviral",
  "Myeloid_M1","Myeloid_M2",
  "Bcell_Activation","Bcell_Plasma",
  "Tcell_TCR_Activation","Tcell_Cytotoxic","Tcell_Treg",
  "NK_Cytotoxic",
  "Platelet_Coagulation"
)
stopifnot(all(module_cols %in% colnames(test@meta.data)))

major_map <- c(
  "CD14_Mono"="Myeloid",
  "CD16_Mono"="Myeloid",
  "cDC1"="Myeloid",
  "cDC2"="Myeloid",
  "pDC"="Myeloid",
  "EMP"="Myeloid",
  "Plasmablast"="Bcell",
  "B_naive"="Bcell",
  "B_memory"="Bcell",
  "B_intermediate"="Bcell",
  "CD4_Naive"="T/NK",
  "CD8_Naive"="T/NK",
  "CD4_Memory"="T/NK",
  "CD8_Memory"="T/NK",
  "CD4_CTL"="T/NK",
  "Treg"="T/NK",
  "MAIT"="T/NK",
  "MKI67+_NKT"="T/NK",
  "gdT"="T/NK",
  "NK"="T/NK",
  "NK_CD56bright"="T/NK",
  "Platelets"="Platelet"
)
major_order <- c("Myeloid","Bcell","T/NK","Platelet")

relevant_modules_by_subtype <- list(
  # Myeloid
  "CD14_Mono"      = c("Myeloid_Inflammation","Myeloid_IFN_Antiviral","Myeloid_M1","Myeloid_M2"),
  "CD16_Mono"      = c("Myeloid_Inflammation","Myeloid_IFN_Antiviral","Myeloid_M1","Myeloid_M2"),
  "cDC1"           = c("Myeloid_Antigen","Myeloid_IFN_Antiviral"),
  "cDC2"           = c("Myeloid_Antigen","Myeloid_IFN_Antiviral"),
  "pDC"            = c("Myeloid_IFN_Antiviral"),
  "EMP"            = c("Myeloid_Inflammation","Myeloid_Antigen","Myeloid_M1","Myeloid_M2"),
  
  # B cells
  "Plasmablast"    = c("Bcell_Plasma"),
  "B_naive"        = c("Bcell_Activation"),
  "B_memory"       = c("Bcell_Activation"),
  "B_intermediate" = c("Bcell_Activation"),
  
  # T / NK
  "CD4_Naive"      = c("Tcell_TCR_Activation"),
  "CD8_Naive"      = c("Tcell_TCR_Activation"),
  "CD4_Memory"     = c("Tcell_TCR_Activation"),
  "CD8_Memory"     = c("Tcell_TCR_Activation","Tcell_Cytotoxic"),
  "CD4_CTL"        = c("Tcell_Cytotoxic","Tcell_TCR_Activation"),
  "Treg"           = c("Tcell_Treg"),
  "MAIT"           = c("Tcell_TCR_Activation"),
  "MKI67+_NKT"     = c("Tcell_TCR_Activation","Tcell_Cytotoxic"),
  "gdT"            = c("Tcell_TCR_Activation"),
  "NK"             = c("NK_Cytotoxic"),
  "NK_CD56bright"  = c("NK_Cytotoxic"),
  
  # Platelet
  "Platelets"      = c("Platelet_Coagulation")
)
allow_pairs <- enframe(relevant_modules_by_subtype, name="subtype", value="module") %>%
  unnest_longer(module)

ignore_subtypes <- c("cDC1","EMP","Late_Eryth","Mast","pro_B", "Plasmablast")


meta <- test@meta.data %>%
  as.data.frame() %>%
  mutate(Annotation.fine = as.character(Annotation.fine),
         major = unname(major_map[Annotation.fine])) %>%
  filter(!is.na(major)) %>%
  filter(!Annotation.fine %in% ignore_subtypes)  


pb_long <- meta %>%
  group_by(patient, condition, Annotation.fine, major) %>%
  summarise(across(all_of(module_cols), mean, .names="{.col}"), .groups="drop") %>%
  pivot_longer(cols = all_of(module_cols), names_to = "module", values_to = "score")

pb_long <- pb_long %>%
  inner_join(allow_pairs, by = c("Annotation.fine" = "subtype", "module" = "module"))

pb_long <- pb_long %>%
  mutate(column_key = paste(Annotation.fine, module, sep="|"))

mat <- pb_long %>%
  dplyr::select(patient, column_key, score) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = column_key, values_from = score) %>%
  as.data.frame()

row_info <- meta %>%
  dplyr::distinct(patient, condition) %>%
  dplyr::right_join(dplyr::select(mat, patient), by = "patient") %>%
  dplyr::mutate(condition = factor(condition, levels = c("epileptic","healthy"))) %>%
  dplyr::arrange(condition, patient)

rownames(mat) <- mat$patient
mat <- mat[row_info$patient, , drop=FALSE]
mat$patient <- NULL

col_df <- tibble(column_key = colnames(mat)) %>%
  separate(column_key, into=c("subtype","module"), sep="\\|", remove=FALSE) %>%
  mutate(major = unname(major_map[subtype]))

ord_cols <- c()
for (M in major_order) {
  subtypes_M <- col_df %>% filter(major==M) %>% pull(subtype) %>% unique()
  for (st in subtypes_M) {
    rel_mods <- relevant_modules_by_subtype[[st]]
    if (is.null(rel_mods)) next
    cols_here <- col_df %>%
      filter(major==M, subtype==st, module %in% rel_mods) %>%
      mutate(mod_rank = match(module, rel_mods)) %>%
      arrange(mod_rank) %>% pull(column_key)
    ord_cols <- c(ord_cols, cols_here)
  }
}
ord_cols <- unique(ord_cols)
mat <- mat[, intersect(ord_cols, colnames(mat)), drop=FALSE]
col_df <- col_df[match(colnames(mat), col_df$column_key), ]

pvals <- sapply(colnames(mat), function(cc) {
  v <- mat[, cc]; g <- row_info$condition
  x <- v[g=="epileptic"]; x <- x[!is.na(x)]
  y <- v[g=="healthy"];   y <- y[!is.na(y)]
  if (length(x) >= 3 && length(y) >= 3) {
    tryCatch(wilcox.test(x, y)$p.value, error=function(e) NA_real_)
  } else NA_real_
})
pstars <- ifelse(is.na(pvals), "",
                 ifelse(pvals < 1e-3, "***",
                        ifelse(pvals < 1e-2, "**",
                               ifelse(pvals < 5e-2, "*", ""))))

mat_scaled <- apply(mat, 2, function(z) {
  if (all(is.na(z))) return(z)
  (z - mean(z, na.rm=TRUE)) / sd(z, na.rm=TRUE)
})
mat_scaled <- as.matrix(mat_scaled)
rownames(mat_scaled) <- rownames(mat)

row_split <- row_info$condition
col_split <- factor(col_df$major, levels=major_order)
bottom_anno <- HeatmapAnnotation(
  " " = anno_text(pstars, gp = gpar(fontsize=6), rot = 0, just = "center"),
  annotation_height = unit(6, "mm")
)

top_anno <- HeatmapAnnotation(
  subtype = anno_text(col_df$subtype, gp=gpar(fontsize=6), rot=90, just="right"),
  which = "column"
)


col_fun <- colorRamp2(c(-2,0,2), c("#4575B4","#F7F7F7","#D73027"))
na_col <- "#EEEEEE"

colnames(mat_scaled) <- col_df$module

ht <- Heatmap(
  mat_scaled,
  name = "Z",
  col = col_fun,
  na_col = na_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  column_split = col_split,
  row_gap = unit(1, "mm"),
  top_annotation = top_anno,      
  bottom_annotation = bottom_anno,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize=6),   # 行名字体 6pt
  show_column_names = TRUE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize=6), # 列名字体 6pt
  heatmap_legend_param = list(title="Z-score", at=c(-2,0,2), 
                              title_gp = gpar(fontsize=6), 
                              labels_gp = gpar(fontsize=6))  # 图例文字 6pt
)

draw(ht, padding=unit(c(5,5,5,5), "mm"))
pdf("C:/Users/88312_8ww6o12/Desktop/heatmap.pdf", width=4.2, height=3); draw(ht); dev.off()


