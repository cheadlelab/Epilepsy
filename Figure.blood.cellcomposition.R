rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(scales)


#################################################################################################################################
############################################################# scRNA #############################################################
#################################################################################################################################
file_path <- "/grid/cheadle/home/qianyu/Epilepsy/output/Blood"
save_path <- "/grid/cheadle/home/qianyu/Epilepsy/figures"


load(paste0(file_path, "/step1.2.umap.fine.RData"))

# Prepare the Seurat object
`%!in%` <- Negate(`%in%`)
test <- subset(test, subset = Annotation.fine %!in% c("Macrophage/platelets", "Myeloid/NKT", "RBC", "Platelets"))


test@active.ident <- as.factor(test@meta.data$orig.ident)
ordered_levels <- c(
  # Myeloid lineage (ordered by differentiation hierarchy)
  "FCGR3A+ Mono",        # Classical monocytes (CD16+)
  "Macrophage",          # Tissue-resident macrophages
  # "Macrophage/platelets",# Macrophage-platelet interacting population
  "DC",                  # Conventional dendritic cells
  "cDCs",                # Classical dendritic cells
  "pDC",                 # Plasmacytoid dendritic cells
  "Mast",                # Mast cells
  # "Myeloid/NKT",         # Myeloid-NKT hybrid population
  
  # Lymphoid lineage
  ## T-cell series (ordered by differentiation state)
  "Naive_CD4+_T",        # Naïve CD4+ T cells
  "Memory_CD4+_T",       # Memory CD4+ T cells
  "regular_T",           # Conventional T cells (unspecified subset)
  "CD8+_T_1",            # CD8+ T cell subset 1 (likely effector)
  "CD8+_T_2",            # CD8+ T cell subset 2 (memory)
  "CD8+_T_3",            # CD8+ T cell subset 3 (exhausted)
  "gamma_delta_T",       # γδ T cells
  
  ## NK-cell series 
  "NK_1",                # NK subset 1 (CD56bright)
  "NK_2",                # NK subset 2 (CD56dim)
  "NK_3",                # NK subset 3 (adaptive)
  
  ## B-cell series (developmental order)
  "Pre_B",               # Pre-B cells (immature)
  "B",                   # Mature B cells
  "Plasma_B",            # Plasma cells (antibody-secreting)
  
  # Proliferating populations
  "MKI67+_NKT"           # Proliferating NKT cells (cell cycle active)
  
  # Terminal functional elements
  # "Platelets"          # Platelets (thrombocytes)
  # "RBC"                # Red blood cells (erythrocytes)
)
test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)
print(levels(test@meta.data$Annotation.fine))

default_colors <- hue_pal()(20)
names(default_colors) <- levels(test@meta.data$Annotation.fine)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
patient_mapping <- c(
  "Epileptic_101923" = "P1",
  "Epileptic_101323" = "P2",
  "Epileptic_22924" = "P3",
  "Epileptic_2849" = "P4",
  
  "Healthy_91724"      = "H1",
  "Healthy_101124"     = "H2",
  "Healthy_092524_607" = "H3",
  "Healthy_092524_400" = "H4",
  "Healthy_102324"     = "H5",
  "Healthy_61023"      = "H6",
  "Healthy_71324"      = "H7"
)

# Add the new 'patient_id' column to metadata
test@meta.data$patient_id <- ifelse(
  test@meta.data$patient %in% names(patient_mapping),
  patient_mapping[test@meta.data$patient],
  as.character(test@meta.data$patient)  # Keep original ID if not in mapping
)
unique(test@meta.data$patient_id)



Idents(test) <- 'Annotation.fine'

pt <- table(Idents(test), test$patient_id)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pdf(file = paste0(save_path, "/Figure.blood.cellcomposition.pdf"), width = 4, height = 4)
p1 <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic() +
  #theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.8) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.position="none") +
  #scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank()) # + coord_flip()
p1
dev.off() # Close the PDF device
