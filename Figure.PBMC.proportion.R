rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(scales)


#################################################################################################################################
############################################################# snRNA #############################################################
#################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/blood"
save_path <- "//grid/cheadle_home/qianyu/Epilepsy/figures/manuscript"

load(paste0(file_path, "/step1.2.umap.fine.RData"))

# Prepare the Seurat object
`%!in%` <- Negate(`%in%`)
test <- subset(test, subset = Annotation.fine %!in% c("B/cDC2", "B/Myeloid", "B/NKT", "B/other", "B/platelets", "Myeloid/NKT", "Myeloid/platelets"))

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

test@active.ident <- as.factor(test@meta.data$patient_id)

ordered_levels <- c(
  ## Myeloid lineage
  "CD14_Mono",      # Classical monocytes
  "CD16_Mono",      # Non-classical monocytes
  "cDC1",           # Conventional DC subset 1
  "cDC2",           # Conventional DC subset 2
  "pDC",            # Plasmacytoid DCs
  "Mast",           # Mast cells
  
  ## Lymphoid lineage
  ### B-cell series (progenitor → mature → effector)
  "pro_B",          
  "B_naive",        
  "B_intermediate", 
  "B_memory",       
  "Plasmablast",    
  
  ### T-cell series
  "CD4_Naive",      # Naïve CD4+
  "CD4_Memory",     # Memory CD4+
  "Treg",           # Regulatory T cells
  "CD4_CTL",        # CD4+ cytotoxic
  "CD8_Naive",      # Naïve CD8+
  "CD8_Memory",     # Memory CD8+
  "MAIT",           # Mucosal-associated invariant T
  "gdT",            # γδ T cells
  
  ### NK-cell series
  "NK_CD56bright",  # CD56bright NK
  "NK",             # CD56dim/adaptive NK
  
  ## Proliferating populations
  "MKI67+_NKT",     # Cycling NKT
  
  ## Erythroid & platelet lineages
  "EMP",            # Erythro-myeloid progenitors
  "Late_Eryth",     # Late erythroid cells
  "Platelets"       # Platelets
)

test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)
print(levels(test@meta.data$Annotation.fine))


default_colors <- hue_pal()(25)
names(default_colors) <- levels(test@meta.data$Annotation.fine)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
Idents(test) <- 'Annotation.fine'
test@meta.data$patient_condition <- test@meta.data$patient_id

pt <- table(Idents(test), test$patient_condition)
pt <- as.data.frame(pt)
colnames(pt) <- c("CellType", "Donor", "Freq")

library(dplyr)

pt <- pt %>%
  group_by(Donor) %>%
  mutate(Percentage = Freq / sum(Freq))

pdf(file = paste0(save_path, "/Figure.PBMC.proportion.pdf"), width = 1.8, height = 1.5)

p1 <- ggplot(pt, aes(x = Donor, y = Percentage, fill = CellType)) +
  geom_col(width = 0.8) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  xlab("Donor") +
  ylab("Proportion") +
  theme_classic(base_size = 6) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 6)
  )

print(p1)
dev.off()
