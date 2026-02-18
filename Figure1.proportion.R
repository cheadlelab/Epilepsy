rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(scales)


#################################################################################################################################
############################################################# snRNA #############################################################
#################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
save_path <- "C:\\Users\\thech\\Desktop"


load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))

# Prepare the Seurat object
test <- JoinLayers(test)
test[["RNA"]] <- as(object = test[["RNA"]], Class = "Assay")
test <- subset(test, subset = Annotation.fine != "Other")
DefaultAssay(test) <- "RNA"

test@active.ident <- as.factor(test@meta.data$orig.ident)

GABAergic <- c("Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier")
Glu <- c("Exc_L234IT_ERBB4",  # L2-3-4
         "Exc_L23IT",         # L2/3
         "Exc_L4IT",          # L4
         "Exc_L4IT_PLCH1",    # L4
         "Exc_L5ET",          # L5
         "Exc_L5IT_GRIN3A",    # L5
         "Exc_L5IT_GRIN3A-",   # L5
         "Exc_L56NP",         # L5-6
         "Exc_L56IT_CAR3",    # L5-6
         "Exc_L6IT",          # L6
         "Exc_L6CT",          # L6
         "Exc_L6b")           # L6
glia <- c("Astrocyte", "OPC", "Olg", "Microglia", "Vascular")
other <- c( "Macrophage", "NK/T")

ordered_levels <- c(GABAergic, Glu, glia, other)
test@meta.data$Annotation.fine <- factor(test@meta.data$Annotation.fine, levels = ordered_levels)
print(levels(test@meta.data$Annotation.fine))


default_colors <- hue_pal()(26)
names(default_colors) <- levels(test@meta.data$Annotation.fine)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
Idents(test) <- 'Annotation.fine'
test@meta.data$patient_condition <- paste(
  test@meta.data$patient,
  substr(test@meta.data$condition, 1, 1),
  sep = ""
)


pt <- table(Idents(test), test$patient_condition)
pt <- as.data.frame(pt)
colnames(pt) <- c("CellType", "Donor", "Freq")

library(dplyr)

pt <- pt %>%
  group_by(Donor) %>%
  mutate(Percentage = Freq / sum(Freq))

pdf(file = paste0(save_path, "/Figure1G.pdf"), width = 2, height = 1.8)

p1 <- ggplot(pt, aes(x = Donor, y = Percentage, fill = CellType)) +
  geom_col(width = 0.8) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  xlab("Donor") +
  ylab("Proportion") +
  theme_classic(base_size = 7) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 7)
  )

print(p1)
dev.off()
