rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(scales)


#################################################################################################################################
############################################################# snRNA #############################################################
#################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output/Blood"
save_path <- "//grid/cheadle/home/qianyu/Epilepsy/figures"


load(paste0(file_path, "/step1.2.umap.fine.RData"))

library(dplyr)
library(tidyr)

meta <- test@meta.data %>%
  select(patient, condition, Annotation.fine) 

counts <- meta %>%
  group_by(patient, condition, Annotation.fine) %>%
  summarise(n = n(), .groups = "drop")

totals <- meta %>%
  group_by(patient, condition) %>%
  summarise(total = n(), .groups = "drop")

prop_df <- counts %>%
  left_join(totals, by = c("patient", "condition")) %>%
  mutate(prop = n / total)

summary_table <- prop_df %>%
  group_by(Annotation.fine) %>%
  summarise(
    mean_epileptic = mean(prop[condition == "epileptic"]),
    mean_healthy  = mean(prop[condition == "healthy"]),
    p.value       = wilcox.test(
      prop[condition == "epileptic"],
      prop[condition == "healthy"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH")
  ) %>%
  arrange(p.adj)

print(summary_table)









