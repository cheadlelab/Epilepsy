rm(list = ls(all = TRUE))
library(ggVennDiagram)
library(ggplot2)
library(dplyr)


##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
file_path <- "//grid/cheadle_home/qianyu/Epilepsy/output"
# file_path <- "/home/qianyu/Desktop/grid/Epilepsy/output"
load(paste0(file_path, "/step2.1.Comparisons.sc2.RData"))

cell_types <- names(final_results)

cell_types <- c(
  "Exc_L23IT", "Exc_L234IT_ERBB4", "Exc_L4IT", "Exc_L4IT_PLCH1", "Exc_L5ET", "Exc_L5IT_GRIN3A",
  "Exc_L5IT_GRIN3A-", "Exc_L56IT_CAR3", "Exc_L56NP", "Exc_L6CT", "Exc_L6IT", "Exc_L6b", 
  "Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", 
  "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier", "Astrocyte", "Microglia", 
  "Olg", "OPC", "Macrophage", "NK/T", "Vascular"
)

# deg <- final_results[[cell_type]][["Focal_vs_Distal"]][["deg"]]
# upregulated <- final_results[[cell_type]][["Focal_vs_Distal"]][["upregulated"]]
# downregulated <- final_results[[cell_type]][["Focal_vs_Distal"]][["downregulated"]]


exc_cell_types <- c("Exc_L23IT", "Exc_L4IT", "Exc_L5IT_GRIN3A", "Exc_L6IT")
inh_cell_types <- c("Inh_SST", "Inh_PVALB", "Inh_VIP", "Inh_CXCL14", "Inh_LAMP5", "Inh_LAMP5_LHX6", "Chandelier")
glia_cell_types <- c("Astrocyte", "Microglia", "Olg", "OPC")


celltype <- c("Exc_L23IT", "Exc_L5IT_GRIN3A", "Exc_L6IT")
direction <- "upregulated" #deg upregulated downregulated
condition <- "Stimulated_vs_Distal" # Focal_vs_Distal Stimulated_vs_Distal


up_list <- list()

# Fill the list with upregulated genes from each inhibitory cell type
for (cell_type in celltype) {
  if (condition %in% names(final_results[[cell_type]])) {
    deg <- final_results[[cell_type]][[condition]][[direction]] %>%
      dplyr::filter(
        abs(avg_log2FC) >= log2(1.5),
        p_val_adj < 0.05
      )
    #up_list[[cell_type]] <- deg$gene
    up_list[[cell_type]] <- deg$gene[!grepl("^MT-", deg$gene)]
  }
}

gene_lists <- up_list


p <- ggVennDiagram(up_list, force_upset = TRUE, relative_height = 2)
p


# pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/figures/manuscript/Figure2B.pdf", width = 5, height = 3)
# print(p)
# dev.off() # Close the PDF device


venn_plot <- ggVennDiagram(
  gene_lists,
  label_alpha = 0,
  edge_size = 0.5,
  label = "count",
  set_size = 4.333,
  category.names = names(gene_lists)
) +
  scale_fill_gradient(low = "#F7FBFF", high = "#E78C8C") +
  theme_void() +
  theme(
    legend.position = "none",
    text = element_text(size = 7),
    plot.title = element_text(size = 7),
    plot.caption = element_text(size = 7)
  ) +
  theme(  # for ggVennDiagram-specific labels
    set.text = element_text(size = 7),
    region.text = element_text(size = 7)
  )

venn_plot
common_genes <- Reduce(intersect, up_list)
common_genes



pdf(file = "//grid/cheadle_home/qianyu/Epilepsy/figures/manuscript/Figure.IT.venn.pdf", width = 2.5, height = 2)
print(venn_plot)
dev.off() # Close the PDF device

