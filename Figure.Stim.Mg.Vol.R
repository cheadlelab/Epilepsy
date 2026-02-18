rm(list = ls(all = TRUE))
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(ggplot2)
library(stringr)

##################################################################################################################################
# Load and Preprocess Data
##################################################################################################################################
os_name <- Sys.info()[["sysname"]]   # "Linux", "Windows", "Darwin" …
file_path <- switch(
  os_name,
  "Linux"   = "/home/qianyu/Desktop/grid/Epilepsy/output", 
  "Windows" = "//grid/cheadle_home/qianyu/Epilepsy/output",
  stop(sprintf("Unsupported OS: %s", os_name))
)

# load(paste0(file_path, "/step1.3.annotation.AC.seurat4.RData"))
load(paste0(file_path, "/step1.2.relabel.cheadle.RData"))
test <- JoinLayers(test)

tmp <- subset(test, subset = Annotation.fine %in% c("Microglia"))


DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- 'condition'
markers <- FindAllMarkers(object = tmp, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0)

markers <- FindMarkers(tmp, 
                       ident.1 = "Stimulated",
                       ident.2 = "Distal",
                       min.pct = 0.1,
                       test.use = "MAST",
                       latent.vars = "patient")

markers$symbol <- rownames(markers)
markers$gene <- rownames(markers)


########################################################################################################################################################
##################################################################### Volcano plot #####################################################################
########################################################################################################################################################

fc_thr   <- log2(1.5)
eps      <- 1e-300
padj_thr <- 0.05

pattern_rm <- "^(AC|AL|AP)[0-9]+(?:\\.[0-9]+)?$|^LINC|.-AS[0-9]+$"

df <- markers %>%
  mutate(log10fdr = -log10(p_val_adj + eps)) %>%
  filter(!grepl(pattern_rm, gene, ignore.case = TRUE)) %>%
  distinct(gene, .keep_all = TRUE) %>% 
  filter(!grepl("^MT-", gene))


force_gene <- c("")

force_row <- df %>%
  dplyr::filter(toupper(gene) %in% toupper(force_gene))

df2 <- df %>%
  mutate(sig = (p_val_adj < padj_thr) & (abs(avg_log2FC) >= fc_thr))

top_up <- df2 %>%
  filter(sig, avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10)

top_down <- df2 %>%
  filter(sig, avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 5)
top_lab <- bind_rows(top_up, top_down, force_row) %>% distinct(gene, .keep_all = TRUE)


mt_in_top <- top_lab$gene[grepl("^MT-", top_lab$gene)]


df$category <- "Other"
df$category[df$gene %in% c(top_up$gene, force_gene)] <- "Top10_Up"
df$category[df$gene %in% top_down$gene] <- "Top10_Down"
# df$category[df$gene %in% mt_in_top] <- "MT_Top"
# cols <- c(Other = "grey75", Top10_Up = "#d62728", Top10_Down = "#1f77b4", MT_Top = "#ff7f0e")
# label_data <- top_lab[!grepl("^MT-", top_lab$gene), , drop = FALSE]
label_data <- top_lab
cols <- c(Other = "grey75", Top10_Up = "#d62728", Top10_Down = "#1f77b4")


p <- ggplot(df, aes(x = avg_log2FC, y = log10fdr)) +
  geom_point(aes(color = category), size = 0.8, alpha = 0.85) +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed",
             size = 0.3, color = "grey50") +
  labs(title = "Microglia (Stim vs non-Focal)",
       x = "avg_log2FC",
       y = expression(-log[10]("adj. p"))) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold"))

pt_to_gg <- function(pt) pt / 2.845

p1 <- p + ggrepel::geom_text_repel(
  data = label_data,
  aes(label = gene),
  size = pt_to_gg(7),
  seed = 233,
  box.padding = 0.2,
  point.padding = 0.1,
  min.segment.length = 0,
  segment.size = 0.2,
  max.overlaps = Inf
)





pdf(file = "//grid/cheadle_home/shared/Qianyu_shared/Epilepsy/manuscript/V7/Figure.Stim.Mg.pdf", width = 2, height = 2)
print(p1)
dev.off() # Close the PDF device
