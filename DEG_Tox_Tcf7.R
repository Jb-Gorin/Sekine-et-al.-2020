library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DESeq2)

#load RDS
lib.integrated <- readRDS(file = "./rds/all_integ_sct_subset_umap.rds")

# Identify cutoffs for Tox+ & Tcf7+ cells based on normalised expression (after integration)
DefaultAssay(lib.integrated) <- "integrated"
plots <- VlnPlot(lib.integrated, features = c("TOX"), pt.size = 0.001, ncol= 2, group.by = 'donor.ID') + geom_hline(yintercept= 0) 
plots & theme(axis.text.x = element_text(angle = 0, hjust=0.5))
plots <- VlnPlot(lib.integrated, features = c("TCF7"), pt.size = 0.001, ncol= 2, group.by = 'donor.ID') + geom_hline(yintercept= -0.65) 
plots & theme(axis.text.x = element_text(angle = 0, hjust=0.5))

## TOX+ TCF1- cells
tox_pos_tcf7_neg <- subset(lib.integrated, subset = TOX > 0 & TCF7 < -0.65)
tox_pos_tcf7_neg[["cell.ID"]] = "Tox+"
## TOX- TCF7+ cells
tox_neg_tcf7_pos <- subset(lib.integrated, subset = TOX < 0 & TCF7 > -0.65)
tox_neg_tcf7_pos[["cell.ID"]] = "Tcf7+"
##TOX+ TCF7+ cells

# Merge Tox+ and Tcf7+ cells in a new Seurat object
tox_tcf_cells <- merge(tox_pos_tcf7_neg, tox_neg_tcf7_pos)
Idents(tox_tcf_cells) <- "cell.ID"
table(tox_tcf_cells@meta.data[["cell.ID"]])

## DESeq2
DefaultAssay(tox_tcf_cells) <- "SCT"
DESeq2.markers <- FindAllMarkers(tox_tcf_cells, test.use= "DESeq2", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(DESeq2.markers, file= "./csv/DESeq2_markers.csv")
top30_DESeq2 <- DESeq2.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)
DoHeatmap(tox_tcf_cells, features = top30_DESeq2$gene) + NoLegend()

#Visualisation
DefaultAssay(tox_tcf_cells) <- "integrated"
features1 <- c("CCR7", "IL7R", "SELL", "NELL2", "LEF1", "BACH2", "MYC", "ID3", "PDCD1", "TIGIT", "CD244", "LAG3")
features2 <- c("GZMB", "GZMA", "GZMM", "PRF1", "GNLY", "CX3CR1", "TBX21", "EOMES", "ZEB2", "PRDM1", "KLRD1", "S100A4")

## Violin plots
plots <- VlnPlot(tox_tcf_cells, features = features1, pt.size = 0, ncol= 4, group.by = "cell.ID")
plots & theme(axis.text.x = element_text(angle = 0, hjust=0.5))
plots <- VlnPlot(tox_tcf_cells, features = features2, pt.size = 0, ncol= 4, group.by = "cell.ID")
plots & theme(axis.text.x = element_text(angle = 0, hjust=0.5))

