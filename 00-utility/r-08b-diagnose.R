# Diagnostic: cluster naming and marker availability in the annotated object

suppressPackageStartupMessages(library(Seurat))
o <- readRDS("codex_output/data/GSE195460-sc-clustered.rds")

cat("cluster level examples:", paste(head(levels(o$seurat_clusters), 5), collapse = ", "), "\n")
cat("number of clusters:", length(levels(o$seurat_clusters)), "\n")
cat("cells per cluster:\n")
print(table(o$seurat_clusters))

cat("\n--- marker presence in the object ---\n")
g <- c("CD3E","CD3D","IL7R","TRAC","MS4A1","CD79A","MZB1","IGHG1",
       "TPSAB1","CPA3","MS4A2","KIT",
       "NPHS1","NPHS2","PODXL","WT1",
       "SLC34A1","LRP2","CUBN","SLC22A6","ALDOB",
       "SLC12A1","UMOD","SLC12A3","CALB1",
       "AQP2","AQP3","SLC4A1",
       "PECAM1","EMCN","CDH5",
       "PDGFRB","RGS5","ACTA2",
       "COL1A1","DCN","LUM",
       "PTPRC","CD68","LYZ","C1QA","CD14",
       "NKG7","GNLY","KLRD1",
       "MKI67","TOP2A",
       "CDKN1A","CDKN2A","TP53","LMNB1","GLB1","SERPINE1","IL6","CXCL8","HMGB1","GDF15","FN1")
present <- g %in% rownames(o)
print(data.frame(gene = g, present = present), row.names = FALSE)

cat("\n--- detection counts (number of nuclei with >0 counts) ---\n")
cnt <- GetAssayData(o, layer = "counts")
for (x in g[present]) cat(sprintf("%-10s %6d cells\n", x, sum(cnt[x, ] > 0)))
