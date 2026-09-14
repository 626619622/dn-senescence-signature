# Step 7b: normalisation, integration, clustering and cell-type annotation of GSE195460

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-08-sc-cluster.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 7b: cluster and annotate GSE195460 ===\n")
cat("date:", format(Sys.time()), "\n\n")

obj <- readRDS("codex_output/data/GSE195460-sc-merged.rds")
if (length(Layers(obj[["RNA"]])) > 1) obj <- JoinLayers(obj)
cat("cells:", ncol(obj), "| genes:", nrow(obj), "\n")

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2500, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 40, verbose = FALSE)

set.seed(2026)
obj <- RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)

cat("clusters:", length(levels(obj$seurat_clusters)), "\n")
print(table(obj$seurat_clusters))
cat("\ncells per group:\n"); print(table(obj$group))

# --- kidney cell-type marker panel ---------------------------------------
panel <- list(
  Podocyte        = c("NPHS1","NPHS2","PODXL","WT1","SYNPO"),
  ProximalTubule  = c("SLC34A1","LRP2","CUBN","SLC22A6","ALDOB","SLC5A2"),
  ThickAscLimb    = c("SLC12A1","UMOD","CASR"),
  DistalTubule    = c("SLC12A3","CALB1","PVALB"),
  CollectingDuct  = c("AQP2","AQP3","SCNN1G","SLC4A1","ATP6V1G3"),
  Endothelium     = c("PECAM1","EMCN","FLT1","CDH5"),
  Mesangium       = c("PDGFRB","RGS5","ACTA2","ITGA8"),
  Fibroblast      = c("COL1A1","DCN","LUM","PDGFRA"),
  Myeloid         = c("PTPRC","CD68","LYZ","C1QA"),
  Tcell           = c("CD3E","CD3D","IL7R"),
  Bcell           = c("MS4A1","IGHG1","MZB1"),
  NKcell          = c("NKG7","GNLY","KLRD1"),
  Mast            = c("TPSAB1","CPA3","MS4A2"),
  Cycling         = c("MKI67","TOP2A","CCNA2")
)
genes <- unique(unlist(panel))
genes <- genes[genes %in% rownames(obj)]

avg <- AverageExpression(obj, features = genes, group.by = "seurat_clusters",
                         assays = "RNA", slot = "data")$RNA
avg <- as.data.frame(avg)
write.csv(avg, "codex_output/results/sc-cluster-average-expression.csv")

# score each cluster for each cell-type panel (mean z-scored average expression)
z <- t(scale(t(as.matrix(avg))))
score <- sapply(panel, function(g) {
  g <- intersect(g, rownames(z))
  if (length(g) == 0) return(rep(NA_real_, ncol(z)))
  colMeans(z[g, , drop = FALSE], na.rm = TRUE)
})
score <- as.data.frame(score)
score$cluster <- rownames(score)
write.csv(score, "codex_output/results/sc-cluster-celltype-scores.csv", row.names = FALSE)

cat("\n=== cell-type score per cluster (rows = clusters) ===\n")
print(round(score[, setdiff(colnames(score), "cluster")], 2))

panel_cols <- setdiff(colnames(score), "cluster")
assign <- colnames(score)[apply(score[, panel_cols], 1, which.max)]
names(assign) <- score$cluster
# AverageExpression prefixes numeric cluster names with "g"; restore the real names
cell_counts <- table(obj$seurat_clusters)
names(assign) <- sub("^g", "", names(assign))
cat("\n=== top scoring cell type per cluster ===\n")
print(data.frame(cluster = names(assign), assigned = assign,
                 cells = as.integer(cell_counts[names(assign)]),
                 row.names = NULL))
write.csv(data.frame(cluster = names(assign), assigned = assign,
                     cells = as.integer(cell_counts[names(assign)])),
          "codex_output/results/sc-cluster-assignment.csv", row.names = FALSE)

# --- figures -------------------------------------------------------------
pdf("codex_output/figures/fig-sc-umap-clusters.pdf", width = 7, height = 6)
print(DimPlot(obj, reduction = "umap", label = TRUE, raster = FALSE) + ggtitle("Clusters"))
dev.off()
png("codex_output/figures/fig-sc-umap-clusters.png", width = 7, height = 6, units = "in",
    res = 300)
print(DimPlot(obj, reduction = "umap", label = TRUE, raster = FALSE) + ggtitle("Clusters"))
dev.off()

p2 <- DimPlot(obj, reduction = "umap", group.by = "group", raster = FALSE,
              cols = c("Control" = "#9ecae1", "DKD" = "#de2d26")) + ggtitle("Condition")
p3 <- DotPlot(obj, features = genes, group.by = "seurat_clusters") +
  RotatedAxis() + theme(axis.text.x = element_text(size = 7))
ggsave("codex_output/figures/fig-sc-umap-condition.png", p2, width = 7, height = 6, dpi = 300)
ggsave("codex_output/figures/fig-sc-marker-dotplot.png", p3, width = 14, height = 6, dpi = 300)

saveRDS(obj, "codex_output/data/GSE195460-sc-clustered.rds")
saveRDS(assign, "codex_output/data/GSE195460-cluster-labels.rds")
cat("\nsaved clustered object\n")
cat("done\n")
sink()
close(log_con)
