# Step 17: independent single-nucleus replication in GSE131882
# Wilson et al. PNAS 2019: 3 early diabetic kidney samples vs 3 controls, Drop-seq nuclei.
# Goal: does the cell-type localisation of the senescence signature reproduce?
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(harmony); library(msigdbr)
  library(ggplot2); library(reshape2); library(AnnotationDbi); library(org.Hs.eg.db)
})
dir.create("codex_output/logs", showWarnings = FALSE, recursive = TRUE)
log_con <- file("codex_output/logs/r-18-sc-replication.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 17: single-nucleus replication (GSE131882) ===\n")
cat("date:", format(Sys.time()), "\n\n")
files <- sort(list.files("codex_output/geo/GSE131882", pattern = "\\.dgecounts\\.rds$", full.names = TRUE))
cat("files:", length(files), "\n")
objs <- list(); qc <- list()
for (f in files) {
  nm <- sub("\\.dgecounts\\.rds$", "", basename(f))
  grp <- ifelse(grepl("diabetes", nm, ignore.case = TRUE), "DKD", "Control")
  cat("--- reading", nm, "(", grp, ")\n")
  x <- readRDS(f)
  m <- x$umicount$inex$all
  rm(x); gc(verbose = FALSE)
  cat("    matrix:", nrow(m), "genes x", ncol(m), "nuclei\n")
  # this dataset stores Ensembl gene ids; convert them to symbols for the marker panels
  ens <- sub("\\..*$", "", rownames(m))
  sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(ens), column = "SYMBOL",
                               keytype = "ENSEMBL", multiVals = "first")
  rn <- unname(sym[match(ens, names(sym))])
  keep_g <- !is.na(rn) & rn != ""
  m <- m[keep_g, , drop = FALSE]
  rn <- rn[keep_g]
  ord <- order(rowSums(m), decreasing = TRUE)
  m <- m[ord, , drop = FALSE]; rn <- rn[ord]
  dedup <- !duplicated(rn)
  m <- m[dedup, , drop = FALSE]; rownames(m) <- rn[dedup]
  cat("    mapped to", nrow(m), "gene symbols\n")
  o <- CreateSeuratObject(counts = m, project = nm, min.cells = 3, min.features = 200)
  o$sample <- nm; o$group <- grp
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  qc[[length(qc) + 1]] <- data.frame(sample = nm, group = grp, nuclei = ncol(o),
                                     median_features = round(median(o$nFeature_RNA), 0),
                                     median_mt = round(median(o$percent.mt), 2))
  objs[[nm]] <- o
  rm(m, o); gc(verbose = FALSE)
}
qcdf <- do.call(rbind, qc)
cat("\n=== per-sample QC ===\n"); print(qcdf, row.names = FALSE)
write.csv(qcdf, "codex_output/results/sc2-qc-per-sample.csv", row.names = FALSE)
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
obj <- JoinLayers(obj)
cat("\nmerged:", ncol(obj), "nuclei\n"); print(table(obj$group))
keep <- obj$nFeature_RNA > 300 & obj$nFeature_RNA < 6000 & obj$percent.mt < 5
cat("passing QC:", sum(keep), "\n")
obj <- subset(obj, cells = colnames(obj)[keep])
print(table(obj$group))
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2500, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 40, verbose = FALSE)
set.seed(2026)
obj <- RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
cat("\nclusters:", length(levels(obj$seurat_clusters)), "\n")
panel <- list(
  Podocyte = c("NPHS1","NPHS2","PODXL","WT1"),
  ProximalTubule = c("SLC34A1","LRP2","CUBN","SLC22A6","ALDOB"),
  ThickAscLimb = c("SLC12A1","UMOD","CASR"),
  DistalTubule = c("SLC12A3","CALB1","PVALB"),
  CollectingDuct = c("AQP2","AQP3","SCNN1G","SLC4A1"),
  Endothelium = c("PECAM1","EMCN","FLT1","CDH5"),
  Mesangium = c("PDGFRB","RGS5","ACTA2","ITGA8"),
  Fibroblast = c("COL1A1","DCN","LUM","PDGFRA"),
  Myeloid = c("PTPRC","CD68","LYZ","C1QA"),
  Tcell = c("CD3E","CD3D","IL7R"),
  Bcell = c("MS4A1","CD79A","MZB1"),
  NKcell = c("NKG7","GNLY","KLRD1"),
  Cycling = c("MKI67","TOP2A","CCNA2")
)
genes_p <- unique(unlist(panel)); genes_p <- genes_p[genes_p %in% rownames(obj)]
avg <- as.data.frame(AverageExpression(obj, features = genes_p, group.by = "seurat_clusters",
                                       assays = "RNA", layer = "data")$RNA)
z <- t(scale(t(as.matrix(avg))))
score <- sapply(panel, function(g) {
  g <- intersect(g, rownames(z))
  if (!length(g)) return(rep(NA_real_, ncol(z)))
  colMeans(z[g, , drop = FALSE], na.rm = TRUE)
})
assign <- colnames(score)[apply(score, 1, which.max)]
names(assign) <- sub("^g", "", rownames(score))
cat("\n=== cluster assignment ===\n")
print(data.frame(cluster = names(assign), celltype = assign,
                 nuclei = as.integer(table(obj$seurat_clusters)[names(assign)]), row.names = NULL))
write.csv(data.frame(cluster = names(assign), celltype = assign),
          "codex_output/results/sc2-cluster-assignment.csv", row.names = FALSE)
obj$celltype <- factor(unname(assign[as.character(obj$seurat_clusters)]),
                       levels = unique(unname(assign)))
print(table(obj$celltype, obj$group))
saveRDS(obj, "codex_output/data/GSE131882-sc-clustered.rds")
cat("\nclustering done\n")
sink()
close(log_con)
