# Smoke test: R version, library paths, availability of packages needed for the pipeline

cat("R version:", R.version.string, "\n")
cat("locale:", Sys.getlocale("LC_CTYPE"), "\n\n")

cat("--- .libPaths()\n")
print(.libPaths())

lib_user <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5")
cat("\npersonal library path:", lib_user, "\n")
cat("personal library exists:", dir.exists(lib_user), "\n")
if (dir.exists(lib_user)) {
  cat("packages installed there:", length(list.files(lib_user)), "\n")
}

pkgs <- c(
  # data download / preprocessing
  "GEOquery", "limma", "sva", "preprocessCore", "affy", "oligo",
  # annotation
  "hgu133a2.db", "hta20sttranscriptcluster.db", "org.Hs.eg.db", "AnnotationDbi",
  # network / enrichment
  "WGCNA", "clusterProfiler", "GSVA", "ConsensusClusterPlus",
  # machine learning
  "glmnet", "randomForest", "e1071", "caret", "pROC", "rms", "Boruta", "xgboost",
  # plotting
  "ggplot2", "pheatmap", "ComplexHeatmap", "ggpubr", "patchwork",
  # single cell
  "Seurat", "harmony", "CellChat"
)

cat("\n--- package availability\n")
res <- data.frame(
  package = pkgs,
  installed = vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1)),
  version = vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
print(res, row.names = FALSE)

cat("\nmissing:", paste(res$package[!res$installed], collapse = ", "), "\n")
