# Install only the packages missing from the local library, without touching existing ones.
# update = FALSE everywhere: never upgrade packages that are already working.

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(timeout = 1800)
options(Ncpus = max(1, parallel::detectCores() - 1))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", update = FALSE)
}

cran_pkgs <- c(
  "glmnet", "randomForest", "pROC", "caret", "Boruta", "xgboost",
  "WGCNA", "forestplot", "ggrepel"
)

bioc_pkgs <- c(
  "sva", "GSVA", "ConsensusClusterPlus", "hgu133a2.db",
  "hta20sttranscriptcluster.db", "preprocessCore"
)

install_if_missing <- function(pkgs, installer) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("[skip]  %-28s %s\n", p, as.character(utils::packageVersion(p))))
      next
    }
    cat(sprintf("[install] %s\n", p))
    tryCatch(
      installer(p, update = FALSE, ask = FALSE),
      error = function(e) cat(sprintf("[FAIL]  %s : %s\n", p, conditionMessage(e)))
    )
    if (requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("  -> ok %s\n", as.character(utils::packageVersion(p))))
    } else {
      cat(sprintf("  -> still missing: %s\n", p))
    }
  }
}

cat("=== CRAN packages ===\n")
install_if_missing(cran_pkgs, function(p, update, ask) install.packages(p, update = update))

cat("\n=== Bioconductor packages ===\n")
install_if_missing(bioc_pkgs, function(p, update, ask) BiocManager::install(p, update = update, ask = ask))

cat("\n=== final check ===\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
ok <- vapply(all_pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
print(data.frame(package = all_pkgs, installed = ok), row.names = FALSE)
cat("\nstill missing:", paste(all_pkgs[!ok], collapse = ", "), "\n")

cat("\npackage library size now:", length(list.files(file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", "4.5"))), "\n")
