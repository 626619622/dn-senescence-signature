# Retry the Bioconductor packages against the official repository (the mirror lacked 3.22)

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://bioconductor.org")
options(timeout = 1800)

need <- c("sva", "GSVA", "ConsensusClusterPlus", "impute",
          "hgu133a2.db", "hta20sttranscriptcluster.db")

for (p in need) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("[skip] %s %s\n", p, as.character(utils::packageVersion(p))))
    next
  }
  cat(sprintf("[install] %s\n", p))
  tryCatch(
    BiocManager::install(p, update = FALSE, ask = FALSE, checkBuilt = FALSE),
    error = function(e) cat(sprintf("[FAIL] %s : %s\n", p, conditionMessage(e)))
  )
  cat(if (requireNamespace(p, quietly = TRUE)) "  -> ok\n" else "  -> still missing\n")
}

cat("\n=== re-check ===\n")
allp <- c(need, "WGCNA")
for (p in allp) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-30s %s\n", p, if (ok) as.character(utils::packageVersion(p)) else "MISSING"))
}

cat("\n=== WGCNA load test ===\n")
res <- tryCatch({ suppressPackageStartupMessages(library(WGCNA)); "ok" },
                error = function(e) paste("FAIL:", conditionMessage(e)))
cat(res, "\n")
