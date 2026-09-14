# hta20sttranscriptcluster.db was retired from current Bioconductor; install it from the
# Bioconductor 3.13 archive (source package, needs Rtools on PATH for the build step).

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout = 1800)

if (requireNamespace("hta20sttranscriptcluster.db", quietly = TRUE)) {
  cat("already installed:", as.character(utils::packageVersion("hta20sttranscriptcluster.db")), "\n")
} else {
  url <- "https://bioconductor.org/packages/3.13/data/annotation/src/contrib/hta20sttranscriptcluster.db_8.8.0.tar.gz"
  dest <- file.path(tempdir(), "hta20sttranscriptcluster.db_8.8.0.tar.gz")
  cat("downloading:", url, "\n")
  tryCatch({
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    cat("downloaded, size:", round(file.info(dest)$size / 1e6, 2), "MB\n")
    install.packages(dest, repos = NULL, type = "source")
  }, error = function(e) cat("FAIL:", conditionMessage(e), "\n"))
}

ok <- requireNamespace("hta20sttranscriptcluster.db", quietly = TRUE)
cat("installed:", ok, "\n")
if (ok) {
  suppressPackageStartupMessages(library(hta20sttranscriptcluster.db))
  cat("probe count:", length(AnnotationDbi::keys(hta20sttranscriptcluster.db, keytype = "PROBEID")), "\n")
}
cat("Rtools on PATH:", nzchar(Sys.which("make")), "\n")
