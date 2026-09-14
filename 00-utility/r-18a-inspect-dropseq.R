# Inspect the Drop-seq dgecounts structure for the GSE131882 replication dataset

suppressPackageStartupMessages({ library(Matrix) })

f <- "codex_output/geo/GSE131882/GSM3823944_diabetes.s3.dgecounts.rds"
cat("file:", basename(f), "| size:", round(file.info(f)$size / 1e6, 1), "MB\n")
x <- readRDS(f)

cat("top-level class:", class(x), "\n")
cat("top-level names:", paste(names(x), collapse = ", "), "\n\n")

if (is.list(x)) {
  for (n in names(x)) {
    cat("--- ", n, " (", class(x[[n]])[1], ")\n", sep = "")
    v <- x[[n]]
    if (is.list(v)) {
      cat("    sub-names:", paste(names(v), collapse = ", "), "\n")
      for (m in names(v)) {
        z <- v[[m]]
        if (is.list(z)) {
          cat(sprintf("    %-8s -> names: %s\n", m, paste(names(z), collapse = ", ")))
          for (k in names(z)) {
            cat(sprintf("        %-8s %-10s dim=%s\n", k, class(z[[k]])[1],
                        paste(dim(z[[k]]), collapse = " x ")))
          }
        } else {
          cat(sprintf("    %-8s %-10s dim=%s\n", m, class(z)[1],
                      paste(dim(z), collapse = " x ")))
        }
      }
    } else {
      cat("    dim:", paste(dim(v), collapse = " x "), "\n")
    }
  }
}
