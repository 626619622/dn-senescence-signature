# Small pilot: can these CEL files be RMA-normalised with the custom ENTREZG CDF in the archive?

suppressPackageStartupMessages({ library(affy) })

cat("=== installed CDF packages ===\n")
ip <- rownames(installed.packages())
print(ip[grepl("cdf$", ip) & grepl("hgu|hta|clariom", ip, ignore.case = TRUE)])

meta <- read.csv("codex_output/geo/GSE104954/GSE104954-samples.csv", stringsAsFactors = FALSE)
cat("\nsamples:", nrow(meta), "| platforms:\n"); print(table(meta$platform))

cel_dir <- "codex_output/geo/GSE104954/cel"
cels <- list.files(cel_dir, pattern = "\\.CEL\\.gz$", full.names = FALSE)
cat("\nCEL files:", length(cels), "\n")
gsm <- sub("_.*$", "", cels)
plat_of <- meta$platform[match(gsm, meta$gsm)]
cat("CEL files per platform:\n"); print(table(plat_of, useNA = "ifany"))

dn <- meta$gsm[grepl("Diabetic nephropathy", meta$characteristics)]
cat("\ndiabetic nephropathy samples:", length(dn), "\n")
cat("platform of DN samples:\n")
print(table(meta$platform[match(dn, meta$gsm)]))

cdf_gz <- file.path(cel_dir, "GPL24120_HGU133A_Hs_ENTREZG.cdf.gz")
cdf_plain <- file.path(cel_dir, "GPL24120_HGU133A_Hs_ENTREZG.cdf")
if (!file.exists(cdf_plain)) {
  con_in <- gzfile(cdf_gz, "rb"); con_out <- file(cdf_plain, "wb")
  repeat {
    chunk <- readBin(con_in, "raw", n = 1e6)
    if (length(chunk) == 0) break
    writeBin(chunk, con_out)
  }
  close(con_in); close(con_out)
  cat("\ncdf uncompressed:", round(file.info(cdf_plain)$size / 1e6, 1), "MB\n")
}

cat("\nbuilding CDF environment...\n")
make.cdf.env(basename(cdf_plain), envname = "cdf24120")
cat("CDF feature count:", length(ls(envir = get("cdf24120", envir = .GlobalEnv))), "\n")

pilot <- cels[which(plat_of == "GPL24120")[1:3]]
cat("\npilot arrays:", paste(pilot, collapse = ", "), "\n")
eset <- justRMA(filenames = file.path(cel_dir, pilot), cdfname = "cdf24120", verbose = FALSE)
cat("pilot RMA ok:", nrow(eset), "probe sets x", ncol(eset), "arrays\n")
cat("value range:", paste(round(range(exprs(eset)), 2), collapse = " - "), "\n")
cat("example feature names:", paste(head(featureNames(eset), 5), collapse = ", "), "\n")
