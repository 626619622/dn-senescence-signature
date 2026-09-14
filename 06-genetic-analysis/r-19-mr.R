# Step 16: cis-eQTL Mendelian randomisation
# Exposure : eQTLGen blood cis-eQTLs (31,684 samples)
# Outcome  : FinnGen R12 DM_NEPHROPATHY (5,579 cases / 90,951 controls)
# Question : do the 49 signature genes show genetic evidence of a causal effect on DKD?

suppressPackageStartupMessages({
  library(data.table)
})

log_con <- file("codex_output/logs/r-19-mr.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 16: cis-eQTL MR ===\n")
cat("date:", format(Sys.time()), "\n\n")

genes <- read.csv("codex_output/results/specificity-filtered-genes.csv",
                  stringsAsFactors = FALSE)$gene
cat("signature genes:", length(genes), "\n\n")

# ---------- 1. exposure: top cis-eQTL per gene ------------------------------
cat("reading eQTLGen cis-eQTLs (this is a large file)...\n")
eq <- fread("codex_output/mr/eqtlgen-cis-eQTL-FDR0.05.txt.gz",
            select = c("SNP","Pvalue","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                       "Zscore","Gene","GeneSymbol","NrSamples","FDR"),
            showProgress = FALSE)
cat("total cis-eQTL rows:", nrow(eq), "\n")

eq_sub <- eq[GeneSymbol %in% genes]
cat("rows for the 49 signature genes:", nrow(eq_sub), "\n")
cat("genes with at least one cis-eQTL:", length(unique(eq_sub$GeneSymbol)), "of", length(genes), "\n")
rm(eq); gc(verbose = FALSE)

eq_sub[, absZ := abs(Zscore)]
top <- eq_sub[order(-absZ), .SD[1], by = GeneSymbol]
setnames(top, c("SNP","AssessedAllele","OtherAllele","Zscore","Pvalue","NrSamples","FDR"),
         c("SNP","ea_exp","oa_exp","z_exp","p_exp","n_exp","fdr_exp"))
cat("\ntop cis-eQTL per gene (first 10):\n")
print(head(top[, .(GeneSymbol, SNP, ea_exp, oa_exp, z_exp, p_exp)], 10), row.names = FALSE)

# ---------- 2. outcome: FinnGen DM_NEPHROPATHY ------------------------------
cat("\nreading FinnGen DM_NEPHROPATHY summary statistics...\n")
fg <- fread("codex_output/mr/finngen_R12_DM_NEPHROPATHY.gz",
            select = c("#chrom","pos","ref","alt","rsids","nearest_genes",
                       "pval","beta","sebeta","af_alt"),
            showProgress = FALSE)
setnames(fg, c("#chrom","ref","alt","rsids","pval","beta","sebeta","af_alt"),
         c("chr","ref","alt","SNP","p_out","b_out","se_out","af_alt"))
cat("FinnGen variants:", nrow(fg), "\n")

fg_sub <- fg[SNP %in% top$SNP]
cat("top cis-eQTL SNPs found in FinnGen:", nrow(fg_sub), "of", nrow(top), "\n")

# ---------- 3. harmonise and compute the Wald ratio --------------------------
dt <- merge(top, fg_sub, by = "SNP")
cat("merged instruments:", nrow(dt), "\n\n")

harmonise <- function(ea, oa, alt, ref, b) {
  same <- (ea == alt & oa == ref)
  flip <- (ea == ref & oa == alt)
  list(keep = same | flip, b = ifelse(flip, -b, b))
}
h <- harmonise(dt$ea_exp, dt$oa_exp, dt$alt, dt$ref, dt$b_out)
dt <- dt[h$keep]; dt$b_out_h <- h$b[h$keep]
cat("instruments after allele harmonisation:", nrow(dt), "\n")

# eQTLGen reports Z-scores; convert to a beta scale using the cohort allele frequency
dt[, maf := pmin(af_alt, 1 - af_alt)]
dt[, b_exp := z_exp / sqrt(2 * maf * (1 - maf) * (n_exp + z_exp^2))]
dt[, b_mr := b_out_h / b_exp]
dt[, se_mr := se_out / abs(b_exp)]
dt[, z_mr := b_mr / se_mr]
dt[, p_mr := 2 * pnorm(-abs(z_mr))]
dt[, OR := exp(b_mr)]
dt[, F_stat := z_exp^2]
dt[, FDR_mr := p.adjust(p_mr, method = "BH")]

setorder(dt, p_mr)
out <- dt[, .(gene = GeneSymbol, snp = SNP, n_exp = n_exp,
              beta_exposure = round(b_exp, 4), z_exposure = round(z_exp, 2),
              F_statistic = round(F_stat, 1),
              beta_MR = round(b_mr, 4), se_MR = round(se_mr, 4),
              OR = round(OR, 3), p_MR = signif(p_mr, 3),
              FDR_MR = signif(FDR_mr, 3))]
cat("\n=== cis-eQTL MR results (all genes with an instrument) ===\n")
print(head(out, 30), row.names = FALSE)
write.csv(out, "codex_output/results/mr-cis-eqtl-results.csv", row.names = FALSE)

cat("\n=== summary ===\n")
cat("genes tested:", nrow(out), "\n")
cat("instruments with F > 10:", sum(out$F_statistic > 10, na.rm = TRUE), "\n")
cat("nominally significant (p < 0.05):", sum(out$p_MR < 0.05, na.rm = TRUE), "\n")
cat("FDR < 0.05:", sum(out$FDR_MR < 0.05, na.rm = TRUE), "\n")
cat("FDR < 0.20:", sum(out$FDR_MR < 0.20, na.rm = TRUE), "\n")

sig <- out[out$p_MR < 0.05, ]
if (nrow(sig)) {
  cat("\nnominally significant genes:\n")
  print(sig, row.names = FALSE)
}

cat("\ndone\n")
sink()
close(log_con)
