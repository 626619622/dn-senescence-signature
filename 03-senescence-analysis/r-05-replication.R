# Step 6d: check which senescence-related DEGs replicate in the validation cohorts

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

dir.create("codex_output/results", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)

log_con <- file("codex_output/logs/r-05-replication.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 6d: replication of senescence DEGs ===\n")
cat("date:", format(Sys.time()), "\n\n")

sen_deg <- read.csv("codex_output/results/GSE96804-senescence-DEGs.csv", stringsAsFactors = FALSE)
cat("senescence DEGs from training cohort:", nrow(sen_deg), "\n")

run_limma <- function(expr_path, grp_path, pos) {
  expr <- readRDS(expr_path)
  grp <- read.csv(grp_path, stringsAsFactors = FALSE)
  keep <- grp$group %in% c("Control", pos)
  grp <- grp[keep, , drop = FALSE]
  expr <- expr[, grp$sample, drop = FALSE]
  ok <- !apply(expr, 1, function(x) any(is.na(x)))
  expr <- expr[ok, , drop = FALSE]
  group <- factor(grp$group, levels = c("Control", pos))
  design <- model.matrix(~ group)
  fit <- eBayes(lmFit(expr, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, P.Value = tt$P.Value,
             adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
}

v1 <- run_limma("codex_output/data/GSE30528-expr.rds",
                "codex_output/data/GSE30528-group.csv", "DKD")
v2 <- run_limma("codex_output/data/GSE99325-expr.rds",
                "codex_output/data/GSE99325-group.csv", "DN")
cat("genes testable in GSE30528:", nrow(v1), "| in GSE99325:", nrow(v2), "\n\n")

rep_tab <- sen_deg[, c("gene", "logFC", "adj.P.Val")]
colnames(rep_tab) <- c("gene", "logFC_training", "FDR_training")
rep_tab$logFC_GSE30528 <- v1$logFC[match(rep_tab$gene, v1$gene)]
rep_tab$P_GSE30528     <- v1$P.Value[match(rep_tab$gene, v1$gene)]
rep_tab$logFC_GSE99325 <- v2$logFC[match(rep_tab$gene, v2$gene)]
rep_tab$P_GSE99325     <- v2$P.Value[match(rep_tab$gene, v2$gene)]

rep_tab$consistent_30528 <- with(rep_tab, !is.na(logFC_GSE30528) & sign(logFC_training) == sign(logFC_GSE30528))
rep_tab$consistent_99325 <- with(rep_tab, !is.na(logFC_GSE99325) & sign(logFC_training) == sign(logFC_GSE99325))
rep_tab$consistent_both  <- rep_tab$consistent_30528 & rep_tab$consistent_99325
rep_tab$nominal_30528    <- rep_tab$consistent_30528 & !is.na(rep_tab$P_GSE30528) & rep_tab$P_GSE30528 < 0.05
rep_tab$nominal_99325    <- rep_tab$consistent_99325 & !is.na(rep_tab$P_GSE99325) & rep_tab$P_GSE99325 < 0.05

cat("direction consistent in GSE30528:", sum(rep_tab$consistent_30528), "\n")
cat("direction consistent in GSE99325:", sum(rep_tab$consistent_99325), "\n")
cat("direction consistent in BOTH    :", sum(rep_tab$consistent_both), "\n")
cat("consistent AND p<0.05 in GSE30528:", sum(rep_tab$nominal_30528, na.rm = TRUE), "\n")
cat("consistent AND p<0.05 in GSE99325:", sum(rep_tab$nominal_99325, na.rm = TRUE), "\n")

core <- rep_tab[rep_tab$consistent_both, ]
core <- core[order(core$FDR_training), ]
write.csv(rep_tab, "codex_output/results/senescence-DEG-replication-all.csv", row.names = FALSE)
write.csv(core, "codex_output/results/senescence-core-genes.csv", row.names = FALSE)

cat("\n=== core genes (consistent direction in all three cohorts), top 40 by training FDR ===\n")
print(head(core[, c("gene", "logFC_training", "FDR_training", "logFC_GSE30528", "P_GSE30528",
                    "logFC_GSE99325", "P_GSE99325")], 40), row.names = FALSE)

# scatter of log fold changes, training vs each validation cohort
d1 <- rep_tab[!is.na(rep_tab$logFC_GSE30528), ]
p1 <- ggplot(d1, aes(x = logFC_training, y = logFC_GSE30528)) +
  geom_point(aes(color = consistent_30528), size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5, alpha = 0.15) +
  scale_color_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey65")) +
  labs(x = expression("log"[2]*"FC (GSE96804, training)"),
       y = expression("log"[2]*"FC (GSE30528, validation)"),
       color = "same direction",
       title = "Replication of senescence-related DEGs") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

d2 <- rep_tab[!is.na(rep_tab$logFC_GSE99325), ]
p2 <- ggplot(d2, aes(x = logFC_training, y = logFC_GSE99325)) +
  geom_point(aes(color = consistent_99325), size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5, alpha = 0.15) +
  scale_color_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey65")) +
  labs(x = expression("log"[2]*"FC (GSE96804, training)"),
       y = expression("log"[2]*"FC (GSE99325, validation)"),
       color = "same direction",
       title = "Replication (independent cohort)") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

ggsave("codex_output/figures/fig-replication-scatter.pdf", p1 + p2, width = 11, height = 5)
ggsave("codex_output/figures/fig-replication-scatter.png", p1 + p2, width = 11, height = 5, dpi = 300)

cat("\nPearson r (training vs GSE30528):",
    round(cor(d1$logFC_training, d1$logFC_GSE30528), 3), "\n")
cat("Pearson r (training vs GSE99325):",
    round(cor(d2$logFC_training, d2$logFC_GSE99325), 3), "\n")
cat("\ndone\n")
sink()
close(log_con)
