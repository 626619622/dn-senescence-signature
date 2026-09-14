# Step 19: manuscript figures 1 and 3, plus the cohort characteristics table
suppressPackageStartupMessages({
  library(ggplot2); library(reshape2); library(data.table)
})
dir.create("codex_output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/results", showWarnings = FALSE, recursive = TRUE)
dir.create("codex_output/manuscript", showWarnings = FALSE, recursive = TRUE)
log_con <- file("codex_output/logs/r-21-figures.log", open = "wt")
sink(log_con, split = TRUE)
cat("=== step 19: manuscript figures ===\n")
cat("date:", format(Sys.time()), "\n\n")

########################################################################
# Figure 1: study design
########################################################################
box <- function(xmin, xmax, ymin, ymax, label, fill, size = 3.0, col = "grey25") {
  list(rect = data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
       text = data.frame(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label, size = size,
                         col = col))
}
b <- list(
  box(2, 98, 88, 97,
      "Five bulk kidney transcriptomic cohorts\nTraining: GSE96804 (41 DKD / 20 controls)  |  Validation: GSE30528 (9 / 13), GSE99325 (18 / 7 / 115 other nephropathies)\nSpecificity cohort: GSE104954 (17 DKD / 173 other nephropathies / 5 controls)",
      "#DEEBF7"),
  box(2, 98, 76, 84,
      "Differential expression and six MSigDB senescence gene sets\n1,575 differentially expressed genes  \u2192  186 senescence-related genes",
      "#DEEBF7"),
  box(2, 98, 64, 72,
      "Cross-cohort replication\n87 genes with consistent direction in all three cohorts",
      "#DEEBF7"),
  box(2, 31, 40, 58,
      "Single-nucleus localisation\n\nGSE195460\n40,869 nuclei\n5 control / 3 DKD\n\nReplication: GSE131882\n26,097 nuclei\n3 / 3\n\n33/87 and 30/87 genes\nfibroblast-dominant",
      "#E2F0D9"),
  box(35, 68, 40, 58,
      "Machine learning\n\n4 algorithms \u2192 23-gene consensus\nExternal AUC 0.79-0.85 vs controls\n\nSpecificity test: AUC 0.495\n\u2193 specificity filtering\n49-gene score\n\nIndependent validation: AUC 0.870",
      "#FCE4D6"),
  box(72, 98, 40, 58,
      "Genetic analysis\n\n32/49 genes instrumented\nITPR3: OR 0.81, FDR = 2.9e-5\n\nColocalisation: PP.H4 = 0.000\n(PP.H3 = 1.000)\nHLA region\n\u2192 not interpreted as causal",
      "#EDEDED"),
  box(2, 98, 24, 34,
      "Conclusion\nSenescence-related programmes are cell-type specific rather than globally activated, and coupled to fibroblast\nextracellular matrix remodelling. Signatures trained against healthy kidney detect diseased kidney, not DKD.",
      "#FFF2CC", size = 3.1)
)
rects <- do.call(rbind, lapply(b, `[[`, "rect"))
texts <- do.call(rbind, lapply(b, `[[`, "text"))
arrows <- data.frame(
  x = c(50, 50, 50, 50, 16.5, 51.5, 85),
  xend = c(50, 50, 16.5, 51.5, 16.5, 51.5, 85),
  y = c(88, 76, 64, 64, 58, 58, 58),
  yend = c(84, 72, 58, 58, 40, 40, 40)
)
arrows2 <- data.frame(x = c(16.5, 51.5, 85), xend = c(16.5, 51.5, 85),
                      y = c(40, 40, 40), yend = c(34, 34, 34))
p1 <- ggplot() +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = "grey35", linewidth = 0.45) +
  scale_fill_identity() +
  geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
               linewidth = 0.45, color = "grey35") +
  geom_segment(data = arrows2, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
               linewidth = 0.45, color = "grey35") +
  geom_text(data = texts, aes(x = x, y = y, label = label, size = size, color = col),
            lineheight = 1.15) +
  scale_size_identity() + scale_color_identity() +
  coord_cartesian(xlim = c(0, 100), ylim = c(20, 100), expand = FALSE) +
  theme_void()
ggsave("codex_output/figures/fig1-study-design.pdf", p1, width = 11, height = 8.6)
ggsave("codex_output/figures/fig1-study-design.png", p1, width = 11, height = 8.6, dpi = 300)
cat("figure 1 saved\n")

########################################################################
# Figure 3: canonical senescence markers
########################################################################
cm <- read.csv("codex_output/results/canonical-senescence-markers.csv", stringsAsFactors = FALSE)
expected <- c(CDKN1A = "up", CDKN2A = "up", CDKN2B = "up", TP53 = "up", LMNB1 = "down",
              GLB1 = "up", MKI67 = "down", SERPINE1 = "up", IL6 = "up", CXCL8 = "up",
              HMGB1 = "up", GDF15 = "up", FN1 = "up", SIRT1 = "up", SIRT6 = "up",
              ATM = "up", TERF2 = "down", RB1 = "down", E2F1 = "down", CCNA2 = "down",
              CCND2 = "down", PTEN = "up")
cm$expected <- expected[cm$gene]
cm <- cm[!is.na(cm$expected), ]
cm$dataset <- factor(cm$dataset, levels = c("GSE96804", "GSE30528", "GSE99325"))
cm$agree <- ifelse(cm$expected == "up", cm$logFC > 0, cm$logFC < 0)
cm$star <- ifelse(cm$p_value < 0.001, "***", ifelse(cm$p_value < 0.01, "**",
                 ifelse(cm$p_value < 0.05, "*", "")))
# order genes: those agreeing with expectation first
ordg <- aggregate(agree ~ gene, data = cm, FUN = function(x) sum(x, na.rm = TRUE))
ordg$expected <- expected[ordg$gene]
ordg <- ordg[order(-ordg$agree, ordg$gene), ]
cm$gene <- factor(cm$gene, levels = rev(ordg$gene))
cm$logFC_plot <- pmax(pmin(cm$logFC, 1.5), -1.5)
p3 <- ggplot(cm, aes(x = dataset, y = gene)) +
  geom_tile(aes(fill = logFC_plot), color = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 3.2, vjust = 0.78) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, limits = c(-1.5, 1.5),
                       name = expression(log[2]~"FC")) +
  labs(x = NULL, y = NULL,
       title = "Canonical senescence markers across three cohorts",
       subtitle = "* P<0.05, ** P<0.01, *** P<0.001 (uncorrected)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("codex_output/figures/fig3-canonical-markers.pdf", p3, width = 5.4, height = 6.4)
ggsave("codex_output/figures/fig3-canonical-markers.png", p3, width = 5.4, height = 6.4, dpi = 300)
cat("figure 3 saved\n")

########################################################################
# Table 1: cohort characteristics
########################################################################
tab <- data.frame(
  Cohort = c("GSE96804", "GSE30528", "GSE99325", "GSE104954", "GSE195460", "GSE131882"),
  Platform = c("Affymetrix HTA 2.0 (GPL17586)", "Affymetrix HG-U133A 2.0 (GPL571)",
               "Affymetrix HG-U133A / HG-U133Plus2, ENTREZG CDF (GPL19184 / GPL19109)",
               "Affymetrix HG-U133A / HG-U133Plus2 (GPL24120 / GPL22945)",
               "10x Genomics single-nucleus RNA", "Drop-seq single-nucleus RNA"),
  Tissue = c("Glomeruli", "Glomeruli", "Tubulointerstitium", "Tubulointerstitium",
             "Kidney cortex", "Kidney"),
  `DKD / controls` = c("41 / 20", "9 / 13", "18 / 7", "17 / 5", "3 / 5", "3 / 3"),
  `Other nephropathies` = c("-", "-", "115", "173", "-", "-"),
  `Units` = c("Samples", "Samples", "Samples", "Samples", "Nuclei (40,869)", "Nuclei (26,097)"),
  Role = c("Training / discovery", "External validation (case-control)",
           "Specificity filter and validation", "Independent specificity validation",
           "Single-nucleus discovery", "Single-nucleus replication"),
  check.names = FALSE
)
write.csv(tab, "codex_output/results/table1-cohort-characteristics.csv", row.names = FALSE)
cat("\n=== Table 1 ===\n")
print(tab, row.names = FALSE)
cat("\ndone\n")
sink()
close(log_con)
