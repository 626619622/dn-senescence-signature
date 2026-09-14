# Organise manuscript supplementary files with journal-friendly names
$src = "codex_output\results"
$fig = "codex_output\figures"
$out = "codex_output\supplementary"
New-Item -ItemType Directory -Force -Path $out | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $out "tables") | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $out "figures") | Out-Null

$tables = [ordered]@{
  "TableS1-senescence-related-DEGs.csv"        = "$src\GSE96804-senescence-DEGs.csv"
  "TableS2-senescence-pathway-scores.csv"      = "$src\senescence-score-comparison.csv"
  "TableS3-core-87-gene-replication.csv"       = "$src\senescence-DEG-replication-all.csv"
  "TableS4-cell-vs-sample-level-tests.csv"     = "$src\sc-cell-vs-sample-level.csv"
  "TableS5-feature-selection.csv"              = "$src\ml-feature-selection.csv"
  "TableS6-specificity-filter.csv"             = "$src\specificity-filter-table.csv"
  "TableS7-mendelian-randomisation.csv"        = "$src\mr-cis-eqtl-results.csv"
  "TableS8-colocalisation.csv"                 = "$src\coloc-results.csv"
  "TableS9-canonical-senescence-markers.csv"   = "$src\canonical-senescence-markers.csv"
  "TableS10a-internal-cross-validation.csv"    = "$src\ml-cv-performance.csv"
  "TableS10b-external-validation.csv"          = "$src\ml-external-validation.csv"
  "TableS10c-specificity-test.csv"             = "$src\ml-specificity.csv"
  "TableS10d-independent-validation.csv"       = "$src\independent-validation-auc.csv"
  "TableS10e-final-model-AUC.csv"              = "$src\final-model-auc.csv"
  "TableS10f-recalibration.csv"                = "$src\recalibration-performance.csv"
  "TableS11-risk-score-weights.csv"            = "$src\final-risk-score-weights.csv"
  "TableS12a-single-cell-QC.csv"               = "$src\sc-qc-per-sample.csv"
  "TableS12b-single-cell-composition.csv"      = "$src\sc-composition-per-sample.csv"
  "TableS12c-replication-cohort-QC.csv"        = "$src\sc2-qc-per-sample.csv"
  "TableS13a-cluster-annotation.csv"           = "$src\sc-cluster-assignment.csv"
  "TableS13b-marker-detection.csv"             = "$src\sc-marker-detection-by-celltype.csv"
  "TableS13c-cell-type-score-tests.csv"        = "$src\sc-senescence-score-tests.csv"
  "TableS14a-replication-cell-type-tests.csv"  = "$src\sc2-senescence-tests.csv"
  "TableS14b-replication-direction-agreement.csv" = "$src\sc2-direction-agreement.csv"
  "TableS14c-replication-pseudobulk.csv"       = "$src\sc2-pseudobulk-tests.csv"
  "TableS14d-replication-core-gene-localisation.csv" = "$src\sc2-core-gene-localisation.csv"
}
$n = 0
foreach ($k in $tables.Keys) {
  if (Test-Path -LiteralPath $tables[$k]) {
    Copy-Item -LiteralPath $tables[$k] -Destination (Join-Path $out "tables\$k") -Force
    $n++
  } else { Write-Output ("MISSING: " + $tables[$k]) }
}
Write-Output ("tables copied: " + $n)

$figs = [ordered]@{
  "FigureS1-single-nucleus-UMAP.png"        = "$fig\fig-sc-umap-clusters.png"
  "FigureS2-canonical-markers-single-cell.png" = "$fig\fig-sc-canonical-markers.png"
  "FigureS3-core-genes-single-cell.png"     = "$fig\fig-sc-core-genes.png"
  "FigureS4-pseudobulk-core87.png"          = "$fig\fig-sc-pseudobulk-core87.png"
  "FigureS5-replication-UMAP.png"           = "$fig\fig-sc2-umap-celltypes.png"
  "FigureS6-replication-scores.png"         = "$fig\fig-sc2-senescence-scores.png"
  "FigureS7-feature-importance.png"         = "$fig\fig-ml-importance.png"
  "FigureS8-SHAP-importance.png"            = "$fig\fig-ml-shap.png"
  "FigureS9-specificity-map.png"            = "$fig\fig-specificity-map.png"
  "FigureS10-calibration.png"               = "$fig\fig-risk-calibration.png"
  "FigureS11-decision-curve.png"            = "$fig\fig-risk-dca.png"
  "FigureS12-recalibration-curves.png"      = "$fig\fig-recalibration-curve.png"
  "FigureS13-recalibration-distribution.png" = "$fig\fig-recalibration-distribution.png"
}
$m = 0
foreach ($k in $figs.Keys) {
  if (Test-Path -LiteralPath $figs[$k]) {
    Copy-Item -LiteralPath $figs[$k] -Destination (Join-Path $out "figures\$k") -Force
    $m++
  } else { Write-Output ("MISSING: " + $figs[$k]) }
}
Write-Output ("figures copied: " + $m)
