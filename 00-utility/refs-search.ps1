# Collect real references for the Introduction and Discussion
$S = Join-Path $env:USERPROFILE ".codex\skills\pubmed\scripts\pubmed.ps1"
if (-not (Test-Path -LiteralPath $S)) { throw "pubmed skill script not found" }
New-Item -ItemType Directory -Force -Path "codex_output\refs" | Out-Null

$queries = [ordered]@{
  "dn-burden"      = "('diabetic kidney disease'[tiab] OR 'diabetic nephropathy'[tiab]) AND ('global burden'[tiab] OR 'prevalence'[tiab] OR 'epidemiology'[tiab])"
  "senescence-dn"  = "('cellular senescence'[tiab] OR 'senescence'[tiab]) AND ('diabetic kidney disease'[tiab] OR 'diabetic nephropathy'[tiab])"
  "senescence-kidney" = "('senescence'[tiab]) AND ('kidney'[tiab]) AND ('fibrosis'[tiab] OR 'fibroblast'[tiab] OR 'tubular'[tiab] OR 'podocyte'[tiab])"
  "bioinfo-signature" = "('diabetic nephropathy'[tiab] OR 'diabetic kidney disease'[tiab]) AND ('bioinformatics'[tiab] OR 'machine learning'[tiab] OR 'diagnostic model'[tiab] OR 'gene signature'[tiab])"
  "singlecell-kidney" = "('single-cell'[tiab] OR 'single cell'[tiab] OR 'single-nucleus'[tiab]) AND ('diabetic kidney'[tiab] OR 'diabetic nephropathy'[tiab] OR 'human kidney'[tiab])"
  "mr-coloc"        = "('Mendelian randomization'[tiab]) AND ('colocalization'[tiab] OR 'colocalisation'[tiab]) AND ('drug target'[tiab] OR 'eQTL'[tiab] OR 'pQTL'[tiab])"
  "mr-pitfall"      = "('Mendelian randomization'[tiab]) AND ('HLA'[tiab] OR 'pleiotropy'[tiab] OR 'linkage disequilibrium'[tiab]) AND ('false positive'[tiab] OR 'pitfall'[tiab] OR 'caveat'[tiab] OR 'interpretation'[tiab])"
  "geo-limitation"  = "('gene expression omnibus'[tiab] OR 'GEO'[tiab]) AND ('biomarker'[tiab] OR 'signature'[tiab]) AND ('reproducibility'[tiab] OR 'batch effect'[tiab] OR 'validation'[tiab] OR 'overfitting'[tiab])"
}

foreach ($k in $queries.Keys) {
  $out = "codex_output\refs\pm-$k.md"
  Write-Output ("=== " + $k)
  powershell -NoProfile -ExecutionPolicy Bypass -File $S -Action search -Query $queries[$k] -YearFrom 2015 -Out $out 2>&1 |
    Select-String -Pattern 'hits|written' | ForEach-Object { "    " + $_.Line }
  Start-Sleep -Milliseconds 400
}
