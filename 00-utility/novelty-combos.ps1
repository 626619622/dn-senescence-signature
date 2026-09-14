# Cross-tabulate angle combinations within the 600 DKD/ML records to locate genuine gaps

$csv = "codex_output\recon\dkd-ml-records-2024-2026.csv"
$rows = Import-Csv $csv
Write-Output ("n = " + $rows.Count)

$nameMap = @{
    "ML model"      = 'machine learning|random forest|lasso|xgboost|lightgbm|nomogram|predictive model|diagnostic model'
    "diagnostic"    = 'diagnos'
    "progression"   = 'progress|prognos|predict'
    "single-cell"   = 'single-cell|single cell|scrna|snrna'
    "Mendelian"     = 'mendelian'
    "external-val"  = 'external validation|validation cohort|independent cohort|external cohort'
    "SHAP"          = 'shap|explainab'
    "drug"          = 'drug|compound|docking|cmap'
    "ferroptosis"   = 'ferroptosis'
    "immune"        = 'immune'
    "podocyte"      = 'podocyte'
    "tubular"       = 'tubul'
    "multi-omics"   = 'multi-omics|multiomics|proteomic|transcriptom.*integrat'
    "lactylation"   = 'lactyl'
    "disulfidptosis"= 'disulfidptosis'
    "senescence"    = 'senescence'
    "deep-learning" = 'deep learning|neural network|transformer'
    "subtype"       = 'subtype|consensus cluster|clustering'
    "microbiome"    = 'microbiota|microbiome'
    "mitochondrial" = 'mitochondri'
}

foreach ($r in $rows) {
    $t = ($r.title + " " + $r.angles_in_text).ToLower()
    $flags = [ordered]@{}
    foreach ($k in $nameMap.Keys) { $flags[$k] = ($t -match $nameMap[$k]) }
    $r | Add-Member -NotePropertyName flags -NotePropertyValue $flags -Force
}

function Count-Combo([string[]]$keys) {
    $n = 0
    foreach ($r in $rows) {
        $ok = $true
        foreach ($k in $keys) { if (-not $r.flags[$k]) { $ok = $false; break } }
        if ($ok) { $n++ }
    }
    return $n
}

$combos = @(
    @("single-cell","Mendelian"),
    @("single-cell","Mendelian","drug"),
    @("single-cell","external-val"),
    @("ML model","external-val"),
    @("ML model","external-val","SHAP"),
    @("single-cell","SHAP"),
    @("Mendelian","external-val"),
    @("lactylation","single-cell"),
    @("disulfidptosis","single-cell"),
    @("senescence","single-cell"),
    @("senescence","Mendelian"),
    @("podocyte","single-cell"),
    @("deep-learning","single-cell"),
    @("subtype","single-cell"),
    @("subtype","Mendelian"),
    @("multi-omics","Mendelian"),
    @("podocyte","Mendelian"),
    @("ferroptosis","Mendelian"),
    @("microbiome","Mendelian"),
    @("microbiome","single-cell"),
    @("drug","Mendelian")
)

Write-Output ""
Write-Output "=== Combination frequency (n = 600) ==="
$res = foreach ($c in $combos) {
    [pscustomobject]@{ combo = ($c -join ' + '); count = (Count-Combo $c) }
}
$res | Sort-Object count | Format-Table -AutoSize | Out-String -Width 200 | Write-Output

Write-Output "=== Example titles for the rarest high-value combos ==="
foreach ($c in @(@("single-cell","Mendelian","drug"), @("ML model","external-val","SHAP"), @("senescence","single-cell"), @("subtype","Mendelian"))) {
    Write-Output ("--- " + ($c -join ' + '))
    $k = 0
    foreach ($r in $rows) {
        $ok = $true
        foreach ($x in $c) { if (-not $r.flags[$x]) { $ok = $false; break } }
        if ($ok) {
            Write-Output ("   [" + $r.year + "] " + $r.journal + " | " + $r.title)
            $k++
            if ($k -ge 6) { break }
        }
    }
}

Write-Output ""
Write-Output "=== Journals of the 600 records (top 20) ==="
$rows | Group-Object journal | Sort-Object Count -Descending | Select-Object -First 20 |
    ForEach-Object { Write-Output ("{0,3}  {1}" -f $_.Count, $_.Name) }
