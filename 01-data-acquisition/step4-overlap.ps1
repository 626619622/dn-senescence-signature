# Check GSM-level sample overlap between candidate cohorts (decisive independence test)

$sets = @{}
foreach ($acc in @('GSE99325','GSE30528','GSE30529','GSE30122','GSE47183','GSE104954','GSE96804','GSE142025','GSE195460','GSE294519','GSE175759')) {
    $p = "codex_output\geo\$acc\$acc-samples.csv"
    if (Test-Path -LiteralPath $p) {
        $sets[$acc] = (Import-Csv $p | ForEach-Object { $_.gsm })
    }
}

Write-Output "=== dataset sizes ==="
foreach ($k in $sets.Keys) { Write-Output ("  {0,-10} {1} samples (first {2} ...)" -f $k, $sets[$k].Count, $sets[$k][0]) }

Write-Output ""
Write-Output "=== pairwise GSM overlap (0 = independent samples) ==="
$keys = $sets.Keys | Sort-Object
for ($i = 0; $i -lt $keys.Count; $i++) {
    for ($j = $i + 1; $j -lt $keys.Count; $j++) {
        $a = $keys[$i]; $b = $keys[$j]
        $inter = @(Compare-Object $sets[$a] $sets[$b] -IncludeEqual -ExcludeDifferent)
        if ($inter.Count -gt 0) {
            Write-Output ("  OVERLAP  {0} vs {1} : {2} shared GSMs" -f $a, $b, $inter.Count)
        }
    }
}
Write-Output "  (pairs not listed share zero samples)"
