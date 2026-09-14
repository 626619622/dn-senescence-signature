# Step 2c: find additional independent human DKD/DN expression cohorts for external validation

$geo = Join-Path $env:USERPROFILE ".codex\skills\geo-search\scripts\geo.ps1"
$outDir = "codex_output\geo"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$queries = @{
    'title-dn'  = "('diabetic nephropathy'[Title] OR 'diabetic kidney disease'[Title] OR 'diabetic glomerulopathy'[Title] OR 'DKD'[Title])"
    'broad-kidney' = "('diabetic nephropathy'[All Fields] OR 'diabetic kidney disease'[All Fields]) AND ('kidney'[Title] OR 'glomerul*'[Title] OR 'renal'[Title] OR 'tubul*'[Title])"
    'dn-tubular'= "('diabetic'[Title] AND ('tubulointerstitial'[All Fields] OR 'proximal tubule'[All Fields] OR 'podocyte'[All Fields]) AND 'Homo sapiens'[Organism])"
}

$all = @()
foreach ($k in $queries.Keys) {
    $out = Join-Path $outDir ("candidates-" + $k + ".csv")
    & powershell -NoProfile -ExecutionPolicy Bypass -File $geo -Action search -Query $queries[$k] -Organism "Homo sapiens" -Retmax 60 -Out $out | Out-Null
    if (Test-Path -LiteralPath $out) {
        $rows = Import-Csv $out
        Write-Output ("$k : " + $rows.Count + " hits")
        $all += $rows
    }
}

$all = $all | Sort-Object gse -Unique
Write-Output ("unique datasets: " + $all.Count)

$expr = $all | Where-Object { $_.data_type -match 'Expression profiling' -and [int]$_.n_samples -ge 20 }
Write-Output ("expression datasets with n >= 20 : " + $expr.Count)

$expr | Sort-Object { [int]$_.n_samples } -Descending |
    Select-Object gse, n_samples, data_type, platform, date, title |
    Format-Table -AutoSize | Out-String -Width 400 | Write-Output

$expr | Sort-Object gse -Unique | Export-Csv -LiteralPath (Join-Path $outDir "dkd-validation-candidates.csv") -NoTypeInformation -Encoding UTF8
Write-Output ("written: " + (Join-Path $outDir "dkd-validation-candidates.csv"))
