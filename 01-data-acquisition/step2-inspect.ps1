# Inspect source_name / title values for datasets whose group field was not auto-detected

foreach ($acc in @('GSE96804','GSE195460','GSE104954')) {
    Write-Output ("=== " + $acc + " : distinct source_name values")
    $rows = Import-Csv ("codex_output\geo\$acc\$acc-samples.csv")
    $rows | Group-Object source | Sort-Object Count -Descending |
        ForEach-Object { Write-Output ("   " + $_.Count.ToString().PadLeft(3) + "  " + $_.Name) }
    Write-Output ("   -- first 6 titles")
    $rows | Select-Object -First 6 | ForEach-Object { Write-Output ("      " + $_.gsm + " | " + $_.title) }
    Write-Output ""
}
