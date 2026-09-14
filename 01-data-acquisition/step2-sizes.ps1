# Step 2e: estimate download sizes via HTTP HEAD (no actual download)

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$H = @{ 'User-Agent' = 'codex-geo-skill' }

function Get-SeriesDir([string]$a) { $a -replace '\d{3}$', 'nnn' }

function Head-Size([string]$url) {
    try {
        $r = Invoke-WebRequest -Uri $url -Method Head -Headers $H -UseBasicParsing -TimeoutSec 40
        $len = $r.Headers['Content-Length']
        if ($len -is [array]) { $len = $len[0] }
        if ($len) { return [int64]$len }
        return -1
    } catch { return -2 }
}

$bulk = @('GSE96804','GSE30528','GSE30529','GSE30122','GSE104954','GSE47183','GSE99325','GSE99340',
          'GSE294519','GSE175759','GSE142153','GSE189005','GSE316126','GSE145747')

$rows = @()
foreach ($a in $bulk) {
    $base = "https://ftp.ncbi.nlm.nih.gov/geo/series/" + (Get-SeriesDir $a) + "/$a/matrix/"
    $html = (Invoke-WebRequest $base -Headers $H -UseBasicParsing -TimeoutSec 40).Content
    $files = [regex]::Matches($html, 'href="(GSE\d+_series_matrix\.txt\.gz)"') | ForEach-Object { $_.Groups[1].Value } | Select-Object -Unique
    foreach ($f in $files) {
        $size = Head-Size ($base + $f)
        $rows += [pscustomobject]@{ dataset = $a; file = $f; mb = [math]::Round($size / 1MB, 1) }
        Start-Sleep -Milliseconds 350
    }
    if ($files.Count -eq 0) { $rows += [pscustomobject]@{ dataset=$a; file='(no series matrix)'; mb=0 } }
}

Write-Output "=== bulk series matrix sizes ==="
$rows | Format-Table -AutoSize | Out-String -Width 200 | Write-Output

Write-Output "=== single-cell supplementary sizes (selected) ==="
$sc = @(
    @('GSE142025', 'GSE142025_RAW.tar'),
    @('GSE131882', 'GSE131882_RAW.tar'),
    @('GSE183276', 'GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Seurat_03282022.h5Seurat'),
    @('GSE183276', 'GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS.gz'),
    @('GSE195460', 'GSE195460_Control1_filtered_feature_bc_matrix.h5'),
    @('GSE286186', 'GSE286186_RAW.tar')
)
foreach ($s in $sc) {
    $a = $s[0]; $f = $s[1]
    $url = "https://ftp.ncbi.nlm.nih.gov/geo/series/" + (Get-SeriesDir $a) + "/$a/suppl/$f"
    $size = Head-Size $url
    $mb = if ($size -gt 0) { [math]::Round($size / 1MB, 1) } else { $size }
    Write-Output ("{0,-10} {1,-70} {2} MB" -f $a, $f, $mb)
    Start-Sleep -Milliseconds 350
}

$total = ($rows | Where-Object { $_.mb -gt 0 } | Measure-Object mb -Sum).Sum
Write-Output ("bulk series matrix total: ~" + [math]::Round($total, 1) + " MB")
