# Step 3b: verify downloaded files are real (size + gzip signature) and fetch the missing ones

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$H = @{ 'User-Agent' = 'codex-geo-skill' }
$root = (Get-Location).Path

Write-Output "=== file integrity check ==="
Get-ChildItem "codex_output\geo" -Recurse -File |
    Where-Object { $_.Name -like '*series_matrix*' -or $_.Name -like '*RAW.tar' } |
    Sort-Object FullName | ForEach-Object {
        $fs = [System.IO.File]::OpenRead($_.FullName)
        $buf = New-Object byte[] 4
        $n = $fs.Read($buf, 0, 4)
        $fs.Close()
        $sig = ($buf[0..($n-1)] | ForEach-Object { $_.ToString('x2') }) -join ''
        $isGzip = $sig.StartsWith('1f8b')
        Write-Output ("{0,10:N0} bytes  gzip={1,-5}  {2}" -f $_.Length, $isGzip, $_.FullName.Replace($root + '\',''))
    }

Write-Output ""
Write-Output "=== remaining GSE195460 h5 files on server ==="
$base195 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195460/suppl/"
$html = (Invoke-WebRequest $base195 -Headers $H -UseBasicParsing -TimeoutSec 60).Content
$h5 = [regex]::Matches($html, 'href="(GSE195460_[^"]*filtered_feature_bc_matrix\.h5)"') |
      ForEach-Object { $_.Groups[1].Value } | Select-Object -Unique | Sort-Object
Write-Output ("server h5 count: " + @($h5).Count)

$need = @()
foreach ($f in $h5) {
    if (-not (Test-Path -LiteralPath ("codex_output\geo\GSE195460\" + $f))) { $need += $f }
}
Write-Output ("missing locally: " + @($need).Count)
foreach ($f in $need) {
    $out = "codex_output\geo\GSE195460\$f"
    try {
        Invoke-WebRequest -Uri ($base195 + $f) -Headers $H -OutFile $out -UseBasicParsing -TimeoutSec 1800
        Write-Output ("downloaded: $f (" + [math]::Round((Get-Item -LiteralPath $out).Length/1MB,1) + " MB)")
    } catch { Write-Output ("FAILED: $f -> " + $_.Exception.Message) }
    Start-Sleep -Milliseconds 300
}

Write-Output ""
Write-Output "=== GSE142025 supplementary ==="
$url142 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142025/suppl/GSE142025_RAW.tar"
$out142 = "codex_output\geo\GSE142025\GSE142025_RAW.tar"
New-Item -ItemType Directory -Force -Path "codex_output\geo\GSE142025" | Out-Null
if (-not (Test-Path -LiteralPath $out142)) {
    try {
        Invoke-WebRequest -Uri $url142 -Headers $H -OutFile $out142 -UseBasicParsing -TimeoutSec 1800
        Write-Output ("downloaded: GSE142025_RAW.tar (" + [math]::Round((Get-Item -LiteralPath $out142).Length/1MB,1) + " MB)")
    } catch { Write-Output ("FAILED: " + $_.Exception.Message) }
}
