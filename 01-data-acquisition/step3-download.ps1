# Step 3a: download bulk series matrices + small single-cell supplementary files

$ErrorActionPreference = 'Continue'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$geo = Join-Path $env:USERPROFILE ".codex\skills\geo-search\scripts\geo.ps1"
$H = @{ 'User-Agent' = 'codex-geo-skill' }
$logDir = "codex_output\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$bulk = @('GSE96804','GSE30528','GSE30529','GSE30122','GSE142153','GSE316126','GSE294519','GSE175759','GSE189005')

Write-Output "=== bulk series matrices ==="
foreach ($a in $bulk) {
    $log = Join-Path $logDir "$a-matrix-dl.log"
    Start-Process -FilePath "powershell" -NoNewWindow -Wait `
        -ArgumentList @("-NoProfile","-ExecutionPolicy","Bypass","-File",$geo,"-Action","matrix","-Accession",$a,"-Download") `
        -RedirectStandardOutput $log | Out-Null
    $ok = (Get-Content $log | Where-Object { $_ -match 'downloaded' }).Count
    Write-Output ("$a : files downloaded = " + $ok)
}

Write-Output ""
Write-Output "=== single-cell supplementary (selected small files) ==="

# GSE142025 full supplementary archive (small)
$targets = @(
    @('GSE142025','GSE142025_RAW.tar')
)

# GSE195460: only the gene-expression filtered matrices (skip large fragment files)
$gse = 'GSE195460'
$base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/$gse/suppl/"
$html = (Invoke-WebRequest $base -Headers $H -UseBasicParsing -TimeoutSec 60).Content
$h5 = [regex]::Matches($html, 'href="(GSE195460_[^"]*filtered_feature_bc_matrix\.h5)"') |
      ForEach-Object { $_.Groups[1].Value } | Select-Object -Unique | Sort-Object
foreach ($f in $h5) { $targets += ,@($gse, $f) }

foreach ($t in $targets) {
    $a = $t[0]; $f = $t[1]
    $dir = "codex_output\geo\$a"
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $out = Join-Path $dir $f
    if (Test-Path -LiteralPath $out) { Write-Output ("skip (exists): " + $f); continue }
    $url = "https://ftp.ncbi.nlm.nih.gov/geo/series/" + ($a -replace '\d{3}$','nnn') + "/$a/suppl/$f"
    try {
        Invoke-WebRequest -Uri $url -Headers $H -OutFile $out -UseBasicParsing -TimeoutSec 1800
        $mb = [math]::Round((Get-Item -LiteralPath $out).Length / 1MB, 1)
        Write-Output ("downloaded: $a/$f ($mb MB)")
    } catch {
        Write-Output ("FAILED: $a/$f -> " + $_.Exception.Message)
    }
    Start-Sleep -Milliseconds 300
}

Write-Output ""
Write-Output "=== downloaded file inventory ==="
Get-ChildItem -Path "codex_output\geo" -Recurse -File |
    Where-Object { $_.Extension -in '.gz','.tar','.h5','.txt' } |
    Sort-Object FullName |
    ForEach-Object { Write-Output ("{0,10:N1} MB  {1}" -f ($_.Length/1MB), $_.FullName.Replace((Get-Location).Path + '\', '')) }
