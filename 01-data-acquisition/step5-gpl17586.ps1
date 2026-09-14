# Download the GPL17586 (HTA 2.0) annotation table used to map transcript clusters to gene symbols

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$H = @{ 'User-Agent' = 'codex-geo-skill' }
$dir = "codex_output\geo\GPL17586"
New-Item -ItemType Directory -Force -Path $dir | Out-Null

$url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL17nnn/GPL17586/annot/GPL17586.annot.gz"
$out = Join-Path $dir "GPL17586.annot.gz"
if (-not (Test-Path -LiteralPath $out)) {
    Invoke-WebRequest -Uri $url -Headers $H -OutFile $out -UseBasicParsing -TimeoutSec 1800
}
Write-Output ("downloaded: " + $out + " (" + [math]::Round((Get-Item -LiteralPath $out).Length/1MB,1) + " MB)")

Add-Type -AssemblyName System.IO.Compression.FileSystem
$fs = [System.IO.File]::OpenRead($out)
$gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
$sr = New-Object System.IO.StreamReader($gz)
$i = 0
while ($i -lt 40 -and -not $sr.EndOfStream) {
    $l = $sr.ReadLine()
    if ($l.Length -gt 180) { $l = $l.Substring(0,180) + '...' }
    Write-Output ("  " + $i.ToString().PadLeft(2) + ": " + $l)
    $i++
}
$sr.Close(); $gz.Close(); $fs.Close()
