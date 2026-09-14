# Verify candidate GEO datasets (human, diabetes-related) actually exist and have usable sample counts

$accs = @(
    'GSE30528','GSE30529','GSE96804','GSE104954','GSE142153','GSE131882','GSE30122','GSE47183',
    'GSE60436','GSE102485','GSE140959','GSE24290','GSE95849',
    'GSE20966','GSE38642','GSE76894','GSE76895','GSE50397','GSE25462','GSE18732','GSE13070',
    'GSE16415','GSE20950','GSE9006','GSE55633'
)

$rows = @()
foreach ($a in $accs) {
    try {
        $u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmode=json&tool=codex&retmax=5&term=" + [uri]::EscapeDataString($a + '[ACCN]')
        $r = Invoke-RestMethod $u -TimeoutSec 40
        $ids = $r.esearchresult.idlist
        if ($ids.Count -gt 0) {
            $s = Invoke-RestMethod ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&tool=codex&id=" + $ids[0]) -TimeoutSec 40
            $d = $s.result.$($ids[0])
            $rows += [pscustomobject]@{
                accession = $a
                organism  = $d.taxon
                n_samples = $d.n_samples
                type      = $d.gdsType
                platform  = ($d.gpl -split ';')[0]
                date      = $d.pdat
                title     = $d.title
            }
        } else {
            $rows += [pscustomobject]@{ accession = $a; organism='NOT FOUND'; n_samples=''; type=''; platform=''; date=''; title='' }
        }
    } catch {
        $rows += [pscustomobject]@{ accession = $a; organism='ERROR'; n_samples=''; type=''; platform=''; date=''; title=$_.Exception.Message }
    }
    Start-Sleep -Milliseconds 400
}

$out = "codex_output\recon\geo-candidates.csv"
New-Item -ItemType Directory -Force -Path (Split-Path $out) | Out-Null
$rows | Export-Csv -Path $out -NoTypeInformation -Encoding UTF8
$rows | Format-Table -AutoSize | Out-String -Width 250 | Write-Output
Write-Output ("saved: " + $out)
