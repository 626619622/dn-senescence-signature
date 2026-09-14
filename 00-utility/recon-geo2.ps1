# Verify a second batch: single-cell / islet / additional kidney datasets

$accs = @(
    'GSE84133','GSE81547','GSE86469','GSE101207','GSE114297','GSE124742','GSE154126',
    'GSE151302','GSE140989','GSE195460','GSE183276','GSE171406','GSE142025','GSE173343',
    'GSE222512','GSE221156','GSE202965','GSE220243'
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
                accession = $a; organism = $d.taxon; n_samples = $d.n_samples; type = $d.gdsType
                platform = ($d.gpl -split ';')[0]; date = $d.pdat; title = $d.title
            }
        } else {
            $rows += [pscustomobject]@{ accession = $a; organism='NOT FOUND'; n_samples=''; type=''; platform=''; date=''; title='' }
        }
    } catch {
        $rows += [pscustomobject]@{ accession = $a; organism='ERROR'; n_samples=''; type=''; platform=''; date=''; title=$_.Exception.Message }
    }
    Start-Sleep -Milliseconds 400
}

$out = "codex_output\recon\geo-candidates-batch2.csv"
$rows | Export-Csv -Path $out -NoTypeInformation -Encoding UTF8
$rows | Format-Table -AutoSize | Out-String -Width 250 | Write-Output
Write-Output ("saved: " + $out)
