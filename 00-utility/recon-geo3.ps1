# Full titles for selected datasets + recount check

$accs = @('GSE221156','GSE173343','GSE142025','GSE183276','GSE195460','GSE104954','GSE96804','GSE131882')
foreach ($a in $accs) {
    $u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmode=json&tool=codex&retmax=1&term=" + [uri]::EscapeDataString($a + '[ACCN]')
    $r = Invoke-RestMethod $u -TimeoutSec 40
    $id = $r.esearchresult.idlist[0]
    $s = Invoke-RestMethod ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&tool=codex&id=" + $id) -TimeoutSec 40
    $d = $s.result.$id
    Write-Output ("{0} | n={1} | {2} | {3}" -f $a, $d.n_samples, $d.summary[0..200], $d.title)
    Start-Sleep -Milliseconds 400
}

Write-Output ""
$qa = '("diabetic kidney disease"[tiab] OR "diabetic nephropathy"[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR "LASSO"[tiab] OR "XGBoost"[tiab]) AND 2024:2026[dp]'
$qb = '(diabetes[tiab] OR diabetic[tiab]) AND (bioinformatics[tiab] OR "gene expression omnibus"[tiab] OR GEO[tiab]) AND ("machine learning"[tiab] OR "random forest"[tiab] OR LASSO[tiab] OR "deep learning"[tiab]) AND 2024:2026[dp]'
foreach ($q in @($qa, $qb)) {
    $u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=0&term=" + [uri]::EscapeDataString($q)
    $r = Invoke-RestMethod $u -TimeoutSec 40
    Write-Output ("recount => " + $r.esearchresult.count + " | " + $q.Substring(0, 60))
    Start-Sleep -Milliseconds 400
}
