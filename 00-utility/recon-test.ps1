$q1 = "diabetes"
$u1 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=0&term=" + [uri]::EscapeDataString($q1)
$r1 = Invoke-RestMethod $u1 -TimeoutSec 40
Write-Output ("simple 'diabetes' count => " + $r1.esearchresult.count)

$q2 = 'diabetes[tiab] AND Mendelian randomization[tiab]'
$u2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=0&term=" + [uri]::EscapeDataString($q2)
$r2 = Invoke-RestMethod $u2 -TimeoutSec 40
Write-Output ("MR query count => " + $r2.esearchresult.count + " warning=" + $r2.esearchresult.warning + " error=" + $r2.esearchresult.error)

$q3 = 'diabetes[tiab] AND "2024"[dp] : "2026"[dp]'
$u3 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&tool=codex&retmax=0&term=" + [uri]::EscapeDataString($q3)
$r3 = Invoke-RestMethod $u3 -TimeoutSec 40
Write-Output ("date query count => " + $r3.esearchresult.count + " warning=" + $r3.esearchresult.warning)
