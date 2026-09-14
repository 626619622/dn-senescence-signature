# Analysis code: senescence-related gene programmes in diabetic nephropathy

This repository contains the analysis code for the manuscript

> Cell-type-specific senescence-related gene programmes in diabetic nephropathy: a multi-cohort
> machine-learning signature validated against non-diabetic kidney disease

The study combines bulk transcriptomic discovery across multiple public cohorts, cross-cohort
replication, disease-specificity filtering, single-nucleus localisation in two independent
cohorts, machine-learning model building with calibration analysis, and cis-eQTL Mendelian
randomisation with colocalisation.

---

## Data sources

All data are public. Nothing in this repository redistributes raw data; the scripts download it.

| Dataset | Accession | Use in the study | Approximate size |
|---|---|---|---|
| Glomerular transcriptome, diabetic nephropathy | GSE96804 | Training / discovery | 14 MB (series matrix) |
| Glomerular transcriptome, diabetic kidney disease | GSE30528 | External validation | 2 MB |
| Tubulointerstitial transcriptome with disease controls | GSE99325 (GSE99340 subseries) | Specificity filter and validation | 25 MB |
| Renal biopsy cohort with multiple nephropathies | GSE104954 | Independent specificity validation | 672 MB (raw CEL) |
| Single-nucleus RNA, kidney cortex | GSE195460 | Single-nucleus discovery | 104 MB |
| Single-nucleus RNA, early diabetic nephropathy | GSE131882 | Single-nucleus replication | 1.4 GB |
| Blood cis-eQTL summary statistics | eQTLGen | Genetic instruments | 308 MB |
| Diabetic nephropathy GWAS | FinnGen R12 DM_NEPHROPATHY | Genetic outcome | 773 MB |

Access points:

- GEO: `https://ftp.ncbi.nlm.nih.gov/geo/series/`
- eQTLGen: `https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/`
- FinnGen R12: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/`

Three additional datasets were deliberately excluded and the reason is documented in the
Methods: GSE30122 (SuperSeries whose samples are identical to GSE30528 and GSE30529),
GSE47183 and GSE47184 (declared sample overlap with GSE99325).

---

## Requirements

- R 4.5.1 with Bioconductor 3.22
- R packages: GEOquery, limma, GSVA, msigdbr, Seurat, harmony, glmnet, randomForest, e1071,
  Boruta, xgboost, pROC, coloc, affy, data.table, AnnotationDbi, org.Hs.eg.db, and the platform
  annotation packages hta20transcriptcluster.db, hgu133a2.db, hgu133a.db, hgu133plus2.db,
  hgu133acdf, hgu133plus2cdf
- PowerShell 5.1 (for the data acquisition scripts)
- Internet access for data download

Exact package versions are listed in `environment/package-versions.txt`.

Installation of the required packages:

```powershell
.\environment\run-r.cmd 00-utility\r-install-deps.R
.\environment\run-r.cmd 00-utility\r-install-bioc.R
.\environment\run-r.cmd 00-utility\r-install-hta-db.R
```

### Windows note

On Windows systems with a non-ASCII user name, R cannot see the personal package library unless
the locale environment variables are cleared before R starts. `environment/run-r.cmd` does this
and should be used to launch every R script. Setting the locale from inside an R session does
not help.

PowerShell scripts in this repository are deliberately ASCII-only. On Windows PowerShell 5.1 a
BOM-less script containing non-ASCII characters is parsed with the local code page and fails.

---

## How to run

Scripts are numbered in execution order. Working directory should be the repository parent, and
all outputs are written to `codex_output/`.

### 1. Data acquisition and verification (`01-data-acquisition`)

```powershell
powershell -File 01-data-acquisition\step2-samples.ps1        # sample metadata for all cohorts
powershell -File 01-data-acquisition\step2-check-groups.ps1   # audit group labels
powershell -File 01-data-acquisition\step4-overlap.ps1        # GSM-level overlap between cohorts
powershell -File 01-data-acquisition\step3-download.ps1       # series matrices and small files
powershell -File 01-data-acquisition\step4-gse99340-dl.ps1    # GSE99325 expression matrices
powershell -File 01-data-acquisition\step3-verify.ps1         # file integrity check
```

`step4-overlap.ps1` is worth running before any analysis: it demonstrates that GSE30122 shares
all of its samples with GSE30528 and GSE30529, which is easy to miss and would invalidate a
validation design.

### 2. Preprocessing (`02-preprocessing`)

```powershell
.\environment\run-r.cmd 02-preprocessing\r-01-prepare.R       # annotation, gene collapse, QC
.\environment\run-r.cmd 02-preprocessing\r-16-gse104954.R     # RMA of 195 CEL files
```

### 3. Senescence analysis (`03-senescence-analysis`)

```powershell
.\environment\run-r.cmd 03-senescence-analysis\r-02-genesets.R
.\environment\run-r.cmd 03-senescence-analysis\r-03-deg.R
.\environment\run-r.cmd 03-senescence-analysis\r-04-senescence-score.R
.\environment\run-r.cmd 03-senescence-analysis\r-05-replication.R
.\environment\run-r.cmd 03-senescence-analysis\r-06-canonical-markers.R
```

### 4. Single-nucleus analysis (`04-single-nucleus`)

```powershell
.\environment\run-r.cmd 04-single-nucleus\r-07-sc-load.R
.\environment\run-r.cmd 04-single-nucleus\r-08-sc-cluster.R
.\environment\run-r.cmd 04-single-nucleus\r-09-sc-senescence.R
.\environment\run-r.cmd 04-single-nucleus\r-10-pseudobulk.R
.\environment\run-r.cmd 04-single-nucleus\r-18-sc-replication.R
.\environment\run-r.cmd 04-single-nucleus\r-18c-sc-scores.R
```

GSE131882 is deposited as double-gzipped dgecounts objects. They must be decompressed twice
before `readRDS` can open them; the decompression is done with the .NET GZipStream in
`r-18-sc-replication.R`.

### 5. Machine learning (`05-machine-learning`)

```powershell
.\environment\run-r.cmd 05-machine-learning\r-11-ml.R
.\environment\run-r.cmd 05-machine-learning\r-12-specificity.R
.\environment\run-r.cmd 05-machine-learning\r-13-specificity-filter.R
.\environment\run-r.cmd 05-machine-learning\r-14c-final-model.R
.\environment\run-r.cmd 05-machine-learning\r-15-recalibration.R
```

`r-11-ml.R` uses nested cross-validation so that feature selection is repeated inside every
training fold. Selecting features on the full dataset and then reporting cross-validated
performance inflates the result; the script would give an area under the curve of 1.000 either
way in this training cohort, but only the nested analysis is defensible.

### 6. Genetic analysis (`06-genetic-analysis`)

```powershell
.\environment\run-r.cmd 06-genetic-analysis\r-19-mr.R
.\environment\run-r.cmd 06-genetic-analysis\r-20-coloc.R
```

`r-20-coloc.R` needs no additional reference panel: `coloc.abf` requires only summary statistics
for shared SNPs, which both datasets provide. An LD reference (for example 1000 Genomes EUR)
would be needed only for `coloc.susie` or for the sensitivity analysis to the single-causal-
variant assumption.

### 7. Figures and tables (`07-figures-tables`)

```powershell
.\environment\run-r.cmd 07-figures-tables\r-21-figures.R
powershell -File 07-figures-tables\organise-supplementary.ps1
```

---

## Key design decisions

Three decisions materially affect the results and are documented here because they are easy to
get wrong.

**Within-cohort standardisation.** The cohorts were generated on different platforms and
GSE30528 had already been batch-corrected before submission. Genes are therefore z-scored
within each cohort before any cross-cohort model is applied. Pooling raw values across cohorts
would be invalid.

**Disease-specificity filtering.** A model trained against tumour-nephrectomy controls
separates diabetic from non-diabetic kidney with an area under the curve of 0.79 to 0.85, but
performs at chance against other biopsy-proven nephropathies. The 115 non-diabetic nephropathy
samples in GSE99325 are used to retain only genes that also differ between diabetic nephropathy
and other renal diseases. GSE104954, which is never used for selection, validates the result.

**Statistic unit in single-nucleus data.** Cell-level tests treat nuclei from the same donor as
independent observations and produce implausibly small p values. Every comparison is repeated
after aggregating to the sample level, and the sample-level result is the one reported.

---

## Outputs

All outputs are written to `codex_output/`:

| Folder | Contents |
|---|---|
| `geo/` | downloaded GEO data and per-sample metadata |
| `data/` | intermediate R objects (expression matrices, Seurat objects, fitted models) |
| `results/` | result tables |
| `figures/` | figures in PDF and PNG |
| `logs/` | full run logs |
| `mr/` | eQTLGen and FinnGen summary statistics |

---

## Licence

MIT Licence. See `LICENSE`.

## Contact

[Corresponding author name and email]
