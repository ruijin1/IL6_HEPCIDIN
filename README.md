# Analysis code: IL-6 signalling, hepcidin-mediated iron restriction and neonatal bilirubin metabolism

Code accompanying the manuscript:

> **IL-6 signalling, hepcidin-mediated iron restriction and neonatal bilirubin metabolism: a triangulation of clinical, Mendelian randomisation and transcriptomic evidence**
> Rui Jin, Xingyu Zhang, Shouhong Yang, Yang Ren and Yan Gao

This repository holds 52 scripts (46 Python, 6 R) that reproduce every number, table and figure
reported in the manuscript, across three evidence layers: two neonatal clinical cohorts
(MIMIC-IV v1.0 and MIMIC-III v1.4), genetic analyses (Mendelian randomisation and colocalisation
for IL6R and the HAMP locus), and human hepatocyte transcriptomes from GEO.

## Start here

The full script-by-script index - `00_Code_Index_and_Run_Order.md` in the repository, shipped as
`README_Code_Inventory.md` inside the code archive `Code_IL6-hepcidin-bilirubin-MR.zip` - contains
the recommended run order, the software versions, the two runtime pitfalls described below, and a
table that maps each number in the manuscript to the script and result file that produces it.

## Repository layout

| Folder | Scripts | Scope |
|---|---|---|
| `01_MR_Genetics` | 16 Python | IL6R cis-pQTL and HAMP instruments; MR on neonatal jaundice, bilirubin and erythrocyte indices; IL6R fine-mapping and credible set; instrument independence checks; FinnGen R12 retrieval |
| `02_Colocalization` | 5 Python, 5 R | Regional data and LD matrices; `coloc.abf` with allele harmonisation and dual windows; SuSiE-based colocalisation |
| `03_Mechanism_GEO` | 14 Python | Search, download, verification and analysis of the human hepatocyte transcriptomic series, including circadian and replication analyses |
| `04_Clinical_MIMIC` | 8 Python | Cohort construction from the local MIMIC databases; bilirubin and phototherapy models; erythrocyte-index models with transfusion adjustment; manuscript tables |
| `05_Figures_Docs` | 2 Python, 1 R | Figure data preparation, all six manuscript figures, and the Markdown-to-Word converter |

### Superseded and exploratory scripts

Two groups deliberately keep scripts whose output is **not** used in the manuscript, so that the
analysis history remains auditable. They are flagged as WARNING in the code index:

- `01_MR_Genetics/20261008_mr_ciseqtl_finngen.py` - an abandoned two-sample cis-eQTL branch. The
  manuscript reports no eQTL analysis, and this script's outputs are not part of the package.
- Several scripts in `02_Colocalization` (`20261009_coloc_run2.R`, `20261009_coloc_susie.R`,
  `20261009_coloc_mch_mcv_jaundice_prep.py`, `20261009_coloc_susie_ogwasld.R` and the two
  `20261008_coloc_*` first-round scripts) - earlier colocalisation paths that were superseded
  because of a genome-build mismatch, an over-wide window, a population mismatch between the pQTL
  and the LD panel, and a saturated exposure dataset.

**All colocalisation values reported in the manuscript come from `run_abf2.R` together with
`20261009_coloc_susie_fetch_ogwasld.py`.** Please do not treat the output of the scripts listed
above as a manuscript result.

## Data access

No patient-level or primary data are redistributed in this repository.

| Source | Access |
|---|---|
| MIMIC-III v1.4, MIMIC-IV v1.0 | Credentialed access on PhysioNet, including training in research with human subjects and acceptance of the data use agreement. Not redistributable |
| GEO (GSE202045, GSE230325, GSE275342, GSE218242, GSE218243, GSE102006) | Public |
| FinnGen Release 12 summary statistics | Public Google Cloud Storage bucket, no authentication required |
| OpenGWAS | Free account; the JWT token must be supplied in the environment variable `OGWAS_JWT` |
| GWAS Catalog, Ensembl REST | Public |

Because MIMIC data are access-controlled, the clinical scripts in `04_Clinical_MIMIC` require a
local PostgreSQL instance loaded with the MIMIC-III and MIMIC-IV schemas. The path to that
instance is configured at the top of each script.

## Environment

| Component | Version |
|---|---|
| Python | 3.13 |
| Python packages | numpy, pandas, scipy, statsmodels, psycopg2, python-docx, Pillow |
| R | 4.4.3 |
| R packages | coloc 5.2.3, susieR, ggplot2, dplyr, patchwork, scales, tidyr |
| Database | local PostgreSQL holding the MIMIC-III and MIMIC-IV schemas |

## Run order

1. `03_Mechanism_GEO` - no external dependencies; large downloads, so run first
2. `01_MR_Genetics` - requires `OGWAS_JWT`
3. `02_Colocalization` - reuses the regional caches written in step 2
4. `04_Clinical_MIMIC` - requires the local PostgreSQL databases
5. `05_Figures_Docs` - requires all result tables produced above

`01_MR_Genetics/finngen_tabix.py` is a shared helper imported by other scripts; it is not meant to
be executed on its own.

## Two runtime pitfalls

1. **R cannot handle non-ASCII paths.** Argument strings passed to `Rscript` from a directory whose
   path contains non-ASCII characters are decoded incorrectly. Route the data through a pure-ASCII
   working directory (for example under `%TEMP%`) and copy the results back afterwards.
2. **R cannot read CSVs that contain non-ASCII characters.** All label mapping is therefore done on
   the Python side, which exports a BOM-free, pure-ASCII copy of every figure-data table before R
   reads it.

## Citation

If you use this code, please cite the manuscript above and the archived release. The code is
archived on Zenodo under the permanent DOI `10.5281/zenodo.23344837`, which always resolves to the
latest release and is also listed in `CITATION.cff`; the release corresponding to the manuscript is
version 1.0.2.

## License

Released under the MIT License. See `LICENSE`.

## Contact

Rui Jin
Department of Neonatal Medical Center
Lianyungang Maternal and Child Health Hospital
669 Qindongmen Street, Lianyungang 222006, China
ORCID: 0009-0001-1803-3420
449378708@njmu.edu.cn
