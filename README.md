# Analysis code: Culture-proven neonatal infection, total bilirubin and an erythrocyte signature

Code accompanying the manuscript:

> **Culture-proven neonatal infection: lower total bilirubin and an erythrocyte signature confined to definite pathogens**
> Rui Jin, Yaozhi Qi, Xingyu Zhang, Shouhong Yang, Yang Ren and Yan Gao

This repository holds 61 scripts (54 Python, 7 R). The clinical analyses (two neonatal cohorts,
MIMIC-IV v1.0 and MIMIC-III v1.4) and the hepatocyte transcriptome analyses (GEO) reproduce every
number, table and figure reported in the current manuscript, including the phototherapy mediation
decomposition added on 24 November 2026 and the MIMIC-IV within-subject replication added on
25 November 2026.

**Scope note (24 November 2026).** The manuscript accompanied by this repository now reports the two
clinical cohorts and the hepatocyte transcriptomes only. The Mendelian randomisation and
colocalisation analyses in `01_MR_Genetics` and `02_Colocalization` belonged to an earlier version
of the work; they are retained so that the analysis history stays auditable, but their results are
**not** reported in the current manuscript and should not be read as manuscript results.

## Start here

The script-by-script index, `00_Code_Index_and_Run_Order.md`, contains the recommended run order,
the software versions, the two runtime pitfalls described below, and a table that maps each number
in the manuscript to the script and result file that produces it.

## Repository layout

| Folder | Scripts | Scope |
|---|---|---|
| `01_MR_Genetics` | 17 Python | IL6R cis-pQTL and HAMP instruments; MR on neonatal jaundice, bilirubin and erythrocyte indices; substitute sensitivity analyses for the correlated four-SNP instrument set (single-instrument Wald ratios, leave-one-instrument-out, weighted median, MR-Egger, LD sensitivity; Table S7); IL6R fine-mapping and credible set; instrument independence checks; FinnGen R12 retrieval |
| `02_Colocalization` | 5 Python, 5 R | Regional data and LD matrices; `coloc.abf` with allele harmonisation and dual windows; SuSiE-based colocalisation |
| `03_Mechanism_GEO` | 14 Python | Search, download, verification and analysis of the human hepatocyte transcriptomic series, including circadian and replication analyses |
| `04_Clinical_MIMIC` | 16 Python | Cohort construction from the local MIMIC databases; bilirubin and phototherapy models; erythrocyte-index models with transfusion adjustment; transfusion adjustment for the within-subject bilirubin outcome; APR-DRG severity-of-illness sensitivity analysis; **phototherapy mediation decomposition (7.2% of the association; 1.6% change when phototherapy periods are excluded)**; **within-subject (48 h before / 72 h after) replication of the paired design in MIMIC-IV (1,582 pairs, 148 positive; -0.42 mg/dL, 95% CI -0.74 to -0.09)**; manuscript tables |
| `05_Figures_Docs` | 2 Python, 2 R, 7 CSV | Figure data preparation, all six manuscript figures, the plotting input data for Figures 2-5, and the Markdown-to-Word converter |

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
local PostgreSQL instance loaded with the MIMIC-III and MIMIC-IV schemas. The connection
parameters (host, port, user, password, database) are defined in each script; the password is a
placeholder for a local development instance and must be replaced with your own configuration
before running. No credential in this repository grants access to any MIMIC data.

## Environment

| Component | Version |
|---|---|
| Python | 3.13 |
| Python packages | numpy, pandas, scipy, statsmodels, psycopg2, python-docx, Pillow |
| R | 4.4.3 |
| R packages | coloc 5.2.3, susieR 0.14.2, ggplot2, dplyr, patchwork, scales, tidyr |
| Database | PostgreSQL 17, holding the MIMIC-III and MIMIC-IV schemas |

## Run order

1. `03_Mechanism_GEO` - no external dependencies; large downloads, so run first
2. `01_MR_Genetics` - requires `OGWAS_JWT`
3. `02_Colocalization` - reuses the regional caches written in step 2
4. `04_Clinical_MIMIC` - requires the local PostgreSQL databases
5. `05_Figures_Docs` - requires all result tables produced above

`01_MR_Genetics/finngen_tabix.py` is a shared helper imported by other scripts; it is not meant to
be executed on its own.

## Where the scripts expect their input

No primary data are redistributed in this repository. Every Python script declares its input
directory as `DIR` (usually in the first 40 lines) and, where relevant, an output directory
as `OUT`. These hold the author's local working paths, renamed to pure ASCII for portability
- for example `D:\il6_hepcidin_work\03-data\MR_finngen` - and are not required to exist on
another machine. Point `DIR`/`OUT` at your own copy of the harmonised tables (or place the
tables at those paths) before running. The input tables of the MR layer are themselves
produced by the earlier scripts in `01_MR_Genetics` from public sources (FinnGen R12,
OpenGWAS), so the run order above is not optional.
`01_MR_Genetics/20261112_mr_sensitivity_substitutes.py` additionally honours the `MR_DIR`
and `MR_OUT` environment variables, so the sensitivity analyses behind Table S7 can be
regenerated without editing the source.

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
version 1.7.6.

## Revision note (26 November 2026)

**Release 1.7.6 (Fig. 2a plotting data).** The plotting-input CSV for Fig. 2a (<code>05_Figures_Docs/fig2A_bilirubin_mimic3_ascii.csv</code>) held the pre-1.7.3 pathogen-class counts (19 definite / 46 CoNS; -2.11 / -0.55 mg/dL) while every analysis script and the manuscript already used the unified rule (20 / 45; -2.08 / -0.61 mg/dL). Corrected so the archived data reproduce the published Fig. 2a and Table 1; no estimate, table or conclusion changes.

**Release 1.7.5 (MIMIC-IV erythrocyte replication).** The erythrocyte signature reported in
MIMIC-III is now replicated cross-sectionally in the MIMIC-IV hospital module
(`20261126_mimiciv_erythrocyte_crosssectional.py`, Table S8, Fig. 3c). In 1,132 admissions with a
birth-weight band and an index on postnatal days 3–14, a definite pathogen was associated with lower
MCV (−2.08 fL, p = 0.004) and MCH (−0.63 pg, p = 0.004), whereas CoNS-only cultures showed no
difference (MCV +1.39 fL, p = 0.13; MCH +0.41 pg, p = 0.16) and RDW did not differ (p = 0.67). The
pattern reproduces at about half the MIMIC-III effect size. Transfusion could not be adjusted for and
no transfusion procedure code was identified in these admissions.

**Release 1.7.4 (culture-to-index interval).** The culture-to-measurement interval analysis
(`20261015_culture_to_index_interval.py`) still used the pre-1.7.3 organism rule. Recomputed under the
unified rule, the definite-pathogen group contains 45 admissions, of which 23 had at least one
post-culture erythrocyte-index measurement in the days 3–14 window and 22 had measurements only before
the culture; the median interval is 2.5 days (IQR 1.7–3.1, range 0.5–12.7). The sensitivity variants
give 2.3 days (all definite episodes as anchors), 2.8 days (mean rather than median post-culture time)
and 2.5 days (four reported indices only). The code index now reports 45 definite admissions and 23 of
45 with a post-culture measurement in place of the earlier 42 and 22 of 42.

**Release 1.7.3 (pathogen-class rule).** The MIMIC-III scripts used to classify a culture as
coagulase-negative when the organism string *contained* "staphylococcus". A culture growing a CoNS
together with another organism was therefore placed in the CoNS stratum, whereas the MIMIC-IV script
already required that every organism be a CoNS. All MIMIC-III scripts now apply that rule, which was
already the definition given in the manuscript:

| Analysis | 1.7.2 | 1.7.3 |
|---|---|---|
| MIMIC-III pairs, definite / CoNS-only | 19 / 46 | 20 / 45 |
| Within-subject fall, definite vs CoNS (mg/dL) | −2.16 vs −0.61 | −2.08 vs −0.61 |
| Within-subject contrast (p) | −1.54 (0.030) | −1.47 (0.032) |
| Erythrocyte cohort, definite / CoNS-only | 42 / 65 | 45 / 62 |
| MCV, definite (fL) | −4.82 | −4.20 |
| MCH, definite (pg) | −1.65 | −1.44 |
| Three-level contrast, MCV / MCH | −5.32 / −1.71 | −4.54 / −1.44 |

The conclusion is unchanged: the erythrocyte signature remains confined to definite pathogens and
the bilirubin fall remains graded only in MIMIC-III.

Two methodological corrections affect the three-level pathogen analyses and are the reason for
release 1.7.2. In both cohorts the two stratum coefficients and the contrast between them now come
from **one** model with culture-negative episodes as the shared reference, and the contrast is
taken from that model's coefficients and covariance matrix:

| Analysis | Before | After |
|---|---|---|
| MIMIC-III, definite vs CoNS | coefficients from one model, contrast from another | −2.16 vs −0.61, contrast −1.54 (p 0.030), all from the three-level model |
| MIMIC-IV, definite vs CoNS | contrast from a positive-pairs-only model | −0.42 vs −0.39, contrast −0.04 (−0.69 to 0.62, p 0.91), from the three-level model; the positive-pairs-only contrast (−0.10) is printed for audit only |

Subtracting two rounded coefficients from different models does not reproduce a model-based
contrast; that is what the earlier versions effectively invited.

## License

Released under the MIT License. See `LICENSE`.

## Contact

Rui Jin
Department of Neonatal Medical Center
Lianyungang Maternal and Child Health Hospital
669 Qindongmen Street, Lianyungang 222006, China
ORCID: 0009-0001-1803-3420
449378708@njmu.edu.cn
