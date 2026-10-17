# Code Index and Run Order

**Project**: IL-6 / hepcidin / bilirubin triangulation (clinical database x transcriptome x Mendelian randomization)
**Packaging date**: 2026-10-09
**Script count**: 53 (Python **47** + R **6**), divided into 5 groups by analysis layer

> This is the English translation of `03_Code/00_代码清单与运行顺序.md` (original Chinese version preserved unchanged). Local working directories were renamed to ASCII equivalents for public release; no analysis logic or numeric value was altered.

> All scripts are the actual production versions used in this project, not example code. Paths, database connections, and API endpoints are configurations verified to work on the local machine.

---

## Runtime Environment

| Component | Version / Location |
|---|---|
| Python | 3.13.12 (`C:\Users\44937\.workbuddy\binaries\python\envs\default\Scripts\python.exe`) |
| Python dependencies | pandas, numpy, scipy, statsmodels, psycopg2, Pillow, urllib (standard library) |
| R | 4.4.3 (`D:\R&Rstudio\R-4.4.3\bin\x64\Rscript.exe`) |
| R packages | **coloc 5.2.3**, **susieR 0.14.2**, ggplot2, dplyr, patchwork, scales, tidyr |
| Database | PostgreSQL `localhost:5432`, `user=postgres`, databases `mimiciii` / `mimiciv10` / `pic` |
| External APIs | OpenGWAS (requires environment variable `OGWAS_JWT`), Ensembl REST, GWAS Catalog v2, NCBI GEO, FinnGen R12 GCS |

**Two runtime constraints that must be noted** (pitfalls encountered in this project; see SKILL):

1. **R cannot handle Chinese (non-ASCII) paths**: `Rscript script.R "D:/workbuddy_workspace/..."` decodes the argument into garbled characters.
   -> Workaround: first copy the data to a pure-ASCII directory (e.g. `%TEMP%\coloc4`), call R with the ASCII path, then copy back after it finishes.
2. **R cannot directly read CSVs containing Chinese characters** (`invalid input found on input connection`, returns 0 rows).
   -> Workaround: perform all label mapping on the Python side, and export a **pure-ASCII, BOM-free** copy (named `*_ascii.csv` in this project) for R.

---

## Recommended Run Order

```
(1) finngen_tabix.py           (library imported by others; no need to run separately; place in the same directory or add to sys.path)
(2) 03_Mechanism_GEO           (no external dependencies, can be run first; large download volume, recommended to run first)
(3) 01_MR_Genetics             (requires OGWAS_JWT; outputs are reused by (4))
(4) 02_Colocalization          (depends on the regional data cached in (3))
(5) 04_Clinical_MIMIC          (requires local PostgreSQL)
(6) 05_Figures_Docs            (depends on all the result tables above)
```

---

## Group 1: `01_MR_Genetics` (16)

### Data Retrieval and Tools

| Script | Function |
|---|---|
| `finngen_tabix.py` | * **Pure-Python BGZF + tabix regional data fetcher** (does not depend on htslib/pysam). FinnGen summary statistics and `lab_values/` are each ~200 MB-10 GB and ship with a `.tbi`; use it to fetch only the target region (measured 3 seconds per region). Includes `read_tbi` / `reg2bins` / `bgzf_decompress` / `fetch_region` / `fetch_all`. **Column positions are derived automatically from the tbi** (`col_seq`/`col_beg` are 1-based; hard-coding them silently returns 0 rows) |
| `20261008_download_finngen_neonatal.py` | Download of FinnGen R12 neonatal phenotype summary statistics (range chunking + resumable download) |

### MR Primary Analysis

| Script | Function | Key Conclusion |
|---|---|---|
| `20261008_mr_crp_bilirubin.py` | CRP <-> total bilirubin **bidirectional** cis-MR | All negative, very small SE (0.0071-0.093) => credible negative |
| `20261008_mr_il6_saa_bilirubin.py` | **IL6R / SAA1 pQTL -> bilirubin** (primary analysis) | IL6R positive (+0.0192, p=3.1e-17); SAA1 negative |
| `20261008_mr_il6_saa_secondary.py` | MR for secondary exposures / outcomes | Direction-consistency check |
| `20261008_mr_ciseqtl_finngen.py` | cis-eQTL MR x FinnGen outcomes | Supplementary exposure |
| `20261009_mr_il6r_iron_erythroid_finngen.py` | * **IL6R pQTL -> FinnGen R12 quantitative iron/erythroid traits** (15 OMOPIDs, N=5,055-447,271) | Serum iron up, TSAT up, sTfR down, MCV up, MCH up, RDW down |

### Instrument Checks and Merging

| Script | Function |
|---|---|
| `20261009_check_il6r_instrument_independence.py` | Pairwise LD among the 4 IL6R instruments (FinnGen official LD matrix) + comparison with the FINEMAP credible set |
| `20261009_mr_il6r_IVW_correlated.py` | * **IVW with correlated instruments** (`β=(gᵀΣ⁻¹b)/(gᵀΣ⁻¹g)`, `Ψ_ij=ρ_ij·se_i·se_j`) + **Cochran Q** heterogeneity screen |
| `20261009_finemap_credible_set_il6r.py` | FinnGen FINEMAP credible sets (IL6R region, multiple traits), regional fetching via tbi |
| `20261009_summarize_iron_mr.py` | **Direction-consistency summary** for the iron/erythroid MR (cell-by-cell determination against the expected sign of the model) |

### FinnGen Neonatal Side

| Script | Function |
|---|---|
| `20261009_mr_il6r_finngen_neonatal_jaundice.py` | IL6R -> neonatal jaundice regional extraction (**scan in chromosome order and stop early once the target is reached**) |
| `20261009_mr_il6r_finngen_full_neonatal.py` | Whole-file MR for 5 neonatal phenotypes (jaundice / sepsis / RDS / BPD / preterm low birth weight) |
| `20261009_finngen_neonatal_power_all5.py` | * **Power comparison across the 5 phenotypes** (threshold model: `i=φ(z)/K`, to avoid an underpowered null being misread as negative) |
| `20261009_finngen_power_and_hamp_review.py` | Power conversion + **re-organization of the HAMP evidence stratified by stimulation type** |
| `20261009_gwascatalog_il6r_hamp_assoc.py` | GWAS Catalog v2 API: pull **all reported associations** at the IL6R / HAMP loci |

---

### WARNING: Exploratory Branch (Not Included in the Manuscript)

| Script | Description |
|---|---|
| `20261008_mr_ciseqtl_finngen.py` | An attempt at **two-sample cis-eQTL MR** on the preterm/low-birth-weight FinnGen data. **The manuscript does not cite this branch** ("eQTL" appears 0 times in the whole manuscript); its outputs `MR_ciseqtl_preterm_result.csv` and `MR_instruments.json` are written to the external data directory `D:\workbuddy_workspace\03-data\MR_finngen`, **not included in the package**. The script is kept only to document this abandoned path. |

---

## Group 2: `02_Colocalization` (10)

> WARNING: This group contains **both new and old paths**. The old paths (`_run2` / `_coloc_susie.R` / `_prep`) have been falsified,
> but are **kept** to document "why they failed at the time". **The final version uses only `run_abf2.R` + `20261009_coloc_susie_fetch_ogwasld.py`**.

### OK Final Path

| Script | Function |
|---|---|
| `20261009_coloc_susie_fetch_ogwasld.py` | Fetch regional data for the exposures (4 IL6R pQTLs) and 11 outcomes from OpenGWAS + **obtain the European LD matrix via `POST /api/ld/matrix`** |
| **`run_abf2.R`** | ** **Final-version coloc.abf**: (1) explicit **allele harmonization** (including flips) (2) **global MAF donation** (3) **dual windows** (whole region vs exposure lead +/-20 kb) (4) saturation check. **All colocalization values in the figures/main text come from this script** |

### Instrument Checks

| Script | Function |
|---|---|
| `20261009_build_ld_matrix.py` | Build a correlation matrix from the FinnGen LD long table (including **nearest positive-definite projection**: eigenvalue truncation + diagonal renormalization) |
| `20261009_effect_ratio_constancy.py` | Effect-ratio constancy test (appears LD-free, but **in practice still requires the same LD source on both sides**) |

### WARNING: Old Paths (Falsified, Kept for Reference)

| Script | Reason It Failed at the Time |
|---|---|
| `20261009_coloc_mch_mcv_jaundice_prep.py` | Merged the exposure (GRCh37) and outcome (GRCh38) by **liftover + by position** => data mismatch |
| `20261009_coloc_run2.R` | Ran abf on the **whole region** => spanned two independent signals, forced an H3 call |
| `20261009_coloc_susie.R` | Fit a **European** pQTL using the **Finnish panel** => `susie_rss` did not converge |
| `20261009_coloc_susie_ogwasld.R` | European LD converged, but the **exposure used a saturated dataset** => arbitrary lead, degenerate fit |
| `20261008_coloc_fetch_regions.py`, `20261008_coloc_il6r_bilirubin.R` | First-round data fetching and coloc (also affected by the saturation and window issues) |

---

## Group 3: `03_Mechanism_GEO` (14)

### Search and Preprocessing

| Script | Function |
|---|---|
| `20261009_geo_mechanism_search.py` | GEO E-utilities **multi-angle search** (not betting everything on a single keyword) |
| `20261009_geo_probe_files.py` | List the file names and sizes in `matrix/` and `suppl/` (obtain Content-Length via HEAD) |
| `20261009_geo_search_fetal_neonatal.py` | Search for fetal/neonatal hepatocytes +/- IL-6 (6 keyword sets, 183 datasets => **no eligible data**) |

### Download and Verification

| Script | Function |
|---|---|
| `20261009_geo_download.py` / `_download2.py` / `_download3.py` | GEO FTP download |
| `20261009_geo_download4_resume.py` | * **Range chunking + per-chunk timeout retry + resumable download** (GEO measured only 3-10 KB/s; `urlopen`'s timeout is ineffective for slow streaming) |
| `20261009_geo_download_gpl6480.py` | Platform annotation (GPL6480) download |
| `20261009_geo_verify_downloads.py` | * Integrity check: compare Content-Length + chunk-by-chunk gzip decompression |

### Analysis

| Script | Data / System | Key Result |
|---|---|---|
| `20261009_geo_mechanism_analysis.py` | **GSE202045** (human primary hepatocytes + IL-6, 3 donors) | **UGT1A1 -1.30 (3/3)**, HAMP +5.44/+3.04 |
| `20261009_geo_circadian_analysis.py` | **GSE218242** (circadian time course) / **GSE218243** (BMAL1/ARNTL knockdown) | UGT1A1 -0.74 (4/4 time points); UGT1A1 **no 24 h rhythm** (cosinor BH q=0.45) |
| `20261009_geo_replication_analysis.py` | **GSE230325** (dHepaRG + LPS) / **GSE275342** (Hep3B/Huh7 + IL-6) | UGT1A1 -1.54/-1.23 (9/9); Hep3B **+0.06 (negative replication)** |
| `20261009_geo_GSE102006_analysis.py` | GSE102006 (HepG2 + cytokine cocktail) | **Platform GPL6480 lacks a UGT1A1 probe => dataset rejected** |
| `20261009_mechanism_summary.py` | **Direction-consistency summary** across all datasets | 5 of 6 systems were negative (primary -1.30, dHepaRG -1.39 (mean of 6 h/24 h), BMAL1-KD -0.74, HepG2 -0.23, Huh7 -0.09), while **dedifferentiated Hep3B was +0.06**; the manuscript uniformly states this as "all systems that retain the hepatic metabolic program were consistently downregulated" |

**General analysis conventions (implemented in every script)**: (1) **hard gate on model validity** (acute-phase genes must be consistently upregulated, otherwise the whole batch is discarded);
(2) **genome-wide log2FC shift + target-gene |FC| percentile** as a specificity test (superior to a p-value with n=8);
(3) **housekeeping genes + liver-identity genes** as a negative reference.

---

## Group 4: `04_Clinical_MIMIC` (11)

| Script | Function | Key Result |
|---|---|---|
| `20261008_neo_infection_bilirubin.py` | MIMIC-IV infection x bilirubin primary analysis | Culture-positive -> bilirubin -0.355 (p=6e-59) |
| `20261008_neo_infection_bilirubin_align.py` | Age-in-days alignment + **five-item confounder test** | All five survived |
| `20261009_mimic_nicu_iron_extract.py` | MIMIC-III NICU data extraction (**first check lab availability**) | Direct iron-metabolism measures were nearly zero; red-cell indices were plentiful |
| `20261009_mimic_nicu_iron_analysis.py` | Red-cell index design A (cross-sectional) / B (within-subject paired) | Non-CoNS positive -> MCV -4.82, MCH -1.65 |
| `20261009_mimic_nicu_iron_transfusion_robust.py` | * **Triple robustness for transfusion** + dose-response | All three versions held: adjusted / un-transfused only / **pre-transfusion labs only** |
| `20261009_mimic_nicu_bilirubin_phototherapy_robust.py` | * **Phototherapy confounding** control | Phototherapy rate 61.5% vs 60.5% (1.0x) => phototherapy is not a confounder in MIMIC-III |
| **`20261010_cohort_flow.py`** (new) | * **Cohort-flow recomputation and audit trail**: connect to the MIMIC-IV database to recompute 17,244 -> 16,623 -> two groups 3,562 -> models A/B/C (8,968 / 8,968 / **7,803**, β(log) -0.3549, SE 0.0219, p 6.06e-59); also reports the structure of the fully adjusted sample (**7,803 measurements from 3,375 admissions in 3,354 infants, 2.31 measurements per admission**) and re-estimates model C with **standard errors clustered by admission and by infant** (-0.3549, 95% CI -0.4282 to -0.2817, p 2.19e-21; -0.4284 to -0.2815, p 2.78e-21), plus the sex of the two compared groups (51.5% male); and recompute the full MIMIC-III pipeline from the already packaged CSVs | Outputs `cohort_flow.csv`; the 17,244 and 7,803 that previously existed only in the exploration records can now be recomputed from within the package |
| **`20261010_make_tables.py`** (new) | * **Table 1 / Table 2 generation and cell-by-cell verification**: rebuild the two main-text tables from `D_bilirubin_phototherapy_robust.csv`, `D_bilirubin_pairs.csv`, `C_transfusion_robust_results.csv`, and compare them cell by cell against the tables in the manuscript docx | Outputs `table1_source_values.csv` / `table2_source_values.csv`; verification result **0 cells differ** |
| **`20261015_culture_to_index_interval.py`** (new) | * **Culture-to-measurement interval (the X quoted in the Discussion)**: reproduces the analysis sample exactly (n = 710; 107 any-positive, 42 definite non-CoNS), then anchors the index measurement on clinical time — the median chronological time of the days 3-14 window measurements obtained after the first definite-pathogen culture — and subtracts the culture collection time | Outputs `culture_to_index_interval_definite.csv`; **22 of 42** admissions had at least one post-culture measurement inside the window, median interval **2.5 days** (IQR 1.7-3.2, range 0.5-12.7). Robust to analytic choices: all definite episodes 2.3 d, mean instead of median 2.6 d, four reported indices only 2.5 d |
| **`20261016_aprdrg_severity_sensitivity.py`** (new) | * **Severity-of-illness sensitivity analysis**: refits the published MIMIC-IV model with the APR-DRG severity-of-illness and risk-of-mortality subclasses (read from the same `hosp.drgcodes` APR row that supplies the birth-weight band, so no extra extraction is needed), restricts to the less severely ill strata, refits across specifications that differ in measurement density, and reports the E-value of the primary estimate | Outputs `aprdrg_severity_sensitivity.csv`; adding severity moves the estimate only from -0.355 to **-0.312** (27% lower, p 2x10-14), and severity subclasses 1-2 still give **-0.191** (17% lower, p 0.019) on 32 exposed admissions; E-value **2.2**; reproduces the published sample exactly (7,803 / 3,375 / 3,354) |
| **`20261017_bilirubin_transfusion_robust.py`** (new) | * **Transfusion confounding for the within-subject bilirubin outcome** — previously handled only for the erythrocyte indices. Rebuilds the 728 culture-anchored pairs and verifies them row by row against `D_bilirubin_pairs.csv`, then layers transfusion in as adjustment (any transfusion during the admission; transfusion inside the −48 h/+72 h window) and as restriction (no transfusion in the window; never transfused). Estimates are OLS with HC3 robust standard errors, the same specification as the phototherapy script | Outputs `E_bilirubin_transfusion_robust.csv` (7 models) and `E_bilirubin_pairs_transfusion.csv`; the difference is attenuated from **−1.07 to −0.89 mg/dL** (95% CI −1.45 to −0.33, p 0.0018) but survives every layer: **−0.99** (−1.70 to −0.29, p 0.0057, 531 pairs, 32 positive) with no transfusion in the window and **−1.23** (−2.18 to −0.29, p 0.011, 389 pairs, 21 positive) among never-transfused admissions. Transfusion is more frequent after a positive culture (68% vs 44%; 51% vs 25% within the window). The rebuild reproduces the published pairs to 8.9e-16 and the published main estimate (−1.0654) exactly |

---

## Group 5: `05_Figures_Docs` (3)

| Script | Function |
|---|---|
| `20261009_figure_data_prep.py` | Figure data preparation (aggregate the result tables into one CSV per figure per panel; includes the dose-response recomputed by querying MIMIC again) |
| `20261009_manuscript_figures.R` | * **Figure 1-6 plotting**: 17 cm wide, each figure output as PDF (vector) + PNG 300 dpi + TIFF 600 dpi. **Pure-ASCII source code**, with label mapping done on the Python side |
| `20261009_md2docx_generic.py` | Markdown -> Word (generic version), **renders `[TBD: ...]` in bold red** to make gaps easy to locate; outputs the paragraph count / table count / number of TBD spots. **New switch `--numbered`**: centered Arabic page numbers in the footer (PAGE field) + continuous line numbers in the body (`w:lnNumType`, inserted per the OOXML schema after `pgMar` and before `cols`). The manuscript was generated with this switch: <br>`md2docx_generic.py <manuscript.md> <manuscript.docx> --numbered` |

---

## Reproduction Notes (Important)

1. **OpenGWAS needs a JWT**: set the environment variable `OGWAS_JWT`. **Do not type the token by hand** -- it must be passed verbatim, and validated for
   `len(token)==546` and a signature segment that decodes to **256 bytes** (RS256/RSA-2048). A 401 has two causes (expired / truncated),
   so decode the payload's `exp` first before deciding.
2. **FinnGen R12 requires no authentication**: `https://storage.googleapis.com/finngen-public-data-r12/`, both `lab_values/` and
   `ld_matrix/` come with `.tbi`; **fetch only the region, do not download the whole package** (medical traits ~200 MB, LD matrix ~10 GB,
   `FINEMAP.snp.bgz` up to 410 MB).
3. **GEO downloads are slow**: must use range chunking + per-chunk timeout retry (see `_download4_resume.py`),
   `urlopen(timeout=)` is **ineffective** for slow streaming.
4. **Check platform probe coverage before downloading large files**: for example, GPL6480 has no UGT1A1 probe, so the entire dataset cannot test the primary endpoint.
5. **MIMIC requires local PostgreSQL**: `mimiciii` (mainly used for NICU) and `mimiciv10` (complete v1.0).
   The MIMIC-IV ICU module **does not contain neonates by design**, so neonatal vital signs / gestational age must be taken from MIMIC-III.
6. **WARNING: old and new versions of result files (found during the 2026-10-10 check; pay close attention)**:
   `04_Results/04_Clinical_MIMIC/B_within_subject_results.csv` is an **earlier version**; its `bili_total` row is
   **n = 726, adj β = -1.0682, 95% CI (-1.663, -0.473), p = 0.000433**.
   The manuscript, Table 1, and Figure 2 use the **later version** (`D_bilirubin_phototherapy_robust.csv` and `fig2A_bilirubin_mimic3.csv`):
   **n = 728, β = -1.0654, 95% CI (-1.660, -0.471), p = 0.000447**.
   -> The two estimates are nearly identical, but the **sample sizes differ**; when citing, **use only the 728 version**; this file serves only as an analysis audit trail.
   The same file also contains the 72 h before/after paired results for `mcv` / `mch` / `rdw` (n = 1166, ncase = 125),
   these are **exploratory outputs not used for any manuscript conclusion** (MCV/MCH/RDW responses are on the scale of weeks, so methodologically the 72 h design does not apply).

---

## Data Traceability (Manuscript Numbers -> Source Files)

**Full check on 2026-10-10**: across the whole manuscript (Abstract / Introduction / Results / Discussion / Methods / figure legends / tables),
all **188 distinct numeric values** were traced back to `04_Results/` one by one, and **187 could be located directly in a source file** (the only exception is the PhysioNet credential number),
the table below lists the key correspondences (for use by reviewers or reproducers):

| Number in Manuscript | Source File |
|---|---|
| MIMIC-IV: 17,244 -> 16,623 -> **3,562 (two groups: 213 culture-positive + 3,349 no infection exposure) -> 7,803 bilirubin measurements** | See item 1 of "Points of Caution" below |
| MIMIC-III: 7,940 / 6,469 / 8,255 episodes / 229 positive / 158 positive admissions / 75 definite pathogens / 710 / 726->728 pairs / 287 no-phototherapy subset | `nicu_base.csv` (7,992 rows; filtering age <= 28 gives 7,940), `nicu_culture.csv` (full set of 8,306 rows; 8,255 episodes after restricting to the 7,940 cohort), `cohort_flow.csv` (step-by-step recomputation), `D_bilirubin_pairs.csv` (728 rows), `fig2A_bilirubin_mimic3.csv` |
| Bilirubin -1.07 (95% CI -1.66~-0.47, p=4e-4), phototherapy-adjusted -1.05, no-phototherapy-only -1.01 | `D_bilirubin_phototherapy_robust.csv` (rows 2, 3, 4) |
| Table 1 definite-pathogen row 728 (19), -2.11 (-3.45, -0.78), p=0.0019; CoNS-only row 728 (46), p=0.049 | `D_bilirubin_phototherapy_robust.csv` (rows (5)(6), **pair level**; see item 5 of "Points of Caution") |
| MCV -4.82 fL, MCH -1.65 pg (definite pathogen, n=710/42) | `A_crosssection_results.csv` (row `positive blood culture (non-CoNS)`) |
| Dose-response (**unified to the primary-analysis sample n=710, groups 603 / 78 / 29**): Hct 42.78 / 39.98 / 37.33%; MCV 104.39/103.49/100.96; MCH 35.59/35.18/34.14; RDW 18.32/18.74/19.29 | `fig3B_iron_dose.csv` (see item 6 of "Points of Caution") |
| All three transfusion robustness versions held | `C_transfusion_robust_results.csv` (32 rows = 4 metrics x 2 exposures x 4 models) |
| Table 1 transfusion rows: −0.89 (−1.45, −0.33) p=0.0018 adjusting for any transfusion; −0.93 adjusting within the window; −0.99 (−1.70, −0.29) p=0.0057 with no transfusion in the window (531 pairs, 32 positive); −1.23 (−2.18, −0.29) p=0.011 never transfused (389 pairs, 21 positive) | `E_bilirubin_transfusion_robust.csv` |
| Genetics: MCH +0.0209, MCV +0.0183, RDW -0.0139, sTfR -0.0263, TSAT +0.0135, total bilirubin +0.0074; Cochran Q p 0.09-0.83; heterogeneity 0.013/0.036/0.0018/0.00097 | `fig4A_genetic_forest.csv` (contains Q and Q_p columns) |
| HAMP locus MCH p=1.7e-17 | `fig4B_hamp_locus.csv` |
| Colocalization narrow-window PP.H4: MCH 0.989, RDW 0.990, MCV 0.934, MCHC 0.973, PLT 0.988; RBC/RETIC/MRV 0.013/0.009/0.015 | `02_Colocalization/coloc_abf_final_summary.csv` (whole-region values also given) |
| Hepatocytes: UGT1A1 -1.30 / -1.54 / -1.23 / -0.74 / +0.06 | `05_Figure_Data/fig5_hepatocyte.csv` (the dHepaRG row is the 6 h/24 h mean -1.3873) |
| Neonatal GWAS power 5.8-11.1% | `finngen_neonatal_power_all5.csv` |
| LD r² 0.00-0.437 (between instruments) | `ld_il6r_tools.csv` (10 of 12 ordered pairs); the missing pair is described in item 4 of "Points of Caution" below |
| r² between leads 0.000-0.058 | `_ld_leads_vs_outcomes.json`, `_coloc_ld_verdict.json` |

**6 Points of Caution (that reviewers may ask about)**

1. **Cohort flow: the MIMIC-IV side is now fully closed (added 2026-10-10)**.
   A new `20261010_cohort_flow.py` -> `04_Results/04_Clinical_MIMIC/cohort_flow.csv`, recomputed against the database and written to disk:
   `17,244` (NICU/SCN) -> `16,623` (>=1 total bilirubin measurement; also equals the number of rows in `observ_infection_bilirubin_adm.csv`)
   -> two groups 3,562 -> n = `8,968` / `8,968` / **`7,803`** for models A/B/C (model C: β(log) -0.3549, SE 0.0219, p 6.06e-59).
   The MIMIC-III side `7,992 / 7,940 / 6,469 / 8,255 / 229 / 158 / 75` was also recomputed by the same script from the packaged CSVs.
   WARNING: **no packaged result file carries the following 2 numbers** (re-verified 2026-10-12: both can be recomputed from the public GEO processed matrices) --
   `60,278` (= number of genes with any sample > 0 in `GSE218243_l2tpm`, out of 114,671 features) and
   `10,821` (= **geometric mean** of CYP3A4 across the three dHepaRG control samples in `GSE230325_normcounts`, 10,820.81) --
   these are mechanistic descriptive denominators and were not involved in any test.
2. **The pQTL dataset file names use the old years**: `exp_IL6R_folkersen2018.csv`, `exp_IL6R_pietzner2020.csv`,
   `exp_IL6R_sun2020.csv` correspond to the papers' official years **Folkersen 2017 / Pietzner 2021 / Sun 2018**
   (the `sun*` file is actually INTERVAL, N=3,301, **not** deCODE). The file names are historical legacy; they were **not renamed** to avoid breaking script references;
   the manuscript and references always use the official years.
3. **Generation scripts for Table 1 / Table 2 have since been added (2026-10-10)**: `20261010_make_tables.py` rebuilds
   the two tables from `D_bilirubin_phototherapy_robust.csv`, `D_bilirubin_pairs.csv`, `C_transfusion_robust_results.csv`,
   `E_bilirubin_transfusion_robust.csv` and
   **compares them cell by cell with the manuscript docx, with 0 cells differing**.
   It also incidentally found and corrected one **double-rounding** error: for the Hct (Any positive, Main) p-value in Table 2, the source data are 0.02648,
   which the manuscript had written as 0.027; the correct value is **0.026**.
   Note also: the "never phototherapy (ICD-9)" subset of Table 1 (125 pairs, 5 positive cases) could not be fit by the model, and the script reports `not estimable` as in the original table.
4. **`ld_il6r_tools.csv` contains only 10 of the 12 ordered pairs**: the pair **rs4601580 <-> rs6690230** is missing
   (FinnGen R12 panel r² < 0.01, i.e. the lower bound of the range **0.00-0.437** cited in the manuscript; the values are in the complete matrix of exploration record §15.5(1)).
   When `20261009_mr_il6r_IVW_correlated.py` builds the correlation matrix Σ, it starts from the identity matrix and fills in only the pairs present in the table,
   **so missing pairs are treated as r = 0** -- this is essentially equivalent to the measured r² < 0.01 (r < 0.1) and does not affect β or p.
   If a reproducer recomputes using only the files in this package, they will obtain a lower bound of 0.083-0.437, so this is stated explicitly here.
5. **Rows (5)/(6) of Table 1 have been recomputed at the pair level (4th verification round, 2026-10-10)**:
   `20261009_mimic_nicu_bilirubin_phototherapy_robust.py` originally used an **admission-level** exposure
   (`P2.hadm_id.isin(cc.hadm_id)`, "whether a definite pathogen occurred during that admission"), which is a different level from the **pair-level** definition of rows (1)-(4) in the same table,
   and 91 + 42 = 133 > 65 positive pairs, which is self-contradictory (the annotation on the left side of Figure 2a is also affected).
   Changed to read the pathogen from that pair's own culture (new columns `cu_is_cons` -> `pos_definite` / `pos_cons_only`), with an assertion added to guarantee their sum equals the number of positive pairs.
   **Result**: (5) 728 (**19**), β = **-2.11** (-3.45, -0.78), p = **0.0019**; (6) 728 (**46**), β = -0.55, p = **0.049**.
   Original values kept for the record: (5) 728 (91) -0.7855 (p = 0.00358); (6) 728 (42) -0.5493 (p = 0.0560).
   **The direction and conclusion are unchanged, with a stronger effect** (definite pathogen changed from -0.79 to -2.11 mg/dL). All other rows (main model / (1) / (2)) are numerically completely unchanged, as verified by reproduction.
6. **The dose-response (Fig 3b) denominator has been unified to the primary-analysis sample (5th verification round, 2026-10-10)**:
   The `fig3B()` function in `20261009_figure_data_prep.py` originally required only "birth weight + the metric itself",
   which made the denominators of the metrics inconsistent with each other and incomparable to the primary cross-sectional analysis (n = 710): MCV/MCH 712, RDW 710, **Hct 957** (120 positive != the 107 in the primary analysis).
   Changed to `core = wide.dropna(subset=["mcv","mch","rdw"])`, the same 710-case sample as the primary analysis,
   with all four metrics grouped as **603 / 78 / 29** (107 positive in total, consistent with the primary analysis).
   **Numerical changes occurred only for Hct**: 42.88 / 39.65 / 37.18% -> **42.78 / 39.98 / 37.33%**;
   the means of MCV, MCH, and RDW are unchanged (group 0 changed from 605 to 603 cases; the mean difference is < 0.01).
   The manuscript was updated accordingly: the dose-response sentence in Results now gives the adjusted β for all four ">=2 vs 0" contrasts
   (MCV -3.64 fL, p = 0.052; MCH -1.40 pg, 95% CI -2.56~-0.25, p = 0.017; Hct -2.86%, p = 0.005; RDW +0.14, p = 0.74, all adjusted for birth weight),
   and the Fig 3 legend (b) now reads "the same 710 cases, 603 / 78 / 29". Original values kept for the record: MCV/MCH n = 605/78/29, RDW 603/78/29, Hct 837/85/35.

7. **The result tables and the main-text Table 1 are now "row-by-row traceable" (7th verification round, 2026-10-10, including code review)**:
   `20261009_mimic_nicu_bilirubin_phototherapy_robust.py` originally had two traceability defects:
   - **(1) Models with insufficient sample were not written to the table**: `rep()` only printed and did not `out.append` when `len(d) < 30 or d.positive.sum() < 8`,
     which meant row 5 of the main-text Table 1, "(3) never phototherapy during that admission only (ICD), **n = 125 (5 positive)**, not estimable",
     was **not traceable** in `D_bilirubin_phototherapy_robust.csv`. Changed to also write the row (`estimable=0`, estimate left blank),
     and re-running confirmed **n = 125, 5 positive**, consistent with the main text.
   - **(2) The redundant duplicate row (4) was discarded**: `rep("(4) non-CoNS positive vs negative", P.assign(positive=lambda d: d.positive), ...)`
     in which `assign` is a no-op, making row (4) **numerically identical** to "+ adjusted for pre/post time interval and baseline value",
     and that row had already been superseded by the pair-level rows (5)/(6) and was not cited in the main text. It has been deleted.
   **Result**: the CSV is now **7 rows, corresponding one-to-one with the 7 rows of the main-text Table 1**, and a new `estimable` column was added.
   `04_Results/05_Figure_Data/table1_source_values.csv` was updated accordingly
   (previously the old values 728 (91) / -0.79 / 0.0036 and 728 (42) / -0.55 / 0.056,
   now corrected to 728 (19) / -2.11 / 0.0019 and 728 (46) / -0.55 / 0.049).
   After re-running, all other rows **fully reproduce** (main model -0.7921 p = 0.0179; +adjusted -1.0654 p = 0.000447;
   (1) -1.0528 p = 0.000535; (2) 287(25) -1.0101 p = 0.0108; (5) 728(19) -2.1137 p = 0.00193; (6) 728(46) -0.5489 p = 0.0494).

8. **Data-compliance note**: `nicu_cbc.csv` (5.4 MB) and `nicu_weight.csv` (180 MB, 2.56 million rows) are raw MIMIC-III extraction data
   subject to the PhysioNet data use agreement, and **must not be publicly released with the manuscript**, so they are not placed in the submission package `04_Results/`;
   the main-text "Data availability" statement already declares that they must be obtained via PhysioNet credentials, and that derived data tables are available from the corresponding author on request.
   What is kept in the submission package is the **result layer** and the **figure-data layer**, which is sufficient to verify all reported values.

