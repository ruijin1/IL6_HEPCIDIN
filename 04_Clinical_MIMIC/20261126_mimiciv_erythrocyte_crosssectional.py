# -*- coding: utf-8 -*-
"""MIMIC-IV hospital module: cross-sectional erythrocyte indices by blood-culture result.

Adds the MIMIC-IV replication of the erythrocyte signature reported in MIMIC-III
(Table 2, Results section 4). The outcome is the mean MCV, MCH or RDW over postnatal
days 3-14; exposure is the three-level pathogen class used throughout the manuscript
(no positive culture / CoNS only / definite pathogen). Adjusted for APR-DRG birth-weight
band (the rank order of the 13 band codes used by the primary MIMIC-IV model) and mean
postnatal age at measurement, with HC3 heteroscedasticity-robust standard errors.

Transfusion cannot be adjusted for in MIMIC-IV: neonatal admissions contribute no ICU
stay, so the ICU chart-event tables carry no transfusion timing, and dated procedure
codes identify transfusion in almost none of them. The script counts how many admissions
in the analysis sample carry a transfusion procedure code, so that the limitation can be
stated quantitatively.

Usage
  python 20261126_mimiciv_erythrocyte_crosssectional.py [output_dir]
  Requires the local PostgreSQL database mimiciv10 with the work.* tables.
"""
import os
import sys

import numpy as np
import pandas as pd
import psycopg2
import statsmodels.formula.api as smf

OUT = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
os.makedirs(OUT, exist_ok=True)

con = psycopg2.connect(host="localhost", port=5432, user="postgres",
                       password="1234", dbname="mimiciv10")


def q(sql):
    return pd.read_sql(sql, con)


BANDS = ["593", "602", "607", "612", "614", "621", "622", "625", "626", "633", "634", "639", "640"]
RANK = {c: i + 1 for i, c in enumerate(BANDS)}

# ---- 1. erythrocyte indices, postnatal days 3-14, mean per admission ----
lab = q("""
SELECT hadm_id,
       avg(valuenum) FILTER (WHERE itemid=51250) AS mcv,
       avg(valuenum) FILTER (WHERE itemid=51248) AS mch,
       avg(valuenum) FILTER (WHERE itemid=51277) AS rdw,
       avg(age_hours)/24.0 AS age_d
FROM work.neo_lab
WHERE itemid IN (51250,51248,51277) AND age_hours BETWEEN 72 AND 336
GROUP BY hadm_id
""")
print("admissions with an index on days 3-14:", len(lab))

# ---- 2. birth-weight band (13 APR-DRG codes, same as the primary model) ----
bw = q("""
SELECT hadm_id, MAX(drg_code::text) AS apr_bw
FROM hosp.drgcodes
WHERE drg_type='APR' AND drg_code::text IN (%s)
GROUP BY hadm_id
""" % ",".join("'%s'" % c for c in BANDS))
bw["bw_rank"] = bw.apr_bw.map(RANK)
print("admissions with a birth-weight band:", len(bw))

# ---- 3. exposure, one rule for both cohorts ----
cx = q("""
WITH p AS (
  SELECT hadm_id, lower(trim(org_name)) AS org
  FROM work.neo_cx
  WHERE culture_positive AND org_name IS NOT NULL AND org_name <> ''
), a AS (
  SELECT hadm_id,
         bool_or(org NOT LIKE '%aureus%' AND org NOT LIKE '%coagulase negative%'
                 AND org NOT IN ('staphylococcus','staph')) AS is_def,
         bool_or(org LIKE '%aureus%' OR org LIKE '%coagulase negative%'
                 OR org IN ('staphylococcus','staph')) AS is_cons
  FROM p GROUP BY hadm_id
)
SELECT hadm_id,
       CASE WHEN is_def THEN 'definite'
            WHEN is_cons THEN 'cons_only'
            ELSE 'none' END AS exposure
FROM a
""")

# ---- 4. NICU/SCN restriction and plausible ranges ----
adm = q("SELECT hadm_id FROM work.nicu_adm WHERE in_nicu OR in_scn")
df = (lab.merge(bw[["hadm_id", "bw_rank"]], on="hadm_id", how="inner")
         .merge(cx, on="hadm_id", how="left")
         .merge(adm, on="hadm_id", how="inner"))
df["exposure"] = df.exposure.fillna("none")
df = df[(df.mcv.between(50, 150)) & (df.mch.between(10, 40)) & (df.rdw.between(10, 40))]
print("\nanalysis sample:", len(df))
print(df.exposure.value_counts().to_dict())

# ---- 5. transfusion code count (limitation, not adjustment) ----
tx = q("""
SELECT DISTINCT hadm_id FROM hosp.procedures_icd
WHERE (icd_version::int = 9  AND icd_code IN ('9904','9901'))
   OR (icd_version::int = 10 AND (icd_code LIKE '30233N%' OR icd_code = '30233P1'))
""")
df["transfused"] = df.hadm_id.isin(tx.hadm_id).astype(int)
print("admissions with a transfusion procedure code:",
      int(df.transfused.sum()), "of", len(df))

# ---- 6. models: three exposures x three indices, HC3 ----
df["definite"] = (df.exposure == "definite").astype(int)
df["anypos"] = (df.exposure != "none").astype(int)
df["cons_only"] = (df.exposure == "cons_only").astype(int)

rows = []
print("\n%-6s %-10s %6s %6s %8s %18s %10s" %
      ("index", "exposure", "n", "cases", "beta", "95% CI", "p"))
for outcome in ["mcv", "mch", "rdw"]:
    for expo in ["definite", "anypos", "cons_only"]:
        d = df.dropna(subset=[outcome, "bw_rank", "age_d"])
        m = smf.ols(f"{outcome} ~ {expo} + bw_rank + age_d", data=d).fit(cov_type="HC3")
        b = m.params[expo]
        lo, hi = m.conf_int().loc[expo]
        p = m.pvalues[expo]
        ncase = int(d[expo].sum())
        print("%-6s %-10s %6d %6d %+8.2f  (%+.2f, %+.2f)  %10.3g"
              % (outcome, expo, int(m.nobs), ncase, b, lo, hi, p))
        rows.append(dict(outcome=outcome, exposure=expo, n=int(m.nobs), ncase=ncase,
                         beta=b, lo=lo, hi=hi, p=p))

# ---- 7. three-level model, no positive culture as shared reference ----
print("\nthree-level model (reference: no positive culture)")
for outcome in ["mcv", "mch", "rdw"]:
    d = df.dropna(subset=[outcome, "bw_rank", "age_d"]).copy()
    d["grp"] = pd.Categorical(d.exposure, categories=["none", "cons_only", "definite"])
    m = smf.ols(f"{outcome} ~ C(grp, Treatment('none')) + bw_rank + age_d",
                data=d).fit(cov_type="HC3")
    for term in m.params.index:
        if "C(grp" in term:
            label = term.split("[T.")[1].rstrip("]")
            b = m.params[term]
            lo, hi = m.conf_int().loc[term]
            print("  %-4s %-10s %+7.2f (%+.2f, %+.2f) p=%.3g"
                  % (outcome, label, b, lo, hi, m.pvalues[term]))

out = pd.DataFrame(rows)
path = os.path.join(OUT, "H_mimiciv_erythrocyte_crosssectional.csv")
out.to_csv(path, index=False, encoding="utf-8-sig")
print("\nwritten:", path)
con.close()
