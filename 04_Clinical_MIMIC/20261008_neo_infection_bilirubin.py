# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""MIMIC-IV v1.0: observational association between neonatal infection and bilirubin trajectory (the step that decides the direction of the whole study)

Design notes
  - Cohort: NICU/SCN admissions with a total bilirubin (50885) record
  - Exposure: (1) positive blood culture (2) culture-negative but >=5 days of systemic antimicrobials (clinical sepsis) (3) no infection
  - Outcome: peak bilirubin / first value / slope over the first 72 h / age at peak
  - * Key confounder: **late-onset infection occurs during the falling phase of bilirubin** (time confounding), so the main analysis is restricted to **early-onset (<3 days)** and adjusts for postnatal age
  - * Phototherapy is a "post-treatment variable" (affected by and affecting bilirubin) -> not used as an ordinary covariate, but handled as a stratified sensitivity analysis
"""
import sys, time, warnings
sys.stdout.reconfigure(encoding="utf-8"); warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import psycopg2

OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"
con = psycopg2.connect(host="localhost", port=5432, user="postgres",
                       password="1234", dbname="mimiciv10")
con.autocommit = True


def q(sql, label=None):
    t = time.time()
    df = pd.read_sql(sql, con)
    if label:
        print("   [%s] %d rows (%.1fs)" % (label, len(df), time.time() - t))
    return df


print("=" * 92)
print("MIMIC-IV v1.0: neonatal infection x bilirubin trajectory")
print("=" * 92)

# -- 1. Cohort ---------------------------------------------
adm = q("""
SELECT hadm_id, subject_id, birth_time, in_nicu, in_scn, los_days, hospital_expire_flag
FROM work.nicu_adm WHERE in_nicu OR in_scn
""", "NICU/SCN admissions")

# -- 2. Bilirubin (in-hospital, by postnatal age) ------------------------
bili = q("""
SELECT l.hadm_id, l.valuenum AS bili, l.age_hours
FROM work.neo_lab l
WHERE l.itemid = 50885 AND l.valuenum > 0 AND l.valuenum < 45
  AND l.age_hours IS NOT NULL
""", "bilirubin records")
bili["age_hours"] = pd.to_numeric(bili["age_hours"], errors="coerce")
bili["bili"] = pd.to_numeric(bili["bili"], errors="coerce")
bili = bili.dropna(subset=["age_hours", "bili"])
bili = bili[bili.age_hours >= 0]


def slope72(g):
    d = g[g.age_hours <= 72].sort_values("age_hours")
    if len(d) < 2:
        return np.nan
    return np.polyfit(d.age_hours.values, d.bili.values, 1)[0] * 24  # mg/dL/day


agg = bili.groupby("hadm_id").apply(lambda g: pd.Series({
    "bili_max": g.bili.max(),
    "bili_first": g.sort_values("age_hours").bili.iloc[0],
    "bili_first_age_h": g.sort_values("age_hours").age_hours.iloc[0],
    "bili_peak_age_h": g.loc[g.bili.idxmax(), "age_hours"],
    "n_bili": len(g),
    "slope72_per_day": slope72(g),
}), include_groups=False).reset_index()

# -- 3. Infection exposure -------------------------------
pos = q("""
SELECT DISTINCT hadm_id FROM work.neo_cx
WHERE spec_class='blood' AND culture_positive AND org_name IS NOT NULL
""", "culture-positive admissions")
sent = q("SELECT DISTINCT hadm_id FROM work.neo_cx WHERE spec_class='blood'", "admissions with a blood culture sent")

# Systemic antimicrobial course (truncated to the admission window), >=5 days
rx = q("""
SELECT hadm_id, starttime::timestamp AS s, stoptime::timestamp AS e, drug, route
FROM work.neo_rx WHERE drug_class='antimicrobial'
""", "antimicrobial orders")
SYS_ROUTES = ("IV", "IM", "PO", "PO/NG", "NG", "PO/OG", "OG", "SC", "SUBCUT",
              "IV DRIP", "IV BOLUS", "PIV", "IVPCA", "PB", "IV/PO", "ID")
rx = rx[rx.route.isin(SYS_ROUTES)].dropna(subset=["s", "e"])
rx = rx[rx.e > rx.s]
abx_days = {}
for h, g in rx.groupby("hadm_id"):
    iv = g[["s", "e"]].sort_values("s").values
    merged = []; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append((cs, ce)); cs, ce = s, e
    merged.append((cs, ce))
    hrs = sum((e - s) / np.timedelta64(1, "h") for s, e in merged)
    abx_days[h] = hrs / 24.0

df = adm.copy()
df["bc_pos"] = df.hadm_id.isin(set(pos.hadm_id))
df["bc_sent"] = df.hadm_id.isin(set(sent.hadm_id))
df["abx_days"] = df.hadm_id.map(abx_days)
df = df.merge(agg, on="hadm_id", how="inner")   # keep only admissions with bilirubin

# Infection groups
def grp(r):
    if r.bc_pos:
        return "1 Blood culture positive"
    if pd.notna(r.abx_days) and r.abx_days >= 5:
        return "2 Culture-negative/antimicrobials >=5d"
    if r.bc_sent or (pd.notna(r.abx_days) and r.abx_days > 0):
        return "3 Infection workup but below the thresholds above"
    return "4 No infection exposure"
df["infect_group"] = df.apply(grp, axis=1)

print()
print("=" * 92)
print("[1] Cohort and groups")
print("=" * 92)
print("   NICU/SCN admissions with bilirubin: %s" % format(len(df), ","))
print()
tab = df.groupby("infect_group").agg(
    n=("hadm_id", "count"),
    bili_max_median=("bili_max", "median"),
    bili_first_median=("bili_first", "median"),
    slope72_median=("slope72_per_day", "median"),
    peak_age_h=("bili_peak_age_h", "median"),
    first_age_h=("bili_first_age_h", "median"),
    n_measurements=("n_bili", "median"),
).round(2)
print(tab.to_string())

print()
print("=" * 92)
print("[2] * Main analysis: early-onset infection (age at first measurement <3 days) -- avoiding onset-time confounding")
print("=" * 92)
early = df[df.bili_first_age_h < 72].copy()
ee = early[early.bc_pos]
ne = early[~early.bc_sent & (early.abx_days.isna() | (early.abx_days == 0))]
print("   Early-onset culture-positive n=%d | no infection exposure n=%d" % (len(ee), len(ne)))
if len(ee) >= 5 and len(ne) >= 5:
    from scipy import stats
    for col, lab in [("bili_max", "peak bilirubin"), ("bili_first", "first bilirubin"),
                     ("slope72_per_day", "slope over first 72 h")]:
        a = ee[col].dropna(); b = ne[col].dropna()
        if len(a) < 5 or len(b) < 5:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        print("   %-14s positive %.2f (IQR %.2f-%.2f) vs no infection %.2f (IQR %.2f-%.2f)   p=%.4g"
              % (lab, a.median(), a.quantile(.25), a.quantile(.75),
                 b.median(), b.quantile(.25), b.quantile(.75), p))

print()
print("=" * 92)
print("[3] Full cohort: infection group vs no infection exposure (unadjusted, includes time confounding)")
print("=" * 92)
g1 = df[df.infect_group == "1 Blood culture positive"]
g4 = df[df.infect_group == "4 No infection exposure"]
from scipy import stats
for col, lab in [("bili_max", "peak bilirubin"), ("bili_first", "first bilirubin"),
                 ("slope72_per_day", "slope over first 72 h"), ("bili_first_age_h", "age at first measurement (h)")]:
    a = g1[col].dropna(); b = g4[col].dropna()
    if len(a) < 5 or len(b) < 5:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    print("   %-14s positive %.2f vs no infection %.2f   p=%.4g  (n=%d / %d)"
          % (lab, a.median(), b.median(), p, len(a), len(b)))

print()
print("=" * 92)
print("[4] Adjusted analysis: log(peak bilirubin) ~ infection + postnatal age + gestational-age band + weight band")
print("=" * 92)
# Gestational age (ICD; weeks taken from the in-database dictionary)
ga = q("""
SELECT d.subject_id, MIN(
  CASE d.icd_code
    WHEN '76521' THEN 23 WHEN '76522' THEN 24 WHEN '76523' THEN 25
    WHEN '76524' THEN 27 WHEN '76525' THEN 29 WHEN '76526' THEN 31
    WHEN 'P0721' THEN 22 WHEN 'P0722' THEN 23 WHEN 'P0723' THEN 24
    WHEN 'P0724' THEN 25 WHEN 'P0725' THEN 26 WHEN 'P0726' THEN 27
    WHEN 'P0731' THEN 28 WHEN 'P0732' THEN 29 WHEN 'P0733' THEN 30
    WHEN 'P0734' THEN 31 WHEN 'P0735' THEN 32 ELSE NULL END) AS ga_wk
FROM hosp.diagnoses_icd d
WHERE d.icd_code LIKE '7652%' OR d.icd_code LIKE 'P072%' OR d.icd_code LIKE 'P073%'
GROUP BY 1
""", "gestational-age codes")
bw = q("""
SELECT g.hadm_id, MAX(g.drg_code::text) AS apr_bw
FROM hosp.drgcodes g WHERE g.drg_type='APR'
  AND g.drg_code::text IN ('593','602','607','612','614','621','622','625','626','633','634','639','640')
GROUP BY 1
""", "APR weight band")
df2 = df.merge(ga, on="subject_id", how="left").merge(bw, on="hadm_id", how="left")
m = df2[(df2.bili_max > 0) & df2.bili_first_age_h.notna()].copy()
m["log_bili_max"] = np.log(m.bili_max)
m["infect"] = (m.infect_group == "1 Blood culture positive").astype(int)
m["ga_wk"] = pd.to_numeric(m.ga_wk, errors="coerce")
m["age_d"] = m.bili_first_age_h / 24.0
mm = m.dropna(subset=["log_bili_max", "age_d"])
import statsmodels.formula.api as smf
print("   Model A: infection term only")
rA = smf.ols("log_bili_max ~ infect", data=mm).fit(cov_type="HC3")
print("      infect coefficient %+.4f (SE %.4f, p=%.4g)  n=%d" % (rA.params.get("infect", np.nan),
      rA.bse.get("infect", np.nan), rA.pvalues.get("infect", np.nan), int(rA.nobs)))
print("   Model B: + age at first measurement")
rB = smf.ols("log_bili_max ~ infect + age_d", data=mm).fit(cov_type="HC3")
print("      infect coefficient %+.4f (SE %.4f, p=%.4g)" % (rB.params.get("infect", np.nan),
      rB.bse.get("infect", np.nan), rB.pvalues.get("infect", np.nan)))
mg = mm.dropna(subset=["ga_wk"])
if len(mg) > 200:
    rC = smf.ols("log_bili_max ~ infect + age_d + ga_wk", data=mg).fit(cov_type="HC3")
    print("   Model C: + gestational age (n=%d)" % int(rC.nobs))
    print("      infect coefficient %+.4f (SE %.4f, p=%.4g) | ga_wk %+.4f (p=%.3g)"
          % (rC.params.get("infect", np.nan), rC.bse.get("infect", np.nan),
             rC.pvalues.get("infect", np.nan), rC.params.get("ga_wk", np.nan),
             rC.pvalues.get("ga_wk", np.nan)))

df.to_csv(OUT + r"\observ_infection_bilirubin_adm.csv", index=False, encoding="utf-8-sig")
print()
print("Saved admission-level dataset:", OUT + r"\observ_infection_bilirubin_adm.csv")
con.close()
