"""APR-DRG severity-of-illness sensitivity analysis for the primary MIMIC-IV model.

Fits the published model (log total bilirubin on culture-proven infection, adjusted for
age at measurement and APR-DRG birth-weight band, standard errors clustered by admission)
and then:
  (1) adds the APR-DRG severity-of-illness subclass and the risk-of-mortality subclass;
  (2) restricts to the less severely ill strata;
  (3) refits across specifications that differ in measurement density;
  (4) reports the E-value for the primary estimate.
Severity and mortality are read from the same hosp.drgcodes APR row that supplies the
birth-weight band, so no additional extraction is required.

Usage: python 20261016_aprdrg_severity_sensitivity.py [output_csv]
"""
import sys
import os
import csv
import warnings

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- configuration
DB = dict(host="localhost", port=5432, user="postgres", password="1234", dbname="mimiciv10")
OUT = sys.argv[1] if len(sys.argv) > 1 else "aprdrg_severity_sensitivity.csv"

APR_BW_CODES = ["593", "602", "607", "612", "614", "621", "622", "625", "626", "633", "634", "639", "640"]
BW_ORDER = {c: i + 1 for i, c in enumerate(APR_BW_CODES)}
SYSTEMIC_ROUTES = ("IV", "IM", "PO", "PO/NG", "NG", "PO/OG", "OG", "SC", "SUBCUT",
                   "IV DRIP", "IV BOLUS", "PIV", "IVPCA", "PB", "IV/PO", "ID")

import psycopg2
import statsmodels.formula.api as smf

con = psycopg2.connect(**DB)
con.autocommit = True
q = lambda s: pd.read_sql(s, con)

# ------------------------------------------------------- rebuild the published sample
adm = q("SELECT hadm_id, subject_id FROM work.nicu_adm WHERE in_nicu OR in_scn")
bili = q("""SELECT hadm_id, valuenum AS bili, age_hours FROM work.neo_lab
            WHERE itemid=50885 AND valuenum>0 AND valuenum<45 AND age_hours IS NOT NULL""")
bili["age_hours"] = pd.to_numeric(bili.age_hours, errors="coerce")
bili["bili"] = pd.to_numeric(bili.bili, errors="coerce")
bili = bili.dropna(subset=["age_hours", "bili"])
bili = bili[(bili.age_hours >= 0) & (bili.age_hours <= 336)]

pos = set(q("""SELECT DISTINCT hadm_id FROM work.neo_cx
               WHERE spec_class='blood' AND culture_positive AND org_name IS NOT NULL""").hadm_id)
sent = set(q("SELECT DISTINCT hadm_id FROM work.neo_cx WHERE spec_class='blood'").hadm_id)

rx = q("""SELECT hadm_id, starttime::timestamp s, stoptime::timestamp e, route
          FROM work.neo_rx WHERE drug_class='antimicrobial'""")
rx = rx[rx.route.isin(SYSTEMIC_ROUTES)].dropna(subset=["s", "e"])
rx = rx[rx.e > rx.s]
abx = {}
for h, g in rx.groupby("hadm_id"):
    iv = g[["s", "e"]].sort_values("s").values
    merged, (cs, ce) = [], tuple(iv[0])
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append((cs, ce)); cs, ce = s, e
    merged.append((cs, ce))
    abx[h] = sum((e - s) / np.timedelta64(1, "h") for s, e in merged) / 24.0

df = adm[adm.hadm_id.isin(set(bili.hadm_id))].copy()
df["bc_pos"] = df.hadm_id.isin(pos)
df["bc_sent"] = df.hadm_id.isin(sent)
df["abx_days"] = df.hadm_id.map(abx)

b = bili.merge(df[["hadm_id", "subject_id", "bc_pos", "bc_sent", "abx_days"]], on="hadm_id", how="inner")
b["group"] = np.where(b.bc_pos, "pos",
                      np.where(~b.bc_sent & (b.abx_days.isna() | (b.abx_days == 0)),
                               "no_infection", "other"))

bw = q("""SELECT hadm_id, MAX(drg_code::text) AS apr_bw FROM hosp.drgcodes
          WHERE drg_type='APR' AND drg_code::text IN (%s) GROUP BY 1"""
       % ",".join("'%s'" % c for c in APR_BW_CODES))
sev = q("""SELECT hadm_id, MAX(NULLIF(drg_severity,''))::int AS sev,
                  MAX(NULLIF(drg_mortality,''))::int AS mort
           FROM hosp.drgcodes WHERE drg_type='APR' GROUP BY 1""")

g = b[b.group.isin(["pos", "no_infection"])].copy()
m = g[g.bili > 0].merge(bw, on="hadm_id", how="left").merge(sev, on="hadm_id", how="left")
m["bw_rank"] = pd.to_numeric(m.apr_bw.map(BW_ORDER), errors="coerce")
m["logb"] = np.log(m.bili)
m["age_d"] = m.age_hours / 24.0
m["infect"] = (m.group == "pos").astype(int)
d = m.dropna(subset=["logb", "age_d", "bw_rank", "sev"]).copy()
d["n_meas"] = d.groupby("hadm_id").hadm_id.transform("size")
d["sev3"] = (d["sev"] >= 3).astype(int)
print("analysis sample: %d measurements / %d admissions / %d infants"
      % (len(d), d.hadm_id.nunique(), d.subject_id.nunique()))

BASE = "logb ~ infect + age_d + bw_rank"
rows = []


def fit(spec, data, family, label):
    r = smf.ols(spec, data=data).fit(cov_type="cluster", cov_kwds={"groups": data["hadm_id"]})
    ci = r.conf_int().loc["infect"]
    row = dict(family=family, model=label, measurements=int(r.nobs),
               admissions=int(data.hadm_id.nunique()),
               exposed_admissions=int(data[data.infect == 1].hadm_id.nunique()),
               beta=round(float(r.params["infect"]), 4),
               change_pct=round(100 * (np.exp(r.params["infect"]) - 1), 1),
               ci_low=round(float(ci[0]), 3), ci_high=round(float(ci[1]), 3),
               p=float(r.pvalues["infect"]))
    rows.append(row)
    print("   %-56s n=%-5d exposed=%-4d beta=%+.4f (%+.1f%%) p=%.3g"
          % (label, row["measurements"], row["exposed_admissions"], row["beta"],
             row["change_pct"], row["p"]))
    return r


print("\n[1] additional adjustment")
fit(BASE, d, "adjustment", "primary model (Table 2)")
fit(BASE + " + sev", d, "adjustment", "+ severity of illness (ordinal)")
fit(BASE + " + sev3", d, "adjustment", "+ severity of illness (binary, subclass >=3)")
fit(BASE + " + mort", d, "adjustment", "+ risk of mortality (ordinal)")
fit(BASE + " + sev + mort", d, "adjustment", "+ severity of illness and risk of mortality")

print("\n[2] restricted to stratum")
for lo, hi, label in [(1, 2, "severity subclass 1-2"), (1, 3, "severity subclass 1-3"),
                      (2, 2, "severity subclass 2"), (3, 3, "severity subclass 3")]:
    sub = d[(d.sev >= lo) & (d.sev <= hi)]
    if sub[sub.infect == 1].hadm_id.nunique() < 3:
        print("   %-56s too few exposed -> not fitted" % label); continue
    fit(BASE, sub, "stratum", label)
sub4 = d[d.sev == 4]
print("   severity subclass 4: %d unexposed admissions -> no internal comparison"
      % sub4[sub4.infect == 0].hadm_id.nunique())

print("\n[3] measurement density")
one = d.sort_values("age_hours").drop_duplicates("hadm_id")
one2 = d[(d.age_d >= 0) & (d.age_d <= 2)].sort_values("age_hours").drop_duplicates("hadm_id")
fit(BASE, d, "density", "all measurements (primary)")
fit(BASE, one, "density", "one measurement per admission")
fit(BASE + " + sev", one, "density", "one measurement per admission, + severity")
fit(BASE, one2, "density", "one measurement per admission, days 0-2")
fit(BASE + " + sev", one2, "density", "one measurement per admission, days 0-2, + severity")
fit(BASE, d[d.n_meas == 1], "density", "admissions contributing a single measurement")
fit(BASE, d[d.n_meas >= 2], "density", "admissions contributing >=2 measurements")

a = d.drop_duplicates("hadm_id")
print("   culture-positive share: single-measurement admissions %.1f%%, >=2 measurements %.1f%%"
      % (100 * a[a.n_meas == 1].infect.mean(), 100 * a[a.n_meas >= 2].infect.mean()))

print("\n[4] E-value for the primary estimate")
r = smf.ols(BASE, data=d).fit(cov_type="cluster", cov_kwds={"groups": d["hadm_id"]})
rr = float(np.exp(r.params["infect"])); lo = float(np.exp(r.conf_int().loc["infect"][0]))


def e_value(x):
    return x + np.sqrt(x * (x - 1)) if x >= 1 else 1 / x + np.sqrt((1 / x) * ((1 / x) - 1))


print("   risk ratio %.4f (CI lower bound %.4f) -> E-value %.2f (CI %.2f)"
      % (rr, lo, e_value(rr), e_value(lo)))

os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for row in rows:
        w.writerow(row)
print("\nwritten:", OUT, "(%d rows)" % len(rows))
con.close()
