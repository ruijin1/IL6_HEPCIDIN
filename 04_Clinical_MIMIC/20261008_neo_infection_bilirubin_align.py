# -*- coding: utf-8 -*-
"""Decisive validation: compare bilirubin between the infection and no-infection groups within the **same postnatal-age window**
  If the association is caused by "measurement at an earlier time", it should disappear once postnatal age is aligned;
  if a difference persists after alignment, this suggests a genuine biological difference.
  We also use the APR-DRG weight band (covering 98.2%) instead of ICD gestational age (only 14.2%) to enlarge the adjusted sample.
"""
import sys, time, warnings
sys.stdout.reconfigure(encoding="utf-8"); warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import psycopg2
from scipy import stats

con = psycopg2.connect(host="localhost", port=5432, user="postgres",
                       password="1234", dbname="mimiciv10")
con.autocommit = True


def q(sql):
    return pd.read_sql(sql, con)


adm = q("""SELECT hadm_id, subject_id, birth_time, in_nicu, in_scn FROM work.nicu_adm
           WHERE in_nicu OR in_scn""")
bili = q("""SELECT hadm_id, valuenum AS bili, age_hours FROM work.neo_lab
            WHERE itemid=50885 AND valuenum>0 AND valuenum<45 AND age_hours IS NOT NULL""")
bili["age_hours"] = pd.to_numeric(bili.age_hours, errors="coerce")
bili["bili"] = pd.to_numeric(bili.bili, errors="coerce")
bili = bili.dropna()
bili = bili[(bili.age_hours >= 0) & (bili.age_hours <= 336)]

pos = set(q("""SELECT DISTINCT hadm_id FROM work.neo_cx
               WHERE spec_class='blood' AND culture_positive AND org_name IS NOT NULL""").hadm_id)
sent = set(q("SELECT DISTINCT hadm_id FROM work.neo_cx WHERE spec_class='blood'").hadm_id)

rx = q("""SELECT hadm_id, starttime::timestamp s, stoptime::timestamp e, route
          FROM work.neo_rx WHERE drug_class='antimicrobial'""")
SYS = ("IV", "IM", "PO", "PO/NG", "NG", "PO/OG", "OG", "SC", "SUBCUT",
       "IV DRIP", "IV BOLUS", "PIV", "IVPCA", "PB", "IV/PO", "ID")
rx = rx[rx.route.isin(SYS)].dropna(subset=["s", "e"])
rx = rx[rx.e > rx.s]
abx = {}
for h, g in rx.groupby("hadm_id"):
    iv = g[["s", "e"]].sort_values("s").values
    mg = []; cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            mg.append((cs, ce)); cs, ce = s, e
    mg.append((cs, ce))
    abx[h] = sum((e - s) / np.timedelta64(1, "h") for s, e in mg) / 24.0

df = adm[adm.hadm_id.isin(set(bili.hadm_id))].copy()
df["bc_pos"] = df.hadm_id.isin(pos)
df["bc_sent"] = df.hadm_id.isin(sent)
df["abx_days"] = df.hadm_id.map(abx)

bili = bili.merge(df[["hadm_id", "bc_pos", "bc_sent", "abx_days"]], on="hadm_id", how="inner")
bili["group"] = np.where(bili.bc_pos, "Culture-positive",
                         np.where(~bili.bc_sent & (bili.abx_days.isna() | (bili.abx_days == 0)),
                                  "No infection exposure", "Other"))
bili = bili[bili.group.isin(["Culture-positive", "No infection exposure"])]

print("=" * 94)
print("* Bilirubin comparison after postnatal-age alignment (key point: if the difference disappears after alignment -> the original association was caused by measurement timing)")
print("=" * 94)
print("   window(h)   positive n   positive median(IQR)       no inf. n  no inf. median(IQR)        p")
print("   " + "-" * 88)
for lo, hi in [(0, 24), (24, 48), (48, 72), (72, 96), (96, 120), (120, 168)]:
    d = bili[(bili.age_hours >= lo) & (bili.age_hours < hi)]
    a = d[d.group == "Culture-positive"].bili
    b = d[d.group == "No infection exposure"].bili
    if len(a) < 10 or len(b) < 10:
        print("   %3d-%-5d   %6d   -                         %6d   -" % (lo, hi, len(a), len(b)))
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    print("   %3d-%-5d   %6d   %5.2f (%.2f-%.2f)   %6d   %5.2f (%.2f-%.2f)   p=%.4g%s"
          % (lo, hi, len(a), a.median(), a.quantile(.25), a.quantile(.75),
             len(b), b.median(), b.quantile(.25), b.quantile(.75), p,
             "  **" if p < 0.05 else ""))

# Average per infant per window, then do a paired window comparison
print()
print("=" * 94)
print("* Also restricting to infants with >=2 measurements (to remove differences in testing patterns)")
print("=" * 94)
cnt = bili.groupby("hadm_id").size()
multi = set(cnt[cnt >= 2].index)
d2 = bili[bili.hadm_id.isin(multi)]
print("   admissions with >=2 measurements: %d (positive %d / no infection %d)"
      % (d2.hadm_id.nunique(),
         d2[d2.group == "Culture-positive"].hadm_id.nunique(),
         d2[d2.group == "No infection exposure"].hadm_id.nunique()))
for lo, hi in [(24, 48), (48, 72), (72, 96)]:
    d = d2[(d2.age_hours >= lo) & (d2.age_hours < hi)]
    a = d[d.group == "Culture-positive"].bili
    b = d[d.group == "No infection exposure"].bili
    if len(a) < 8 or len(b) < 8:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    print("   %3d-%-5dh  positive %5.2f (n=%d) vs no infection %5.2f (n=%d)   p=%.4g%s"
          % (lo, hi, a.median(), len(a), b.median(), len(b), p, "  **" if p < 0.05 else ""))

# Adjust for the APR-DRG weight band
print()
print("=" * 94)
print("* Adjusted analysis: log(bilirubin) ~ infection + postnatal age + APR-DRG weight band (covering 98.2%, replacing ICD gestational age)")
print("=" * 94)
bw = q("""SELECT hadm_id, MAX(drg_code::text) AS apr_bw FROM hosp.drgcodes
          WHERE drg_type='APR' AND drg_code::text IN
          ('593','602','607','612','614','621','622','625','626','633','634','639','640')
          GROUP BY 1""")
# Approximate order of the DRG weight bands: 593=<750g ... 640=>2499g (rough ordering by code)
BW_ORDER = {"593": 1, "602": 2, "607": 3, "612": 4, "614": 5, "621": 6,
            "622": 7, "625": 8, "626": 9, "633": 10, "634": 11, "639": 12, "640": 13}
m = bili[bili.bili > 0].merge(bw, on="hadm_id", how="left")
m["bw_rank"] = pd.to_numeric(m.apr_bw.map(BW_ORDER), errors="coerce")
m["logb"] = np.log(m.bili)
m["age_d"] = m.age_hours / 24.0
m["infect"] = (m.group == "Culture-positive").astype(int)
import statsmodels.formula.api as smf
for f, lab in [("logb ~ infect", "A infection only"),
               ("logb ~ infect + age_d", "B + postnatal age"),
               ("logb ~ infect + age_d + bw_rank", "C + postnatal age + weight band")]:
    d = m.dropna(subset=["logb", "age_d"]) if "bw_rank" not in f else m.dropna(subset=["logb", "age_d", "bw_rank"])
    r = smf.ols(f, data=d).fit(cov_type="HC3")
    print("   %-22s n=%-7d infect %+.4f (SE %.4f, p=%.4g)" % (
        lab, int(r.nobs), r.params.get("infect", np.nan),
        r.bse.get("infect", np.nan), r.pvalues.get("infect", np.nan)))
con.close()
