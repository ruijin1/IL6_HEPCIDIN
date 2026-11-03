# -*- coding: utf-8 -*-
r"""MIMIC-III NICU: erythrocyte indices by pathogen class (CoNS-only vs definite pathogen vs none).

Added during revision of the manuscript, as a companion to the bilirubin negative-control analysis.
Cross-sectional design identical to 20261009_mimic_nicu_iron_analysis.py (Design A):
  mean of each index over postnatal days 3-14 per admission, adjusted for birth weight and mean postnatal age,
  OLS with HC3 robust standard errors.
Exposure is three-level: no positive culture (reference) / CoNS-only positive / definite pathogen positive.
The published definite-pathogen estimates (MCV -4.82, MCH -1.65) are reproduced first.
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

DIR = (sys.argv[1] if len(sys.argv) > 1 else
       os.environ.get("MIMIC_CSV_DIR", r"D:\workbuddy_workspace\03-data\MIMIC_NICU_iron"))
OUT = os.path.join(DIR, "F_erythrocyte_conscocci_negative_control.csv")

base = pd.read_csv(os.path.join(DIR, "nicu_base.csv"), encoding="utf-8-sig",
                   parse_dates=["dob", "admittime", "intime", "outtime"])
lab = pd.read_csv(os.path.join(DIR, "nicu_cbc.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
cul = pd.read_csv(os.path.join(DIR, "nicu_culture.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
w = pd.read_csv(os.path.join(DIR, "nicu_weight.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
b = w[w.label.astype(str).str.contains("Birth Weight", case=False, na=False)]
b = b[(b.valuenum > 0.2) & (b.valuenum < 7)]
bw = b.sort_values("charttime").groupby("hadm_id").valuenum.first().rename("bw_kg")

base = base[base.age_days_at_icu <= 28].copy().merge(bw, left_on="hadm_id", right_index=True, how="left")
cul = cul[cul.hadm_id.isin(base.hadm_id)]

pos = cul[cul.positive == 1]
cons_mask = pos.orgs.astype(str).str.lower().str.contains("staphylococcus|coagulase negative", na=False)
pos_clean = pos[~cons_mask]

idx = base.set_index("hadm_id").index
exp = pd.DataFrame(index=idx)
exp["pos_blood"] = idx.isin(pos.hadm_id).astype(int)
exp["pos_noncons"] = idx.isin(pos_clean.hadm_id).astype(int)
exp["pos_cons_any"] = idx.isin(pos[cons_mask].hadm_id).astype(int)
exp["cons_only"] = ((exp.pos_blood == 1) & (exp.pos_noncons == 0)).astype(int)

lab = lab[lab.hadm_id.isin(base.hadm_id)].merge(base[["hadm_id", "dob"]], on="hadm_id", how="left")
lab["dol"] = (lab.charttime - lab.dob).dt.total_seconds() / 86400.0
lab = lab[(lab.dol >= 0) & (lab.dol <= 60)]

win = lab[(lab.dol >= 3) & (lab.dol <= 14) & (lab["var"].isin(["mcv", "mch", "mchc", "rdw", "hb", "hct"]))]
agg = win.groupby(["hadm_id", "var"]).agg(val=("valuenum", "mean"), dol=("dol", "mean")).reset_index()
wide = agg.pivot(index="hadm_id", columns="var", values="val")
wide["dol_mean"] = agg.groupby("hadm_id").dol.mean()
wide = wide.join(exp).join(base.set_index("hadm_id")[["bw_kg"]])
wide = wide.dropna(subset=["mcv", "mch", "rdw", "bw_kg"])
print("cross-sectional sample = %d ; positive (any) = %d ; definite pathogen = %d ; CoNS-only = %d"
      % (len(wide), int(wide.pos_blood.sum()), int(wide.pos_noncons.sum()), int(wide.cons_only.sum())))

rows = []

def fit(tag, y, d, formula, term):
    d = d.dropna(subset=[y, term, "bw_kg", "dol_mean"])
    m = smf.ols(formula, data=d).fit(cov_type="HC3")
    ci = m.conf_int().loc[term]
    print("  %-46s %-4s n=%3d (cases %2d)  %+8.3f (%+7.3f, %+7.3f)  p=%.3g"
          % (tag, y, len(d), int(d[term].sum()), m.params[term], ci[0], ci[1], m.pvalues[term]))
    rows.append(dict(outcome=y, analysis=tag, term=term, n=len(d), ncase=int(d[term].sum()),
                     b=m.params[term], lo=ci[0], hi=ci[1], p=m.pvalues[term]))
    return m

print("\n-- reproduction of the published definite-pathogen estimates --")
for y in ["mcv", "mch", "rdw"]:
    fit("Definite pathogen (published specification)", y, wide, "%s ~ pos_noncons + bw_kg + dol_mean" % y, "pos_noncons")

print("\n-- CoNS-only stratum, reference = no positive culture --")
for y in ["mcv", "mch", "rdw"]:
    fit("CoNS-only positive culture", y, wide, "%s ~ cons_only + bw_kg + dol_mean" % y, "cons_only")

print("\n-- three-level model within one specification --")
wide["g_def"] = wide.pos_noncons
wide["g_cons"] = wide.cons_only
for y in ["mcv", "mch", "rdw"]:
    m = fit("Three-level: definite pathogen", y, wide, "%s ~ g_def + g_cons + bw_kg + dol_mean" % y, "g_def")
    fit("Three-level: CoNS-only", y, wide, "%s ~ g_def + g_cons + bw_kg + dol_mean" % y, "g_cons")
    dfi = m.params["g_def"] - m.params["g_cons"]
    se = np.sqrt(m.cov_params().loc["g_def", "g_def"] + m.cov_params().loc["g_cons", "g_cons"]
                 - 2 * m.cov_params().loc["g_def", "g_cons"])
    p = 2 * (1 - stats.norm.cdf(abs(dfi / se)))
    print("     contrast definite - CoNS (%s) = %+.3f (SE %.3f), p = %.3f" % (y.upper(), dfi, se, p))
    rows.append(dict(outcome=y, analysis="Contrast: definite minus CoNS-only", term="difference", n=len(wide),
                     ncase=np.nan, b=dfi, lo=dfi - 1.96 * se, hi=dfi + 1.96 * se, p=p))

print("\n-- unadjusted means --")
wide["grp"] = np.where(wide.pos_noncons == 1, "definite", np.where(wide.cons_only == 1, "cons", "none"))
for g, gg in wide.groupby("grp"):
    print("  %-9s n=%3d  MCV %6.2f  MCH %5.2f  RDW %5.2f" % (g, len(gg), gg.mcv.mean(), gg.mch.mean(), gg.rdw.mean()))

pd.DataFrame(rows).to_csv(OUT, index=False, encoding="utf-8-sig")
print("\nsaved:", OUT)
