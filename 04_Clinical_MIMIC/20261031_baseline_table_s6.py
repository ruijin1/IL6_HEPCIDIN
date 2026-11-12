# Build a baseline-characteristics table for the MIMIC-III neonatal cohort (n = 710),
# stratified by blood-culture result (no positive / CoNS only / definite pathogen).
# Reads only the cached MIMIC-III extracts; writes G_baseline_by_exposure.csv next to them.
# Usage: python 20261031_baseline_table_s6.py [csv_dir]
import os
import sys

import numpy as np
import pandas as pd

DIR = (sys.argv[1] if len(sys.argv) > 1
       else os.environ.get("MIMIC_CSV_DIR", r"D:\il6_hepcidin_work\03-data\MIMIC_NICU_iron"))

# Classification identical to the published analyses: a positive culture row is CoNS when the
# organism string contains "staphylococcus" (MIMIC-III reports "STAPHYLOCOCCUS, COAGULASE NEGATIVE");
# S. aureus ("STAPH AUREUS COAG +") therefore counts as a definite pathogen.
base = pd.read_csv(os.path.join(DIR, "nicu_base.csv"), encoding="utf-8-sig",
                   parse_dates=["dob", "admittime", "intime", "outtime"])
base = base[base.age_days_at_icu <= 28].copy()

w = pd.read_csv(os.path.join(DIR, "nicu_weight.csv"), encoding="utf-8-sig",
                usecols=["hadm_id", "label", "valuenum", "charttime"], parse_dates=["charttime"])
bw = w[w.label.astype(str).str.contains("Birth Weight", case=False, na=False)]
bw = bw[(bw.valuenum > 0.2) & (bw.valuenum < 7)]
bw = bw.sort_values("charttime").groupby("hadm_id").valuenum.first().rename("bw_kg")
base = base.merge(bw, left_on="hadm_id", right_index=True, how="left")

cul = pd.read_csv(os.path.join(DIR, "nicu_culture.csv"), encoding="utf-8-sig",
                  parse_dates=["charttime"])
cul = cul[cul.hadm_id.isin(base.hadm_id)]
pos = cul[cul.positive == 1].copy()
# Definition of the published analyses (20261009_mimic_nicu_iron_analysis.py):
# a definite pathogen is any organism row whose name does not contain "staphylococcus" or
# "coagulase negative"; CoNS-only admissions are positive admissions with no such row.
org = pos.orgs.astype(str).str.lower()
clean = pos[~org.str.contains("staphylococcus|coagulase negative", na=False)]
definite_hadm = set(clean.hadm_id)
cons_hadm = set(pos.hadm_id) - definite_hadm

idx = base.set_index("hadm_id").index
exp = pd.DataFrame(index=idx)
exp["exposure"] = np.where(idx.isin(definite_hadm), "definite",
                           np.where(idx.isin(cons_hadm), "cons", "none"))
pos_count = pos.groupby("hadm_id").orgs.nunique().rename("n_pos_cultures_org_rows")
episodes = cul.groupby("hadm_id").size().rename("n_culture_episodes")
npos = pos.groupby("hadm_id").charttime.nunique().rename("n_pos_culture_times")

lab = pd.read_csv(os.path.join(DIR, "nicu_cbc.csv"), encoding="utf-8-sig",
                  parse_dates=["charttime"])
lab = lab[lab.hadm_id.isin(base.hadm_id)].merge(base[["hadm_id", "dob"]], on="hadm_id", how="left")
lab["dol"] = (lab.charttime - lab.dob).dt.total_seconds() / 86400.0
lab = lab[(lab.dol >= 0) & (lab.dol <= 60)]
win = lab[(lab.dol >= 3) & (lab.dol <= 14) & (lab["var"].isin(["mcv", "mch", "rdw"]))]
agg = win.groupby(["hadm_id", "var"]).agg(val=("valuenum", "mean"), dol=("dol", "mean")).reset_index()
wide = agg.pivot(index="hadm_id", columns="var", values="val")
wide["dol_mean"] = agg.groupby("hadm_id").dol.mean()

d = base.set_index("hadm_id").join(exp).join(wide).join(episodes).join(npos)
d["exposure"] = d.exposure.fillna("none")
coh = d.dropna(subset=["mcv", "mch", "rdw", "bw_kg"]).copy()
print("cohort n =", len(coh), "| by exposure:", coh.exposure.value_counts().to_dict())


def med_iqr(s):
    s = s.dropna()
    if len(s) == 0:
        return "—"
    return "%.2f (%.2f–%.2f)" % (s.median(), s.quantile(.25), s.quantile(.75))

def n_pct(mask, denom):
    return "%d (%.0f%%)" % (mask.sum(), 100 * mask.sum() / denom)

rows = []
order = [("none", "No positive culture"), ("cons", "CoNS only"), ("definite", "Definite pathogen"),
         ("__any__", "Any positive culture")]


def row_for(g, label):
    n = len(g)
    return (label, n,
            med_iqr(g.bw_kg), int(g.bw_kg.isna().sum()),
            n_pct(g.gender == "M", n),
            med_iqr(g.age_days_at_icu),
            med_iqr(g.dol_mean),
            med_iqr(g.los),
            n_pct(g.hospital_expire_flag == 1, n),
            "%d / %d / %d" % ((g.n_pos_culture_times.fillna(0) == 0).sum(),
                              (g.n_pos_culture_times.fillna(0) == 1).sum(),
                              (g.n_pos_culture_times.fillna(0) >= 2).sum()))


for key, label in order:
    g = coh[coh.exposure == key] if key != "__any__" else coh[coh.exposure != "none"]
    rows.append(row_for(g, label))

out = pd.DataFrame(rows, columns=["Group", "n", "Birth weight, kg, median (IQR)",
                                  "Birth weight missing, n", "Male, n (%)",
                                  "Age at NICU admission, days, median (IQR)",
                                  "Postnatal age at index measurement, days, median (IQR)",
                                  "NICU length of stay, days, median (IQR)",
                                  "In-hospital death, n (%)",
                                  "0 / 1 / >=2 positive culture time points, n"])
out.to_csv(os.path.join(DIR, "G_baseline_by_exposure.csv"), index=False, encoding="utf-8-sig")

pd.set_option("display.width", 250)
print("\n=== baseline ===")
print(out.to_string(index=False))
print("\nwhole-cohort: n=%d, weight median %s, male %s, death %s" %
      (len(coh), med_iqr(coh.bw_kg), n_pct(coh.gender == "M", len(coh)),
       n_pct(coh.hospital_expire_flag == 1, len(coh))))
