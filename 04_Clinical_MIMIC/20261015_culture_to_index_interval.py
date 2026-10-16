# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MIMIC_NICU_iron). Analysis logic is unchanged.
r"""Culture -> erythrocyte-index measurement interval (MIMIC-III NICU, definite-pathogen group).

Purpose
-------
The Discussion argues that RDW is a slow-moving index and therefore cannot be expected to
respond within the days-3-14 observation window. That argument needs the elapsed time between
the positive blood culture and the erythrocyte index actually being measured. This script
computes that interval and reproduces the analysis sample it refers to.

Definitions (identical to those used for Table 2)
-------------------------------------------------
* Postnatal day      = (laboratory collection time - date of birth) expressed in days.
* Measurement window = postnatal days 3-14 (inclusive), using the same variable list as the
                       cross-sectional aggregation: MCV, MCH, MCHC, RDW and haematocrit.
* Definite pathogen  = positive blood culture in which none of the isolated organisms matched
                       "staphylococcus" or "coagulase negative".
* Index time         = median chronological time of those window measurements obtained AFTER
                       the first definite-pathogen culture of the admission.
* Interval           = index time - culture collection time (days).

Output: culture_to_index_interval_definite.csv (one row per definite-pathogen admission).
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd

DIR = r"D:\workbuddy_workspace\03-data\MIMIC_NICU_iron"
OUTDIR = DIR

# The variable list used verbatim by 20261009_mimic_nicu_iron_analysis.py. Note that "hb"
# does not match any label in the extracted table (haemoglobin arrives as hb_heme / hb_bg),
# so the window aggregates five indices: mcv, mch, mchc, rdw, hct.
WINVARS = ["mcv", "mch", "mchc", "rdw", "hb", "hct"]


def load():
    base = pd.read_csv(os.path.join(DIR, "nicu_base.csv"), encoding="utf-8-sig",
                       parse_dates=["dob", "admittime", "intime", "outtime"])
    lab = pd.read_csv(os.path.join(DIR, "nicu_cbc.csv"), encoding="utf-8-sig",
                      parse_dates=["charttime"])
    cul = pd.read_csv(os.path.join(DIR, "nicu_culture.csv"), encoding="utf-8-sig",
                      parse_dates=["charttime"])
    return base, lab, cul


def build_birthweight():
    w = pd.read_csv(os.path.join(DIR, "nicu_weight.csv"), encoding="utf-8-sig",
                    parse_dates=["charttime"])
    b = w[w.label.astype(str).str.contains("Birth Weight", case=False, na=False)]
    b = b[(b.valuenum > 0.2) & (b.valuenum < 7)]
    return b.sort_values("charttime").groupby("hadm_id").valuenum.first().rename("bw_kg")


def main():
    base, lab, cul = load()
    bw = build_birthweight()

    # ---- cohort and analysis sample, exactly as in the main analysis ----
    base = base[base.age_days_at_icu <= 28].copy()
    base = base.merge(bw, left_on="hadm_id", right_index=True, how="left")

    cul = cul[cul.hadm_id.isin(base.hadm_id)]
    pos = cul[cul.positive == 1]
    pos_clean = pos[~pos.orgs.astype(str).str.lower().str.contains(
        "staphylococcus|coagulase negative", na=False)]

    exp = pd.DataFrame(index=base.hadm_id.values)
    exp["pos_blood"] = base.set_index("hadm_id").index.isin(pos.hadm_id).astype(int)
    exp["pos_blood_noncons"] = base.set_index("hadm_id").index.isin(pos_clean.hadm_id).astype(int)

    lab = lab[lab.hadm_id.isin(base.hadm_id)].copy()
    lab = lab.merge(base[["hadm_id", "dob"]], on="hadm_id", how="left")
    lab["dol"] = (lab.charttime - lab.dob).dt.total_seconds() / 86400.0
    lab = lab[(lab.dol >= 0) & (lab.dol <= 60)]

    win = lab[(lab.dol >= 3) & (lab.dol <= 14) & (lab["var"].isin(WINVARS))]
    agg = (win.groupby(["hadm_id", "var"])
              .agg(val=("valuenum", "mean"), n=("valuenum", "size"), dol=("dol", "mean"))
              .reset_index())
    wide = agg.pivot(index="hadm_id", columns="var", values="val")
    wide["dol_mean"] = agg.groupby("hadm_id").dol.mean()
    wide = wide.join(exp).join(base.set_index("hadm_id")[["bw_kg"]])
    wide = wide.dropna(subset=["mcv", "mch", "rdw", "bw_kg"])

    print("analysis sample n = %d (any-positive %d, definite non-CoNS %d)"
          % (len(wide), int(wide.pos_blood.sum()), int(wide.pos_blood_noncons.sum())))
    assert len(wide) == 710 and int(wide.pos_blood_noncons.sum()) == 42, \
        "sample does not reproduce the published n=710 / 42"

    # ---- interval ----
    definite = wide.index[wide.pos_blood_noncons == 1].tolist()
    culture_time = (pos_clean[pos_clean.hadm_id.isin(definite)]
                    .sort_values(["hadm_id", "charttime"])
                    .groupby("hadm_id").charttime.first())
    winrows = lab[(lab.dol >= 3) & (lab.dol <= 14) & (lab["var"].isin(WINVARS))]

    rows = []
    for h in definite:
        wl = winrows[winrows.hadm_id == h]
        if h not in culture_time.index:
            rows.append(dict(hadm_id=h, n_window=len(wl), n_post=0, interval_days=np.nan))
            continue
        ct = culture_time.loc[h]
        post = wl[wl.charttime > ct]
        if len(post) == 0:
            rows.append(dict(hadm_id=h, n_window=len(wl), n_post=0, interval_days=np.nan))
            continue
        rows.append(dict(hadm_id=h, n_window=len(wl), n_post=len(post),
                         interval_days=(post.charttime.median() - ct).total_seconds() / 86400.0))

    R = pd.DataFrame(rows)
    have = R.dropna(subset=["interval_days"])
    q = have.interval_days
    print("definite admissions                  = %d" % len(definite))
    print("with >=1 post-culture window measure = %d" % len(have))
    print("without any post-culture measure     = %d" % (len(definite) - len(have)))
    print("median interval (days)               = %.3f" % q.median())
    print("IQR                                  = %.3f - %.3f"
          % (q.quantile(.25), q.quantile(.75)))
    print("range                                = %.3f - %.3f" % (q.min(), q.max()))

    out = os.path.join(OUTDIR, "culture_to_index_interval_definite.csv")
    R.to_csv(out, index=False, encoding="utf-8-sig")
    print("saved ->", out)


if __name__ == "__main__":
    main()
