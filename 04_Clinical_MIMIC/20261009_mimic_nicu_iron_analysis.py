# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""MIMIC-III NICU: infection x iron-restricted erythropoiesis -- observational test (filling the clinical side of the argument)

Hypothesis (derived from the genetic conclusion in section 15):
  infection -> IL-6 up -> hepcidin(HAMP) up -> iron utilisation restricted -> iron-restricted erythropoiesis
  => clinically visible: **MCV down, MCH down, RDW up, reticulocytes down, Hb down**
  If the reverse holds (infection -> bilirubin down driven by "iron restriction -> less haem turnover"), then **total bilirubin down** should appear at the same time (the validation endpoint of this project).

Two complementary designs:
  A. Cross-sectional (primary): in NICU neonates, compare the erythrocyte-index means on postnatal days 3-14 between the "positive-culture / negative-culture" groups;
     adjust for **birth weight** and **postnatal age** (MIMIC-III has no numerical gestational age, so weight serves as a maturity proxy), and stratify by weight.
  B. Pre/post-culture pairing (secondary): anchoring on each blood culture, take the nearest erythrocyte indices / reticulocytes / bilirubin within 72 h before and after,
     and compare the delta between "positive vs negative cultures" -- the negative-culture group naturally controls for postnatal age, care intensity and number of blood draws.

* Three points that must go into the limitations (write them down first, don't add them afterwards):
  1. Direct iron-metabolism assays are almost zero in neonates (Iron 5 rows, Ferritin 7) -> one can only detour through erythrocyte indices; **not a direct iron marker**
  2. **CoNS accounts for 53% of positive cultures** (a typical NICU contaminant) -> a "CoNS excluded" sensitivity analysis is mandatory
  3. The biological response of MCV/MCH/RDW is measured in **weeks**; the 72 h window is only meaningful for **reticulocytes/bilirubin**
"""
import os, sys, math
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

DIR = r"D:\workbuddy_workspace\03-data\MIMIC_NICU_iron"
OUT = DIR

CONT = ["mcv", "mch", "mchc", "rdw", "hb", "hct", "retic", "bili_total"]
COMS = {"staphylococcus": 1, "coagulase negative": 1}   # contaminant-prone organisms


def load():
    base = pd.read_csv(os.path.join(DIR, "nicu_base.csv"), encoding="utf-8-sig", parse_dates=["dob", "admittime", "intime", "outtime"])
    lab = pd.read_csv(os.path.join(DIR, "nicu_cbc.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    cul = pd.read_csv(os.path.join(DIR, "nicu_culture.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    return base, lab, cul


def build_birthweight():
    w = pd.read_csv(os.path.join(DIR, "nicu_weight.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    b = w[w.label.astype(str).str.contains("Birth Weight", case=False, na=False)]
    b = b[(b.valuenum > 0.2) & (b.valuenum < 7)]
    bw = b.sort_values("charttime").groupby("hadm_id").valuenum.first().rename("bw_kg")
    return bw


def main():
    base, lab, cul = load()
    bw = build_birthweight()

    # ---- Cohort: age at ICU entry <= 28 days ----
    base = base[base.age_days_at_icu <= 28].copy()
    base = base.merge(bw, left_on="hadm_id", right_index=True, how="left")
    base["bw_lt1500"] = (base.bw_kg < 1.5).astype("float")
    base.loc[base.bw_kg.isna(), "bw_lt1500"] = np.nan
    print("=" * 108)
    print("(1) Cohort (NICU and age at ICU entry <=28 days)")
    print("=" * 108)
    print("   hadm = %d; with birth weight = %d (%.1f%%); birth weight <1500 g = %d"
          % (len(base), base.bw_kg.notna().sum(), 100 * base.bw_kg.notna().mean(),
             int((base.bw_lt1500 == 1).sum())))
    print("   birth weight kg: median %.2f  IQR %.2f-%.2f"
          % (base.bw_kg.median(), base.bw_kg.quantile(.25), base.bw_kg.quantile(.75)))

    # ---- Exposure ----
    cul = cul[cul.hadm_id.isin(base.hadm_id)]
    pos = cul[cul.positive == 1]
    pos_orgs = ", ".join(pos.orgs.dropna().astype(str))
    pos_clean = pos[~pos.orgs.astype(str).str.lower().str.contains("staphylococcus|coagulase negative", na=False)]
    exp = pd.DataFrame(index=base.hadm_id.values)
    exp["any_culture"] = base.set_index("hadm_id").index.isin(cul.hadm_id).astype(int)
    exp["pos_blood"] = base.set_index("hadm_id").index.isin(pos.hadm_id).astype(int)
    exp["pos_blood_noncons"] = base.set_index("hadm_id").index.isin(pos_clean.hadm_id).astype(int)
    print("\n   exposure distribution (hadm): blood culture done %d | positive %d (%.1f%%) | positive and non-CoNS %d"
          % (exp.any_culture.sum(), exp.pos_blood.sum(),
             100 * exp.pos_blood.mean(), exp.pos_blood_noncons.sum()))
    print("   example of the positive organism string: %s" % pos_orgs[:180])

    # ---- Merge CBC ----
    lab = lab[lab.hadm_id.isin(base.hadm_id)].copy()
    lab = lab.merge(base[["hadm_id", "dob"]], on="hadm_id", how="left")
    lab["dol"] = (lab.charttime - lab.dob).dt.total_seconds() / 86400.0
    lab = lab[(lab.dol >= 0) & (lab.dol <= 60)]

    # Reticulocytes: merge automated + manual (two methods for the same variable)
    ret = lab[lab["var"].isin(["retic_auto", "retic_man"])].copy()
    ret["var"] = "retic"

    # =====================================================================
    print("\n" + "=" * 108)
    print("(2) Design A (primary): erythrocyte-index means on postnatal days 3-14 -- infection vs no infection")
    print("=" * 108)
    win = lab[(lab.dol >= 3) & (lab.dol <= 14) & (lab["var"].isin(["mcv", "mch", "mchc", "rdw", "hb", "hct"]))]
    agg = win.groupby(["hadm_id", "var"]).agg(val=("valuenum", "mean"), n=("valuenum", "size"),
                                              dol=("dol", "mean")).reset_index()
    wide = agg.pivot(index="hadm_id", columns="var", values="val")
    wide["dol_mean"] = agg.groupby("hadm_id").dol.mean()
    wide = wide.join(exp).join(base.set_index("hadm_id")[["bw_kg", "bw_lt1500"]])
    wide = wide.dropna(subset=["mcv", "mch", "rdw", "bw_kg"])
    print("   hadm available for the cross-sectional analysis = %d (positive %d / non-CoNS positive %d)"
          % (len(wide), int(wide.pos_blood.sum()), int(wide.pos_blood_noncons.sum())))

    rows = []
    for y in ["mcv", "mch", "mchc", "rdw", "hb", "hct"]:
        if y not in wide.columns:
            continue
        for xname, xcol in [("Positive blood culture", "pos_blood"), ("Positive blood culture (non-CoNS)", "pos_blood_noncons")]:
            d = wide.dropna(subset=[y, xcol, "bw_kg", "dol_mean"]).copy()
            if d[xcol].sum() < 10:
                continue
            # Crude model
            m1 = smf.ols("%s ~ %s" % (y, xcol), data=d).fit(cov_type="HC3")
            # Adjust for birth weight + postnatal age
            m2 = smf.ols("%s ~ %s + bw_kg + dol_mean" % (y, xcol), data=d).fit(cov_type="HC3")
            b, p = m2.params[xcol], m2.pvalues[xcol]
            ci = m2.conf_int().loc[xcol]
            rows.append(dict(outcome=y, exposure=xname, n=int(len(d)), ncase=int(d[xcol].sum()),
                             crude_b=m1.params[xcol], crude_p=m1.pvalues[xcol],
                             adj_b=b, adj_lo=ci[0], adj_hi=ci[1], adj_p=p,
                             mean_case=d.loc[d[xcol] == 1, y].mean(), mean_ctrl=d.loc[d[xcol] == 0, y].mean()))
    A = pd.DataFrame(rows)
    lab_map = {"mcv": "MCV (fL)", "mch": "MCH (pg)", "mchc": "MCHC (g/dL)", "rdw": "RDW (%)",
               "hb": "Hb (g/dL)", "hct": "Hct (%)"}
    print("\n   Expected directions (hepcidin-iron-restriction model): MCV down  MCH down  MCHC down(-)  RDW up  Hb down  Hct down")
    print("   %-11s %-20s %5s %5s %9s %9s %22s %9s %6s" % ("Outcome", "Exposure", "n", "cases", "case mean", "control mean", "adjusted β (95%CI)", "p", "crude p"))
    print("   " + "-" * 104)
    for _, r in A.iterrows():
        print("   %-11s %-20s %5d %5d %9.3f %9.3f  %+8.3f (%+.3f,%+.3f) %9.3g %6.3f"
              % (lab_map.get(r.outcome, r.outcome), r.exposure, r.n, r.ncase,
                 r.mean_case, r.mean_ctrl, r.adj_b, r.adj_lo, r.adj_hi, r.adj_p, r.crude_p))
    A.to_csv(os.path.join(OUT, "A_crosssection_results.csv"), index=False, encoding="utf-8-sig")

    # ---- Stratified by weight (<1500 g) ----
    print("\n   Stratified (birth weight <1500 g subgroup, by positive blood culture):")
    for y in ["mcv", "mch", "rdw"]:
        d = wide.dropna(subset=[y, "bw_lt1500", "dol_mean"])
        d = d[d.bw_lt1500 == 1]
        if d.pos_blood.sum() < 5 or (d.pos_blood == 0).sum() < 5:
            print("      %-6s subgroup too small (cases %d) -> skip" % (y, int(d.pos_blood.sum()))); continue
        m = smf.ols("%s ~ pos_blood + bw_kg + dol_mean" % y, data=d).fit(cov_type="HC3")
        print("      %-6s n=%3d cases=%2d  adjusted beta=%+.3f (p=%.3f)"
              % (y, len(d), int(d.pos_blood.sum()), m.params.pos_blood, m.pvalues.pos_blood))

    # =====================================================================
    print("\n" + "=" * 108)
    print("(3) Design B (secondary): 72 h before/after pairing anchored on each blood culture (including reticulocytes and bilirubin)")
    print("=" * 108)
    pair_rows = []
    series = pd.concat([lab, ret], ignore_index=True)
    for var in ["retic", "bili_total", "hb", "mcv", "mch", "rdw"]:
        s = series[series["var"] == var][["hadm_id", "charttime", "valuenum"]].sort_values("charttime")
        for _, c in cul.iterrows():
            g = s[s.hadm_id == c.hadm_id]
            if not len(g):
                continue
            pre = g[(g.charttime <= c.charttime) & (g.charttime >= c.charttime - pd.Timedelta(hours=48))]
            post = g[(g.charttime > c.charttime) & (g.charttime <= c.charttime + pd.Timedelta(hours=72))]
            if not len(pre) or not len(post):
                continue
            pair_rows.append(dict(var=var, hadm_id=c.hadm_id, positive=c.positive,
                                  delta=post.valuenum.iloc[0] - pre.valuenum.iloc[-1],
                                  pre=pre.valuenum.iloc[-1], post=post.valuenum.iloc[0],
                                  gap_h=float((post.charttime.iloc[0] - pre.charttime.iloc[-1]).total_seconds() / 3600)))
    P = pd.DataFrame(pair_rows)
    print("   pairs %d; by index:" % len(P))
    resB = []
    for var in ["retic", "bili_total", "hb", "mcv", "mch", "rdw"]:
        d = P[P["var"] == var].dropna(subset=["delta"])
        if len(d) < 30 or d.positive.sum() < 8:
            print("      %-12s pairs %d (positive %d) -> insufficient" % (var, len(d), int(d.positive.sum()) if len(d) else 0)); continue
        m1 = smf.ols("delta ~ positive", data=d).fit(cov_type="HC3")
        m2 = smf.ols("delta ~ positive + gap_h + pre", data=d).fit(cov_type="HC3")
        ci = m2.conf_int().loc["positive"]
        print("      %-12s pairs %4d (positive %3d)  Delta positive-negative = %+.4f (%+.4f,%+.4f) p=%.3g | crude %+.4f (p=%.3f)"
              % (var, len(d), int(d.positive.sum()), m2.params.positive, ci[0], ci[1],
                 m2.pvalues.positive, m1.params.positive, m1.pvalues.positive))
        resB.append(dict(var=var, n=len(d), ncase=int(d.positive.sum()),
                         adj_b=m2.params.positive, adj_lo=ci[0], adj_hi=ci[1], adj_p=m2.pvalues.positive,
                         crude_b=m1.params.positive, crude_p=m1.pvalues.positive,
                         mean_delta_pos=d.loc[d.positive == 1, "delta"].mean(),
                         mean_delta_neg=d.loc[d.positive == 0, "delta"].mean()))
    B = pd.DataFrame(resB)
    if len(B):
        B.to_csv(os.path.join(OUT, "B_within_subject_results.csv"), index=False, encoding="utf-8-sig")
    P.to_csv(os.path.join(OUT, "B_pairs_long.csv"), index=False, encoding="utf-8-sig")

    print("\n   Expected directions: within 72 h around a positive culture, **reticulocytes down** and **bilirubin down** (the latter being the pipeline validation of an existing conclusion of this project)")


if __name__ == "__main__":
    main()
