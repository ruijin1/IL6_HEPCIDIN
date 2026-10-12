# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""MIMIC-III NICU: infection x erythrocyte indices -- robustness check **excluding transfusion confounding**

Why this step is mandatory:
  The main analysis (`20261009_mimic_nicu_iron_analysis.py`) found in the non-CoNS positive blood-culture group
  **MCV -4.82 (p=0.0004), MCH -1.65 (p=0.0001), Hct -1.45 (p=0.027)**.
  But **transfusion** is the number-one confounder of MCV: adult donor blood MCV~=90 fL vs neonate~=104 fL,
  infected neonates are transfused more -> transfusion alone can produce a fall in MCV.
  Rough estimate: replacing ~40 mL of adult blood in a 2.5 kg neonate (~=16% of blood volume) can lower MCV by about **2 fL**, the same order as observed.

Triple robustness:
  (1) `transfused_any` (ICD-9 9904 red-cell transfusion / 9901 exchange transfusion) as an adjustment variable
  (2) **using only non-transfused neonates** (sensitivity)
  (3) * using only labs **before the "transfusion consent" record** (charttime of chartevents itemid 3667 as a transfusion-time proxy)
  Plus a **dose-response**: number of positive culture episodes 0 / 1 / >=2
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
import psycopg2
import statsmodels.formula.api as smf

DIR = r"D:\workbuddy_workspace\03-data\MIMIC_NICU_iron"
CONN = dict(host="localhost", port=5432, user="postgres", password="1234", dbname="mimiciii")


def main():
    c = psycopg2.connect(**CONN)
    base = pd.read_csv(os.path.join(DIR, "nicu_base.csv"), encoding="utf-8-sig",
                       parse_dates=["dob", "admittime", "intime", "outtime"])
    lab = pd.read_csv(os.path.join(DIR, "nicu_cbc.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    cul = pd.read_csv(os.path.join(DIR, "nicu_culture.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    w = pd.read_csv(os.path.join(DIR, "nicu_weight.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    base = base[base.age_days_at_icu <= 28].copy()
    ids = tuple(base.hadm_id.astype(int).tolist())
    ids_sql = ",".join(str(x) for x in ids)

    # ---- Birth weight ----
    b = w[w.label.astype(str).str.contains("Birth Weight", case=False, na=False)]
    b = b[(b.valuenum > 0.2) & (b.valuenum < 7)]
    bw = b.sort_values("charttime").groupby("hadm_id").valuenum.first().rename("bw_kg")

    # ---- Transfusion ----
    pr = pd.read_sql("select hadm_id, icd9_code from mimiciii.procedures_icd where hadm_id in (%s)" % ids_sql, c)
    transf_icd = set(pr[pr.icd9_code.isin(["9904", "9901"])].hadm_id)
    tx = pd.read_sql("""select hadm_id, charttime from mimiciii.chartevents
                        where hadm_id in (%s) and itemid=3667 and charttime is not null""" % ids_sql, c)
    tx["charttime"] = pd.to_datetime(tx.charttime)
    txflag = tx.groupby("hadm_id").charttime.min().rename("tx_consent_time")
    print("=" * 108)
    print("(1) Transfusion availability")
    print("=" * 108)
    print("   hadm flagged by ICD-9 9904/9901 (red-cell transfusion/exchange transfusion): %d (%.1f%%)"
          % (len(transf_icd & set(ids)), 100 * len(transf_icd & set(ids)) / len(ids)))
    print("   hadm with a 'Transfusion Consent' charttime: %d" % len(txflag))
    print("   either = %d" % len((transf_icd & set(ids)) | set(txflag.index)))

    # ---- Exposure ----
    cul = cul[cul.hadm_id.isin(base.hadm_id)]
    pos = cul[cul.positive == 1]
    pos_clean = pos[~pos.orgs.astype(str).str.lower().str.contains("staphylococcus|coagulase negative", na=False)]
    npos = pos.groupby("hadm_id").size().rename("n_pos")

    # ---- Labs ----
    lab = lab[lab.hadm_id.isin(base.hadm_id)].copy()
    lab = lab.merge(base[["hadm_id", "dob"]], on="hadm_id")
    lab["dol"] = (lab.charttime - lab.dob).dt.total_seconds() / 86400.0
    lab = lab[(lab.dol >= 3) & (lab.dol <= 14) & lab["var"].isin(["mcv", "mch", "rdw", "hb", "hct", "mchc"])]
    lab = lab.merge(txflag, left_on="hadm_id", right_index=True, how="left")
    lab["pre_tx"] = (lab.tx_consent_time.isna()) | (lab.charttime < lab.tx_consent_time)

    agg = lab.groupby(["hadm_id", "var"]).agg(val=("valuenum", "mean"), dol=("dol", "mean")).reset_index()
    wide = agg.pivot(index="hadm_id", columns="var", values="val")
    wide["dol_mean"] = agg.groupby("hadm_id").dol.mean()
    wide = wide.join(bw).join(npos).join(txflag)
    wide["n_pos"] = wide.n_pos.fillna(0).astype(int)
    wide["pos_blood"] = (wide.n_pos > 0).astype(int)
    wide["pos_noncons"] = wide.index.isin(pos_clean.hadm_id).astype(int)
    wide["transfused"] = wide.index.isin(transf_icd | set(txflag.index)).astype(int)
    wide["any_icd_tx"] = wide.index.isin(transf_icd).astype(int)
    wide["any_tx_consent"] = wide.tx_consent_time.notna().astype(int)

    # Use only pre-transfusion labs (re-aggregated within the same hadm)
    labpre = lab[lab.pre_tx]
    agg2 = labpre.groupby(["hadm_id", "var"]).agg(val=("valuenum", "mean")).reset_index()
    wide_pre = agg2.pivot(index="hadm_id", columns="var", values="val")
    wide_pre.columns = [str(x) + "_preTx" for x in wide_pre.columns]

    D = wide.join(wide_pre)
    D = D.dropna(subset=["mcv", "mch", "rdw", "bw_kg"])
    print("\n   hadm available for analysis = %d (positive %d, non-CoNS positive %d, transfused %d)"
          % (len(D), D.pos_blood.sum(), D.pos_noncons.sum(), D.transfused.sum()))
    print("   relation between infection and transfusion: transfusion rate %.0f%% in the positive group vs %.0f%% in the negative group"
          % (100 * D.loc[D.pos_blood == 1, "transfused"].mean(),
             100 * D.loc[D.pos_blood == 0, "transfused"].mean()))

    def fit(data, y, x, extra=""):
        d = data.dropna(subset=[y, x, "bw_kg", "dol_mean"])
        if len(d) < 40 or d[x].sum() < 8:
            return None
        m = smf.ols("%s ~ %s + bw_kg + dol_mean %s" % (y, x, extra), data=d).fit(cov_type="HC3")
        ci = m.conf_int().loc[x]
        return dict(outcome=y, exposure=x, n=len(d), ncase=int(d[x].sum()),
                    b=m.params[x], lo=ci[0], hi=ci[1], p=m.pvalues[x])

    print("\n" + "=" * 108)
    print("(2) Triple robustness (MCV / MCH / RDW; adjusted for birth weight + postnatal age)")
    print("=" * 108)
    rows = []
    for y in ["mcv", "mch", "rdw", "hct"]:
        for x in ["pos_blood", "pos_noncons"]:
            r = fit(D, y, x); 
            if r: r["model"] = "Main (all)"; rows.append(r)
            r = fit(D, y, x, "+ transfused")
            if r: r["model"] = "(1) adjusted for transfusion"; rows.append(r)
            r = fit(D[D.transfused == 0], y, x)
            if r: r["model"] = "(2) never transfused only"; rows.append(r)
        # (3) Use only pre-transfusion labs (x uses the main exposure, y uses _preTx)
        for x in ["pos_blood", "pos_noncons"]:
            yp = y + "_preTx"
            if yp in D.columns:
                r = fit(D, yp, x)
                if r: r["model"] = "(3) pre-transfusion measurements only"; r["outcome"] = y + "(preTx)"; rows.append(r)
    R = pd.DataFrame(rows)
    print("%-11s %-16s %-18s %5s %5s %22s %9s" % ("Outcome", "Exposure", "Model", "n", "cases", "β (95%CI)", "p"))
    print("-" * 108)
    for _, r in R.iterrows():
        print("%-11s %-16s %-18s %5d %5d  %+7.3f (%+.3f,%+.3f) %9.3g"
              % (r.outcome, r.exposure, r.model, r.n, r.ncase, r.b, r.lo, r.hi, r.p))
    R.to_csv(os.path.join(DIR, "C_transfusion_robust_results.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 108)
    print("(3) Dose-response (number of positive culture episodes 0 / 1 / >=2), adjusted for birth weight + postnatal age")
    print("=" * 108)
    D["npos_cat"] = np.where(D.n_pos >= 2, 2, D.n_pos)
    for y in ["mcv", "mch", "rdw", "hct"]:
        d = D.dropna(subset=[y, "bw_kg", "dol_mean"])
        try:
            m = smf.ols("%s ~ C(npos_cat) + bw_kg + dol_mean" % y, data=d).fit(cov_type="HC3")
        except Exception as e:
            print("  %-6s failed %s" % (y, e)); continue
        cells = []
        for k in [1, 2]:
            key = "C(npos_cat)[T.%d]" % k
            if key in m.params.index:
                cells.append("%d cultures β=%+.2f(p=%.3g)" % (k, m.params[key], m.pvalues[key]))
        grp = d.groupby("npos_cat")[y].agg(["mean", "size"])
        print("  %-6s group mean 0 cultures=%.2f(n=%d) 1 culture=%.2f(n=%d) >=2 cultures=%.2f(n=%d) | %s"
              % (y, grp.loc[0, "mean"], grp.loc[0, "size"],
                 grp.loc[1, "mean"] if 1 in grp.index else float("nan"), grp.loc[1, "size"] if 1 in grp.index else 0,
                 grp.loc[2, "mean"] if 2 in grp.index else float("nan"), grp.loc[2, "size"] if 2 in grp.index else 0,
                 "  ".join(cells)))

    print("\n* Interpretation:")
    print("  (1) If the decline in MCV and MCH **persists** after adjusting for / excluding transfusion -> it cannot be attributed entirely to transfusion, supporting iron-restricted erythropoiesis")
    print("  (2) If it is significant only in the main model and disappears in all others -> **honestly write it up as 'explained by transfusion'**, and do not keep the conclusion")
    print("  (3) The absence of any change in RDW (which was expected to rise) must be stated honestly in every version")


if __name__ == "__main__":
    main()
