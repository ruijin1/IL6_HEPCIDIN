# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""MIMIC-III NICU: culture-positive -> bilirubin decline -- robustness check **controlling for phototherapy**

Background: the project's **first observational finding came from MIMIC-IV** (culture-positive -> bilirubin -0.355, p=6x10^-59),
      but a key confounder was found there: **phototherapy rate 59.6% vs 18.3% (3.3-fold) in the infection group**,
      i.e. the "lower bilirubin" may be partly **treatment-induced**.
      This section does the same thing **completely independently in the MIMIC-III NICU** (72 h before/after pairing),
      and controls this confounder directly using **chartevents phototherapy lamp/blanket (3570/3569, with charttime)**.

Quantifying the confounder:
  - How much the phototherapy rate differs between positive vs negative cultures
  - Keep only pairs with "**no** phototherapy record within the +/-72 h window" (cleanest)
  - Adjust for "phototherapy within the window"
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
    base = pd.read_csv(os.path.join(DIR, "nicu_base.csv"), encoding="utf-8-sig", parse_dates=["dob", "intime"])
    lab = pd.read_csv(os.path.join(DIR, "nicu_cbc.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    cul = pd.read_csv(os.path.join(DIR, "nicu_culture.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    base = base[base.age_days_at_icu <= 28]
    ids = tuple(base.hadm_id.astype(int).tolist())
    ids_sql = ",".join(str(x) for x in ids)

    # Phototherapy (chartevents lamp/blanket) + ICD-9 9983
    ph = pd.read_sql("""select hadm_id, charttime, itemid from mimiciii.chartevents
                        where hadm_id in (%s) and itemid in (3569,3570) and charttime is not null""" % ids_sql, c)
    ph["charttime"] = pd.to_datetime(ph.charttime)
    pr = pd.read_sql("select hadm_id, icd9_code from mimiciii.procedures_icd where hadm_id in (%s)" % ids_sql, c)
    ph_icd = set(pr[pr.icd9_code == "9983"].hadm_id)
    print("=" * 108)
    print("(1) Phototherapy availability (NICU, age at ICU entry <=28d)")
    print("=" * 108)
    print("   chartevents phototherapy records: hadm %d; ICD-9 9983: hadm %d; either: %d"
          % (ph.hadm_id.nunique(), len(ph_icd & set(ids)), len(set(ph.hadm_id) | (ph_icd & set(ids)))))

    # Bilirubin time series
    bili = lab[(lab["var"] == "bili_total")][["hadm_id", "charttime", "valuenum"]].sort_values("charttime")
    cul = cul[cul.hadm_id.isin(base.hadm_id)]

    rows = []
    for _, cu in cul.iterrows():
        g = bili[bili.hadm_id == cu.hadm_id]
        if not len(g):
            continue
        t0, t1 = cu.charttime, cu.charttime + pd.Timedelta(hours=72)
        pre = g[(g.charttime <= t0) & (g.charttime >= t0 - pd.Timedelta(hours=48))]
        post = g[(g.charttime > t0) & (g.charttime <= t1)]
        if not len(pre) or not len(post):
            continue
        # Phototherapy within the window
        hp = ph[(ph.hadm_id == cu.hadm_id) & (ph.charttime >= t0) & (ph.charttime <= t1)]
        hp_pre = ph[(ph.hadm_id == cu.hadm_id) & (ph.charttime >= t0 - pd.Timedelta(hours=48)) & (ph.charttime <= t1)]
        rows.append(dict(hadm_id=cu.hadm_id, positive=cu.positive,
                         delta=post.valuenum.iloc[0] - pre.valuenum.iloc[-1],
                         pre=pre.valuenum.iloc[-1], post=post.valuenum.iloc[0],
                         gap_h=(post.charttime.iloc[0] - pre.charttime.iloc[-1]).total_seconds() / 3600,
                         ph_in_win=int(len(hp) > 0), ph_around=int(len(hp_pre) > 0),
                         any_ph_icd=int(cu.hadm_id in ph_icd),
                         # * whether this pair's own culture is CoNS (for pair-level stratification, see (5)/(6))
                         cu_is_cons=int(('staphylococcus' in str(cu.orgs).lower())
                                        or ('coagulase negative' in str(cu.orgs).lower()))
                         ))
    P = pd.DataFrame(rows)
    print("\n   pairs = %d (positive %d)" % (len(P), int(P.positive.sum())))
    print("\n" + "=" * 108)
    print("(2) * Quantifying the phototherapy confounder")
    print("=" * 108)
    print("   Positive culture: phototherapy rate within window %.1f%% (phototherapy around %.1f%%)"
          % (100 * P.loc[P.positive == 1, "ph_in_win"].mean(), 100 * P.loc[P.positive == 1, "ph_around"].mean()))
    print("   Negative culture: phototherapy rate within window %.1f%% (phototherapy around %.1f%%)"
          % (100 * P.loc[P.positive == 0, "ph_in_win"].mean(), 100 * P.loc[P.positive == 0, "ph_around"].mean()))
    r = P.groupby("positive").ph_in_win.mean()
    print("   => the phototherapy rate in the positive group is **%.1f-fold** that in the negative group" % (r.get(1, np.nan) / max(r.get(0, np.nan), 1e-9)))

    print("\n" + "=" * 108)
    print("(3) Delta bilirubin (post - pre, mg/dL) -- controlling for phototherapy layer by layer")
    print("=" * 108)
    out = []
    def rep(tag, d, formula):
        # * 2026-10-10: rows with too few samples must still be written to the results table, otherwise in the main-text Table 1
        #   the row "(3) admission never receiving phototherapy (ICD) only, n=125 (positive 5), not estimable"
        #   would have no traceable record in the CSV for reviewers. Changed to record n / ncase + an empty estimate.
        if len(d) < 30 or d.positive.sum() < 8:
            print("   %-40s n=%d (positive %d) -> insufficient" % (tag, len(d), int(d.positive.sum()) if len(d) else 0))
            out.append(dict(model=tag, n=len(d), ncase=int(d.positive.sum()) if len(d) else 0,
                            b=np.nan, lo=np.nan, hi=np.nan, p=np.nan, estimable=0))
            return
        m = smf.ols(formula, data=d).fit(cov_type="HC3")
        ci = m.conf_int().loc["positive"]
        print("   %-40s n=%4d (positive %3d)  delta-diff = %+.4f (%+.4f,%+.4f)  p=%.3g"
              % (tag, len(d), int(d.positive.sum()), m.params.positive, ci[0], ci[1], m.pvalues.positive))
        out.append(dict(model=tag, n=len(d), ncase=int(d.positive.sum()), b=m.params.positive,
                        lo=ci[0], hi=ci[1], p=m.pvalues.positive, estimable=1))
    rep("Main model (positive vs negative)", P, "delta ~ positive")
    rep("+ adjusted for interval and baseline value", P, "delta ~ positive + gap_h + pre")
    rep("(1) adjusted for phototherapy in window", P, "delta ~ positive + ph_in_win + gap_h + pre")
    rep("(2) no phototherapy in window only", P[P.ph_in_win == 0], "delta ~ positive + gap_h + pre")
    # * 2026-10-10: the original "(4) non-CoNS positive vs negative" row was removed. That row `P.assign(positive=lambda d: d.positive)`
    #   was a no-op, making (4) numerically identical to "+ adjusted for interval and baseline value"; it was a deprecated duplicate
    #   superseded by (5)/(6), and was not referenced in the main-text Table 1. After removal the 7 CSV rows
    #   correspond one-to-one with the 7 rows in the main-text Table 1.
    rep("(3) admission never receiving phototherapy (ICD) only", P[P.any_ph_icd == 0], "delta ~ positive + gap_h + pre")

    # Positive stratification: non-CoNS -- * must be at the **pair level** (the pathogen of that pair's own culture),
    # consistent with the exposure level of the other rows of this table (main model/(1)/(2)); using the admission level would
    # count pairs with negative cultures as "cases" and would make the sum of the "definite pathogen" and "CoNS only" rows
    # (91+42=133) exceed the number of positive pairs (65) -- self-contradictory.
    P2 = P.copy()
    P2["pos_definite"] = ((P2.positive == 1) & (P2.cu_is_cons == 0)).astype(int)
    P2["pos_cons_only"] = ((P2.positive == 1) & (P2.cu_is_cons == 1)).astype(int)
    assert int(P2.pos_definite.sum() + P2.pos_cons_only.sum()) == int(P2.positive.sum()), \
        "the case count of (5)+(6) should equal the number of positive pairs"
    for tag, x in [("(5) definite pathogen (non-CoNS) positive vs negative", "pos_definite"), ("(6) CoNS-only positive vs negative", "pos_cons_only")]:
        d = P2.dropna(subset=["delta", "gap_h", "pre"])
        m = smf.ols("delta ~ %s + gap_h + pre" % x, data=d).fit(cov_type="HC3")
        ci = m.conf_int().loc[x]
        print("   %-40s n=%4d (cases %3d)  delta-diff = %+.4f (%+.4f,%+.4f)  p=%.3g"
              % (tag, len(d), int(d[x].sum()), m.params[x], ci[0], ci[1], m.pvalues[x]))
        out.append(dict(model=tag, n=len(d), ncase=int(d[x].sum()), b=m.params[x],
                        lo=ci[0], hi=ci[1], p=m.pvalues[x], estimable=1))

    pd.DataFrame(out).to_csv(os.path.join(DIR, "D_bilirubin_phototherapy_robust.csv"), index=False, encoding="utf-8-sig")
    P.to_csv(os.path.join(DIR, "D_bilirubin_pairs.csv"), index=False, encoding="utf-8-sig")
    print("\n* Interpretation: if 'no phototherapy in window only' remains significantly negative -> infection -> bilirubin down is not explained by phototherapy; if it disappears -> the interpretation must honestly be changed to phototherapy.")


if __name__ == "__main__":
    main()
