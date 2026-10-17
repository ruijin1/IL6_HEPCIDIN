# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MIMIC_NICU_iron). Analysis logic is unchanged.
r"""MIMIC-III NICU: within-subject bilirubin decline -- robustness check **controlling for red-cell transfusion**

Why this step matters:
  The within-subject analysis (`20261009_mimic_nicu_bilirubin_phototherapy_robust.py`) found that in the 72 h after a
  positive blood culture bilirubin fell more than after a negative one (-1.07 mg/dL, 95% CI -1.66 to -0.47).
  Infection is associated with a much higher transfusion rate, and adult donor blood **dilutes** neonatal bilirubin.
  Transfusion is therefore the leading competing explanation that has NOT yet been excluded for bilirubin
  (it *was* excluded for the erythrocyte indices, in `20261009_mimic_nicu_iron_transfusion_robust.py`).

Design:
  Step 0  Rebuild the 728 culture-anchored pairs and verify they reproduce the published CSV **row for row**.
  Step 1  Build transfusion variables from two sources:
            (a) ICD-9 procedure codes 9904/9901/9903/9905/9906/9907 (red-cell / exchange transfusion) -- no date
            (b) chartevents itemid 3667 "Transfusion Consent" -- carries charttime, used as the time proxy
  Step 2  Layer the confounder in, mirroring the phototherapy script:
            main              delta ~ positive + gap_h + pre
            (1) adjust        + tx_any          (admission-level, either source)
            (2) adjust        + tx_in_win       (transfusion consent inside the -48 h / +72 h window)
            (3) subgroup      no transfusion in window
            (4) subgroup      never transfused (admission level, either source)
  Estimator: OLS with HC3 robust SEs, identical to the published phototherapy robustness script.
  Sample-size guard, identical to that script: n < 30 or positive < 8 -> recorded but not estimated.

Interpretation (written into the script before seeing the numbers, so the reading cannot be fitted afterwards):
  If the decline persists after adjusting for / excluding transfusion -> dilution does not explain it.
  If the decline is present only in the main model and vanishes in every layer -> it must honestly be
  written up as explained by transfusion, and the conclusion must be changed.
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
import psycopg2
import statsmodels.formula.api as smf

DIR = r"D:\workbuddy_workspace\03-data\MIMIC_NICU_iron"
CONN = dict(host="localhost", port=5432, user="postgres", password="1234", dbname="mimiciii")

TX_ICD = ("9904", "9901", "9903", "9905", "9906", "9907")   # red-cell / exchange / other transfusion
TX_ITEM = 3667                                             # 'Transfusion Consent' (time proxy)


def build_pairs():
    """Rebuild the culture-anchored before/after pairs exactly as the published script does."""
    base = pd.read_csv(os.path.join(DIR, "nicu_base.csv"), encoding="utf-8-sig",
                       parse_dates=["dob", "intime"])
    lab = pd.read_csv(os.path.join(DIR, "nicu_cbc.csv"), encoding="utf-8-sig",
                      parse_dates=["charttime"])
    cul = pd.read_csv(os.path.join(DIR, "nicu_culture.csv"), encoding="utf-8-sig",
                      parse_dates=["charttime"])
    base = base[base.age_days_at_icu <= 28]
    cul = cul[cul.hadm_id.isin(base.hadm_id)]
    bili = lab[lab["var"] == "bili_total"][["hadm_id", "charttime", "valuenum"]].sort_values("charttime")

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
        # pair-level pathogen class, exactly as the published script derives it
        cu_is_cons = int(("staphylococcus" in str(cu.orgs).lower())
                         or ("coagulase negative" in str(cu.orgs).lower()))
        rows.append(dict(hadm_id=cu.hadm_id, positive=int(cu.positive),
                         cu_is_cons=cu_is_cons,
                         t0=t0, t1=t1,
                         pre_t=pre.charttime.iloc[-1], post_t=post.charttime.iloc[0],
                         delta=post.valuenum.iloc[0] - pre.valuenum.iloc[-1],
                         pre=pre.valuenum.iloc[-1], post=post.valuenum.iloc[0],
                         gap_h=(post.charttime.iloc[0] - pre.charttime.iloc[-1]).total_seconds() / 3600))
    return pd.DataFrame(rows)


def transfusion():
    """Transfusion flags per hadm_id: any-source flag plus charttimes of the consent proxy."""
    c = psycopg2.connect(**CONN)
    icd = pd.read_sql("select hadm_id from mimiciii.procedures_icd where icd9_code in %s" % (TX_ICD,), c)
    ce = pd.read_sql("select hadm_id, charttime from mimiciii.chartevents "
                     "where itemid = %s and charttime is not null" % TX_ITEM, c)
    c.close()
    return set(icd.hadm_id.astype(int)), ce.assign(charttime=pd.to_datetime(ce.charttime))


def main():
    P = build_pairs()

    # ---- Step 0: verify the rebuild reproduces the published 728 pairs -------------------------
    ref = pd.read_csv(os.path.join(DIR, "D_bilirubin_pairs.csv"), encoding="utf-8-sig")
    print("=" * 108)
    print("(0) Reproduction check against the published D_bilirubin_pairs.csv")
    print("=" * 108)
    print("   rebuilt pairs = %d ; published = %d" % (len(P), len(ref)))
    # Note: a simple merge over-counts here (770 > 728) because several admissions contribute
    # identical (positive, pre, post) tuples. Compare the two multisets row by row instead.
    key = ["hadm_id", "positive", "pre", "post", "delta", "gap_h"]
    a = ref[key].sort_values(key).reset_index(drop=True)
    b = P[key].sort_values(key).reset_index(drop=True)
    ok = len(a) == len(b) and np.allclose(a.select_dtypes("number"), b.select_dtypes("number"),
                                          atol=1e-9, equal_nan=True)
    print("   row-by-row multiset comparison on %s" % ", ".join(key))
    print("   max |delta difference| = %.2e" % np.nanmax(np.abs(a.delta - b.delta)))
    print("   max |gap_h difference| = %.2e" % np.nanmax(np.abs(a.gap_h - b.gap_h)))
    assert ok, "rebuild did not reproduce the published pairs"
    print("   -> reproduction OK (%d rows, identical to floating-point tolerance)\n" % len(a))

    # ---- Step 1: transfusion variables ---------------------------------------------------------
    tx_icd, ce = transfusion()
    ids = set(P.hadm_id.astype(int))
    P["tx_icd"] = P.hadm_id.isin(tx_icd).astype(int)
    ce = ce[ce.hadm_id.isin(ids)]
    by_hadm = {h: g.charttime.values for h, g in ce.groupby("hadm_id")}

    def in_win(r):
        t = by_hadm.get(r.hadm_id)
        if t is None:
            return 0
        lo = np.datetime64(r.t0 - pd.Timedelta(hours=48))
        hi = np.datetime64(r.t1)
        return int(np.any((t >= lo) & (t <= hi)))
    P["tx_consent_any"] = P.hadm_id.isin(set(ce.hadm_id)).astype(int)
    P["tx_in_win"] = P.apply(in_win, axis=1)
    P["tx_any"] = ((P.tx_icd == 1) | (P.tx_consent_any == 1)).astype(int)

    print("=" * 108)
    print("(1) Transfusion availability and the size of the competing explanation")
    print("=" * 108)
    print("   pairs = %d (positive %d) ; admissions = %d" % (len(P), int(P.positive.sum()), P.hadm_id.nunique()))
    for col, tag in [("tx_icd", "ICD-9 transfusion code"),
                     ("tx_consent_any", "Transfusion Consent (itemid 3667)"),
                     ("tx_any", "either source"),
                     ("tx_in_win", "consent inside the -48 h / +72 h window")]:
        print("   %-42s %4d pairs (%.0f%%)" % (tag, int(P[col].sum()), 100 * P[col].mean()))
    for v, g in P.groupby("positive"):
        print("   %-14s n=%-4d  any transfusion %.0f%%  in-window %.0f%%"
              % ("positive" if v == 1 else "negative", len(g), 100 * g.tx_any.mean(), 100 * g.tx_in_win.mean()))
    r = P.groupby("positive").tx_any.mean()
    print("   => any-transfusion rate ratio positive/negative = %.2f" % (r.get(1, np.nan) / max(r.get(0, np.nan), 1e-9)))
    print()

    # ---- Step 2: layer the confounder in -------------------------------------------------------
    print("=" * 108)
    print("(2) Delta bilirubin (post - pre, mg/dL) -- controlling for transfusion layer by layer")
    print("=" * 108)
    out = []

    def rep(tag, d, formula, estimable_guard=True, term="positive"):
        d = d.dropna(subset=["delta", "gap_h", "pre"] + ([term] if term != "positive" else []))
        if estimable_guard and (len(d) < 30 or d.positive.sum() < 8):
            print("   %-46s n=%d (positive %d) -> insufficient" % (tag, len(d), int(d.positive.sum())))
            out.append(dict(model=tag, n=len(d), ncase=int(d.positive.sum()),
                            b=np.nan, lo=np.nan, hi=np.nan, p=np.nan, estimable=0))
            return
        m = smf.ols(formula, data=d).fit(cov_type="HC3")
        ci = m.conf_int().loc[term]
        print("   %-46s n=%4d (cases %4d)  %+.4f (%+.4f,%+.4f)  p=%.3g"
              % (tag, len(d), int(d[term].sum()), m.params[term], ci[0], ci[1], m.pvalues[term]))
        out.append(dict(model=tag, n=len(d), ncase=int(d[term].sum()), b=m.params[term],
                        lo=ci[0], hi=ci[1], p=m.pvalues[term], estimable=1))

    rep("Main model (positive vs negative)", P, "delta ~ positive")
    rep("+ adjusted for interval and baseline value", P, "delta ~ positive + gap_h + pre")
    rep("(1) + adjusted for transfusion (any source)", P, "delta ~ positive + tx_any + gap_h + pre")
    rep("(2) + adjusted for transfusion in window", P, "delta ~ positive + tx_in_win + gap_h + pre")
    rep("(3) no transfusion in window only", P[P.tx_in_win == 0], "delta ~ positive + gap_h + pre")
    rep("(4) never transfused (any source) only", P[P.tx_any == 0], "delta ~ positive + gap_h + pre")

    # Definite-pathogen stratum, matching Table 1 row (5) of the published analysis
    P["pos_definite"] = ((P.positive == 1) & (P.cu_is_cons == 0)).astype(int)
    print("   definite-pathogen pairs = %d (sanity: positive %d)" % (int(P.pos_definite.sum()), int(P.positive.sum())))
    rep("(5) definite pathogen, + adjusted for transfusion", P,
        "delta ~ pos_definite + tx_any + gap_h + pre", estimable_guard=False, term="pos_definite")

    R = pd.DataFrame(out)
    R.to_csv(os.path.join(DIR, "E_bilirubin_transfusion_robust.csv"), index=False, encoding="utf-8-sig")
    P.drop(columns=["t0", "t1", "pre_t", "post_t"]).to_csv(
        os.path.join(DIR, "E_bilirubin_pairs_transfusion.csv"), index=False, encoding="utf-8-sig")
    print("\n   saved: E_bilirubin_transfusion_robust.csv / E_bilirubin_pairs_transfusion.csv")


if __name__ == "__main__":
    main()
