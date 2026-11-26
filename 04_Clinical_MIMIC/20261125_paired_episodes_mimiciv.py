# -*- coding: utf-8 -*-
"""MIMIC-IV v1.0: within-subject (-48 h / +72 h) bilirubin change around every blood culture.

Replaces the earlier episode construction (`paired_replication_mimiciv.py`, 2026-10-24),
which iterated over `work.neo_cx` row by row. One row of `microbiologyevents` is one
organism x one antibiotic susceptibility result, and every culture also contributes a
Gram-stain row that carries no organism. Row-by-row iteration therefore counted a single
positive culture up to 20 times, split it across several "episodes", and mixed Gram-stain
rows into the negative group. The result was 2,261 pairs / 845 positive pairs, and a
positive episode was defined at a different level from the one the Methods describe.

The Methods define an episode as admission x collection time x specimen type, with
positivity assigned AFTER aggregation: a culture is positive when at least one organism is
reported, and Gram-stain rows (org_name IS NULL) are not cultures. This script implements
that definition, refits every number reported in Table S5 and in the Results, and writes a
pair-level dataset.

Design
  Cohort, exposure and outcome as in `20261008_neo_infection_bilirubin.py`.
  Each blood culture is an anchor irrespective of result. A pair is kept when a total
  bilirubin (itemid 50885) is available in the 48 h before and in the 72 h after collection;
  the last pre-culture and the first post-culture values are compared (delta = post - pre).
  Negative-culture pairs are the internal control.
  OLS of delta on culture positivity, adjusted for the sampling interval (gap_h), the
  baseline value (pre) and postnatal age at collection (age_d), with standard errors
  clustered by admission.

Usage
  python 20261125_paired_episodes_mimiciv.py [output_dir]
  Requires the local PostgreSQL database mimiciv10 with the work.* views.
"""
import os
import sys

import numpy as np
import pandas as pd
import psycopg2
import statsmodels.formula.api as smf

OUT = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()
os.makedirs(OUT, exist_ok=True)

con = psycopg2.connect(host="localhost", port=5432, user="postgres",
                       password="1234", dbname="mimiciv10")


def q(sql):
    return pd.read_sql(sql, con)


# ---------------------------------------------------------------- 1. episodes
# Every blood-culture row is read, including Gram-stain rows, so that the aggregation can
# be shown explicitly rather than assumed.
cx = q("""SELECT hadm_id, charttime, spec_type_desc, org_name
          FROM work.neo_cx
          WHERE spec_class = 'blood' AND charttime IS NOT NULL
            AND hadm_id IN (SELECT hadm_id FROM work.nicu_adm WHERE in_nicu OR in_scn)""")
cx["org"] = cx.org_name.fillna("").str.strip().str.lower()

n_gram = int((cx.org == "").sum())
print("blood-culture rows read: %d, of which Gram-stain rows (no organism): %d" % (len(cx), n_gram))

CONS = ("coagulase negative", "staphylococcus epidermidis", "staphylococcus hominis",
        "staphylococcus haemolyticus", "staphylococcus capitis", "staphylococcus warneri",
        "staphylococcus lugdunensis", "staphylococcus simulans",
        "staphylococcus saccharolyticus", "staphylococcus pettenkoferi",
        "staphylococcus pasteuri")


def is_cons(o):
    return any(k in o for k in CONS) or o in ("staphylococcus", "staph")


def classify(g):
    """Collapse every row of one episode (admission x collection time x specimen type)."""
    orgs = [o for o in g.org.unique() if o]          # Gram-stain rows drop out here
    if not orgs:
        return pd.Series({"positive": 0, "definite": 0})
    return pd.Series({"positive": 1,
                      "definite": int(any(not is_cons(o) for o in orgs))})


ep = cx.groupby(["hadm_id", "charttime", "spec_type_desc"]).apply(classify).reset_index()
print("episodes after aggregation: %d (positive %d, definite pathogen %d, CoNS only %d)"
      % (len(ep), int(ep.positive.sum()), int(ep.definite.sum()),
         int(((ep.positive == 1) & (ep.definite == 0)).sum())))

# ---------------------------------------------------------------- 2. bilirubin and pairing
bili = q("""SELECT hadm_id, charttime, valuenum AS bili, age_hours FROM work.neo_lab
            WHERE itemid = 50885 AND valuenum > 0 AND valuenum < 45
              AND age_hours >= 0""").dropna().sort_values("charttime")
bt = q("SELECT hadm_id, birth_time FROM work.nicu_adm").set_index("hadm_id").birth_time


def build_pairs(bili_df):
    """One pair per episode with a bilirubin in the 48 h before and the 72 h after."""
    bg = {h: g for h, g in bili_df.groupby("hadm_id")}
    rows = []
    for h, g in ep.groupby("hadm_id"):
        bb = bg.get(h)
        if bb is None:
            continue
        b0 = bt.get(h, pd.NaT)
        for r in g.itertuples():
            pre = bb[(bb.charttime <= r.charttime) &
                     (bb.charttime >= r.charttime - pd.Timedelta(hours=48))]
            post = bb[(bb.charttime >= r.charttime) &
                      (bb.charttime <= r.charttime + pd.Timedelta(hours=72))]
            if pre.empty or post.empty:
                continue
            rows.append(dict(hadm_id=h, t0=r.charttime, positive=int(r.positive),
                             definite=int(r.definite),
                             pre=float(pre.iloc[-1].bili), post=float(post.iloc[0].bili),
                             gap_h=(post.iloc[0].charttime -
                                    pre.iloc[-1].charttime).total_seconds() / 3600,
                             age_d=(r.charttime - b0).total_seconds() / 86400
                             if pd.notna(b0) else np.nan))
    p = pd.DataFrame(rows)
    p["delta"] = p.post - p.pre
    p["cls"] = np.where(p.positive == 0, "neg", np.where(p.definite == 1, "def", "cons"))
    return p


P = build_pairs(bili)

print("paired episodes: %d (positive %d = %.1f%%; definite %d, CoNS %d; %d admissions, %d positive)"
      % (len(P), int(P.positive.sum()), 100 * P.positive.mean(),
         int((P.cls == "def").sum()), int((P.cls == "cons").sum()),
         P.hadm_id.nunique(), P.loc[P.positive == 1, "hadm_id"].nunique()))


def fit(formula, data, label, terms=("positive",)):
    r = smf.ols(formula, data=data).fit(cov_type="cluster",
                                        cov_kwds={"groups": data.hadm_id})
    ci = r.conf_int()
    for t in terms:
        print("  %-46s n=%4d pos=%3d  %+.2f (%+.2f to %+.2f)  p=%.3g"
              % (label, int(r.nobs), int(data.positive.sum()), r.params[t],
                 ci.loc[t, 0], ci.loc[t, 1], r.pvalues[t]))
    return r


print("\nTable S5 - all episodes")
fit("delta ~ positive", P, "unadjusted")
fit("delta ~ positive + gap_h + pre", P, "+ interval and baseline value")
fit("delta ~ positive + gap_h + pre + age_d", P, "+ age at the episode")

# The two pathogen strata are read off one three-level model that shares the negative
# reference, so that the three groups are compared on the same footing.
print("\nTable S5 - pathogen stratum (three-level exposure, negative reference)")
P["cls3"] = pd.Categorical(P.cls, categories=["neg", "def", "cons"])
r3 = smf.ols("delta ~ C(cls3, Treatment('neg')) + gap_h + pre + age_d", data=P).fit(
    cov_type="cluster", cov_kwds={"groups": P.hadm_id})
ci3 = r3.conf_int()
for key, lab in [("def", "definite pathogen (non-CoNS)"),
                 ("cons", "CoNS only")]:
    k = [c for c in r3.params.index if c.endswith("[T.%s]" % key)][0]
    print("  %-46s n=%4d pos=%3d  %+.2f (%+.2f to %+.2f)  p=%.3g"
          % (lab, int(r3.nobs), int((P.cls == key).sum()),
             r3.params[k], ci3.loc[k, 0], ci3.loc[k, 1], r3.pvalues[k]))

# Contrast reported in Table S5: the model-based difference of the two coefficients in the
# three-level fit above, so that the two stratum rows and the contrast share one reference and
# one covariance matrix. Differencing rounded coefficients, or taking the contrast from a
# different (positive-only) model, gives a number that does not reconcile with the rows.
from scipy import stats as _st

kd = [c for c in r3.params.index if c.endswith("[T.def]")][0]
kc = [c for c in r3.params.index if c.endswith("[T.cons]")][0]
diff = r3.params[kd] - r3.params[kc]
V = r3.cov_params()
se_d = np.sqrt(V.loc[kd, kd] + V.loc[kc, kc] - 2 * V.loc[kd, kc])
p_d = 2 * (1 - _st.norm.cdf(abs(diff / se_d)))
print("  %-46s n=%4d        %+.3f (%+.3f to %+.3f)  p=%.3g"
      % ("contrast: definite minus CoNS (three-level)", int(r3.nobs),
         diff, diff - 1.96 * se_d, diff + 1.96 * se_d, p_d))

# Audit only, not tabulated: the same contrast within the 148 culture-positive pairs.
d = P[P.cls.isin(["def", "cons"])].copy()
d["isdef"] = (d.cls == "def").astype(int)
rc = smf.ols("delta ~ isdef + gap_h + pre + age_d", data=d).fit(
    cov_type="cluster", cov_kwds={"groups": d.hadm_id})
cic = rc.conf_int().loc["isdef"]
print("  %-46s n=%4d        %+.3f (%+.3f to %+.3f)  p=%.3g"
      % ("audit: within-positive pairs only", int(rc.nobs),
         rc.params.isdef, cic[0], cic[1], rc.pvalues.isdef))

print("\nTable S5 - baseline-bilirubin tertiles")
P["tert"] = pd.qcut(P.pre, 3, labels=["T1", "T2", "T3"])
for t, g in P.groupby("tert", observed=True):
    r = smf.ols("delta ~ positive", data=g).fit(cov_type="cluster",
                                                cov_kwds={"groups": g.hadm_id})
    ci = r.conf_int().loc["positive"]
    print("  %-46s n=%4d pos=%3d  %+.2f (%+.2f to %+.2f)  p=%.3g   pre %.1f-%.1f mg/dL"
          % (t + " (unadjusted)", len(g), int(g.positive.sum()), r.params.positive,
             ci[0], ci[1], r.pvalues.positive, g.pre.min(), g.pre.max()))
fit("delta ~ positive + C(tert)", P, "adjusted for tertile (fixed effects)")

print("\nSensitivity checks (not tabulated)")
s = build_pairs(bili[bili.age_hours <= 336])
fit("delta ~ positive + gap_h + pre + age_d", s,
    "measurements within the first 14 days only")
big = P.loc[P.positive == 1].reindex(P.loc[P.positive == 1].delta.abs()
                                     .sort_values(ascending=False).index).iloc[0]
inf_adm = big.hadm_id
s = P[P.hadm_id != inf_adm]
fit("delta ~ positive + gap_h + pre + age_d", s,
    "excluding the most influential admission (%d)" % inf_adm)

P.drop(columns=["tert"]).to_csv(os.path.join(OUT, "D_paired_episodes_MIMICIV_clean.csv"),
                                index=False, encoding="utf-8-sig")
P.to_csv(os.path.join(OUT, "D_paired_episodes_MIMICIV_clean_tertile.csv"),
         index=False, encoding="utf-8-sig")
print("\nwritten: D_paired_episodes_MIMICIV_clean.csv in", OUT)
con.close()
