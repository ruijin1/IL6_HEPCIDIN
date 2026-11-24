#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phototherapy as a mediator of the infection-bilirubin association (MIMIC-III).

Companion to the manuscript sentence:

    "Quantitatively, phototherapy accounted for only 7.2% of the association in a
     measurement-level model of the MIMIC-III cohort covering the full 14-day window,
     and excluding measurements taken while phototherapy was running changed the
     estimate by 1.6% ..."

Design
------
Measurement-level model of the MIMIC-III NICU cohort (age <= 28 days at ICU entry):

    log(total bilirubin) ~ infected + age_at_measurement        (i)  total effect
    log(total bilirubin) ~ infected + age + in_phototherapy     (iii) direct effect
    in_phototherapy      ~ infected + age                       (a path)

The share of the association carried by the phototherapy route is the product of
coefficients, a * b, divided by the total effect. Because phototherapy is started
*because* bilirubin is high, this quantity is a descriptive decomposition of the
observed association, not an estimate of a treatment effect; the manuscript labels
the analysis as post hoc for that reason.

Phototherapy is identified from chartevents itemids 3569 (blanket) and 3570 (light).
A measurement is counted as "during phototherapy" if its charttime falls between the
first and the last phototherapy charting of the same admission. That interval is a
deliberate over-approximation: phototherapy is charted repeatedly while running, and
using the first-to-last span avoids under-counting gaps in charting.

Outputs
-------
phototherapy_mediation.csv   one row per model plus the decomposition
phototherapy_mediation.log   the same numbers as plain text

Environment
-----------
PGHOST, PGPORT, PGUSER, PGPASSWORD, PGDATABASE   PostgreSQL holding the MIMIC-III schema
OUTDIR                                           output directory (default: .)
"""

import os

import numpy as np
import pandas as pd
import psycopg2
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get("OUTDIR", HERE)
SCHEMA = os.environ.get("MIMIC3_SCHEMA", "mimiciii")

NEGATIVE_PATTERN = (
    r"NO GROWTH|^NEGATIVE|NEGATIVE FOR|NO BACTERIAL|NO ANAEROBIC|"
    r"PENDING|CANCELLED|NORMAL FLORA"
)
# "COAGULASE NEGATIVE" must NOT be treated as a negative culture.
PHOTOTHERAPY_ITEMS = ("3569", "3570")
BILIRUBIN_ITEM = "50885"
BIRTH_WEIGHT_ITEM = "3723"


def q(sql, con=None):
    con = con or psycopg2.connect(
        host=os.environ.get("PGHOST", "localhost"),
        port=os.environ.get("PGPORT", "5432"),
        user=os.environ.get("PGUSER", "postgres"),
        password=os.environ.get("PGPASSWORD", ""),
        dbname=os.environ.get("PGDATABASE", "mimiciii"),
    )
    return pd.read_sql(sql, con)


def main():
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(str(s))

    log("Phototherapy mediation analysis (MIMIC-III NICU, age <= 28 days)")
    log("=" * 78)

    stay = q(f"""
        select i.icustay_id, i.hadm_id, i.subject_id, i.intime, i.outtime,
               p.dob, a.hospital_expire_flag
        from {SCHEMA}.icustays i
        join {SCHEMA}.patients p using (subject_id)
        join {SCHEMA}.admissions a using (hadm_id)
        where i.first_careunit = 'NICU'
          and (i.intime::date - p.dob::date) between 0 and 28
    """)
    log("ICU stays: %d   admissions: %d" % (len(stay), stay.hadm_id.nunique()))

    mic = q(f"""
        select hadm_id, charttime, spec_type_desc, org_name
        from {SCHEMA}.microbiologyevents
        where spec_type_desc ~* 'BLOOD'
    """)
    mic = mic[mic.hadm_id.isin(stay.hadm_id)].copy()
    mic["org"] = mic.org_name.fillna("").str.upper().str.strip()
    mic["is_negative"] = (mic.org == "") | mic.org.str.contains(
        NEGATIVE_PATTERN, na=False, regex=True
    )
    events = (
        mic.groupby(["hadm_id", "charttime", "spec_type_desc"])
        .agg(positive=("is_negative", lambda s: int((~s).any())))
        .reset_index()
    )
    infected = set(events.loc[events.positive == 1, "hadm_id"])
    tested = set(events.hadm_id)
    negative_only = tested - infected
    log(
        "Blood-culture events: %d   admissions tested: %d   positive: %d   "
        "tested but never positive: %d"
        % (len(events), len(tested), len(infected), len(negative_only))
    )

    lab = q(f"""
        select hadm_id, charttime, valuenum
        from {SCHEMA}.labevents
        where itemid = {BILIRUBIN_ITEM}
          and valuenum is not null and valuenum between 0.1 and 40
    """)
    lab = lab[lab.hadm_id.isin(stay.hadm_id)].copy()
    dob = stay.drop_duplicates("hadm_id").set_index("hadm_id")["dob"]
    lab["age_d"] = (lab.charttime - lab.hadm_id.map(dob)).dt.total_seconds() / 86400.0
    lab = lab[(lab.age_d >= 0) & (lab.age_d <= 14)]
    log(
        "Total-bilirubin measurements (days 0-14): %d   admissions: %d   median: %.2f mg/dL"
        % (len(lab), lab.hadm_id.nunique(), lab.valuenum.median())
    )

    ph = q(f"""
        select icustay_id, charttime
        from {SCHEMA}.chartevents
        where itemid in ({", ".join(PHOTOTHERAPY_ITEMS)})
    """)
    ph = ph[ph.icustay_id.isin(stay.icustay_id)]
    span = ph.groupby("icustay_id").charttime.agg(["min", "max"]).reset_index()
    span.columns = ["icustay_id", "ph_start", "ph_end"]
    span["hadm_id"] = span.icustay_id.map(stay.set_index("icustay_id")["hadm_id"])
    log(
        "Phototherapy charting: %d records   stays covered: %d (%.1f%%)"
        % (len(ph), span.icustay_id.nunique(), 100.0 * span.icustay_id.nunique() / len(stay))
    )

    lab = lab.merge(span[["hadm_id", "ph_start", "ph_end"]], on="hadm_id", how="left")
    lab["in_ph"] = (
        (lab.charttime >= lab.ph_start) & (lab.charttime <= lab.ph_end)
    ).fillna(False).astype(int)
    log(
        "Measurements taken while phototherapy was running: %d (%.1f%%)"
        % (int(lab.in_ph.sum()), 100.0 * lab.in_ph.mean())
    )

    d = lab[lab.hadm_id.isin(infected | negative_only)].copy()
    d["infected"] = d.hadm_id.isin(infected).astype(int)
    d["log_bili"] = np.log(d.valuenum)
    log(
        "Analysis set: %d measurements from %d admissions (%d infected, %d tested negative)"
        % (
            len(d),
            d.hadm_id.nunique(),
            d.loc[d.infected == 1, "hadm_id"].nunique(),
            d.loc[d.infected == 0, "hadm_id"].nunique(),
        )
    )

    rows = []

    def fit(formula, data, label, keep=None):
        m = smf.ols(formula, data=data).fit(
            cov_type="cluster", cov_kwds={"groups": data["hadm_id"]}
        )
        ci = m.conf_int().loc["infected"]
        rows.append(
            dict(
                model=label,
                term="infected",
                coefficient=m.params["infected"],
                ci_low=ci[0],
                ci_high=ci[1],
                p=m.pvalues["infected"],
                n_measurements=len(data),
                n_admissions=data.hadm_id.nunique(),
            )
        )
        log(
            "  %-42s %+.4f (95%% CI %+.4f to %+.4f)  p=%.2e  n=%d"
            % (label, m.params["infected"], ci[0], ci[1], m.pvalues["infected"], len(data))
        )
        return m

    base = "log_bili ~ infected + age_d"
    log("")
    log("Models (dependent variable: log total bilirubin; SE clustered by admission)")
    total = fit(base, d, "(i) all measurements, not adjusted for phototherapy").params["infected"]
    fit(base, d[d.in_ph == 0], "(ii) measurements outside phototherapy only")
    direct_model = fit(base + " + in_ph", d, "(iii) adjusted for measurement during phototherapy")
    direct = direct_model.params["infected"]
    b_path = direct_model.params["in_ph"]
    a_path = smf.ols("in_ph ~ infected + age_d", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": d["hadm_id"]}
    ).params["infected"]

    indirect = a_path * b_path
    log("")
    log("Decomposition of the phototherapy route")
    log("  a (infection -> measurement during phototherapy): %+.4f" % a_path)
    log("  b (measurement during phototherapy -> log bilirubin): %+.4f" % b_path)
    log("  a*b (indirect): %+.4f" % indirect)
    log("  share of the total association: %.1f%%" % (100.0 * indirect / total))
    log("  change in the estimate when phototherapy measurements are excluded: %+.1f%%"
        % (100.0 * (rows[1]["coefficient"] - total) / abs(total)))

    rows.append(
        dict(
            model="decomposition",
            term="a_path (infection -> phototherapy)",
            coefficient=a_path,
            ci_low=np.nan,
            ci_high=np.nan,
            p=np.nan,
            n_measurements=len(d),
            n_admissions=d.hadm_id.nunique(),
        )
    )
    rows.append(
        dict(
            model="decomposition",
            term="b_path (phototherapy -> log bilirubin)",
            coefficient=b_path,
            ci_low=np.nan,
            ci_high=np.nan,
            p=np.nan,
            n_measurements=len(d),
            n_admissions=d.hadm_id.nunique(),
        )
    )
    rows.append(
        dict(
            model="decomposition",
            term="indirect (a*b)",
            coefficient=indirect,
            ci_low=np.nan,
            ci_high=np.nan,
            p=np.nan,
            n_measurements=len(d),
            n_admissions=d.hadm_id.nunique(),
        )
    )
    rows.append(
        dict(
            model="decomposition",
            term="share of total association (%)",
            coefficient=100.0 * indirect / total,
            ci_low=np.nan,
            ci_high=np.nan,
            p=np.nan,
            n_measurements=len(d),
            n_admissions=d.hadm_id.nunique(),
        )
    )

    out = pd.DataFrame(rows)
    os.makedirs(OUTDIR, exist_ok=True)
    out.to_csv(os.path.join(OUTDIR, "phototherapy_mediation.csv"), index=False)
    with open(os.path.join(OUTDIR, "phototherapy_mediation.log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    log("")
    log("Written: phototherapy_mediation.csv, phototherapy_mediation.log")


if __name__ == "__main__":
    main()
