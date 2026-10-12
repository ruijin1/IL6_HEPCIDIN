# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""MIMIC-III NICU: extract the "infection x erythrocyte index/reticulocyte" observational data

Why this step (the gap of the triangulation):
  Section 15.6 used FinnGen genetics to show IL6R -> iron utilisation down + iron-restricted erythropoiesis (MCV/MCH down, RDW up, sTfR up).
  But the **clinical side is completely missing** -- the observational MIMIC work only did "infection -> bilirubin down" (-0.355).
  This section fills in **the clinical side of the same mechanism**: in NICU neonates, is infection accompanied by
  **MCV down / MCH down / RDW up / reticulocytes down** -- i.e. the textbook combination of iron-restricted erythropoiesis.

* Availability explored first (this step decides whether the analysis can hold up; it cannot be skipped):
  - Direct iron-metabolism assays in neonates are **almost zero**: Iron 5 rows, Ferritin 7, TIBC 4, Transferrin 4, CRP 23
    -> one **cannot** directly measure "serum iron decline"; one can only detour through **erythrocyte indices + reticulocytes**
  - Available: MCV 51250 (10,270 rows), MCH 51248 (10,270), MCHC 51249 (10,275),
    RDW 51277 (10,259), Hb 51222 (10,318), Hct 51221 (13,301),
    reticulocytes 51283 automated (980) + 51284 manual (958), total bilirubin 50885 (24,192)
  - NICU cohort: `icustays.first_careunit='NICU'` -> **8,100** stays (consistent with the previous BPD work)

Output (03-data\MIMIC_NICU_iron\):
  nicu_cbc.csv       erythrocyte-index time series in NICU neonates
  nicu_culture.csv   blood cultures (positive/negative after aggregation by time x site)
  nicu_weight.csv    birth weight / weight at ICU entry (used for adjustment)
  nicu_base.csv      cohort baseline (DOB, age, gestational-age related)
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd
import psycopg2

OUT = r"D:\workbuddy_workspace\03-data\MIMIC_NICU_iron"
os.makedirs(OUT, exist_ok=True)
CONN = dict(host="localhost", port=5432, user="postgres", password="1234", dbname="mimiciii")

# Index itemid -> variable name (with unit-conversion notes)
CBC = {
    51250: "mcv",        # fL
    51248: "mch",        # pg
    51249: "mchc",       # g/dL
    51277: "rdw",        # %
    51222: "hb_heme",    # g/dL  (Hematology)
    50811: "hb_bg",      # g/dL  (Blood Gas, alternative)
    51221: "hct",        # %
    51283: "retic_auto", # %  automated
    51284: "retic_man",  # %  manual
    50885: "bili_total", # mg/dL (bilirubin, as the control endpoint of the previous section)
    50889: "crp",        # mg/dL (very few)
    50952: "iron",
    50924: "ferritin",
    50953: "tibc",
    50998: "transferrin",
}


def main():
    c = psycopg2.connect(**CONN)
    print("=" * 104)
    print("(1) NICU cohort (icustays.first_careunit='NICU')")
    print("=" * 104)
    base = pd.read_sql("""
        select i.subject_id, i.hadm_id, i.icustay_id, i.intime, i.outtime, i.los,
               p.dob, p.gender, a.admittime, a.dischtime, a.admission_type, a.hospital_expire_flag,
               extract(epoch from (i.intime - p.dob))/86400.0 as age_days_at_icu
        from mimiciii.icustays i
        join mimiciii.patients p on i.subject_id=p.subject_id
        join mimiciii.admissions a on i.hadm_id=a.hadm_id
        where i.first_careunit='NICU'
    """, c)
    # The same hadm with multiple NICU stays -> take the first
    base = base.sort_values("intime").drop_duplicates("hadm_id", keep="first")
    print("   NICU hadm count = %d (unique subject %d)" % (len(base), base.subject_id.nunique()))
    print("   age at ICU entry (days): median %.1f, IQR %.1f-%.1f, min %.2f, max %.1f"
          % (base.age_days_at_icu.median(), base.age_days_at_icu.quantile(.25),
             base.age_days_at_icu.quantile(.75), base.age_days_at_icu.min(), base.age_days_at_icu.max()))
    base.to_csv(os.path.join(OUT, "nicu_base.csv"), index=False, encoding="utf-8-sig")

    ids = tuple(base.hadm_id.astype(int).tolist())
    print("\n" + "=" * 104)
    print("(2) Erythrocyte index / iron-metabolism assay time series")
    print("=" * 104)
    ids_sql = ",".join(str(x) for x in ids)
    items = ",".join(str(k) for k in CBC)
    lab = pd.read_sql("""
        select le.hadm_id, le.itemid, le.charttime, le.value, le.valuenum, le.valueuom, le.flag
        from mimiciii.labevents le
        where le.hadm_id in (%s) and le.itemid in (%s)
    """ % (ids_sql, items), c)
    lab["var"] = lab.itemid.map(CBC)
    lab["valuenum"] = pd.to_numeric(lab.valuenum, errors="coerce")
    lab = lab.dropna(subset=["valuenum", "charttime"])
    print("   total rows %d" % len(lab))
    print("   %-12s %8s %8s   %s" % ("Index", "rows", "non-null", "hadm covered"))
    for v, g in lab.groupby("var"):
        print("   %-12s %8d %8d   %d" % (v, len(g), g.valuenum.notna().sum(), g.hadm_id.nunique()))
    lab.to_csv(os.path.join(OUT, "nicu_cbc.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 104)
    print("(3) Blood cultures (positive/negative determined after aggregation by hadm x collection time x site)")
    print("=" * 104)
    mb = pd.read_sql("""
        select m.hadm_id, m.charttime, m.spec_type_desc, m.org_name, m.ab_name
        from mimiciii.microbiologyevents m
        where m.hadm_id in (%s) and upper(m.spec_type_desc) like '%%BLOOD%%'
          and m.charttime is not null
    """ % ids_sql, c)
    print("   raw rows %d (hadm %d)" % (len(mb), mb.hadm_id.nunique()))
    agg = (mb.groupby(["hadm_id", "charttime", "spec_type_desc"])
             .agg(n_org=("org_name", lambda s: s.notna().sum()),
                  orgs=("org_name", lambda s: ",".join(sorted({str(x) for x in s.dropna()}))),
                  n_ab=("ab_name", lambda s: s.notna().sum()))
             .reset_index())
    agg["positive"] = (agg.n_org > 0).astype(int)
    print("   culture episodes after aggregation %d" % len(agg))
    print("   positive %d episodes / negative %d episodes; positive involves %d hadm"
          % (agg.positive.sum(), (agg.positive == 0).sum(), agg[agg.positive == 1].hadm_id.nunique()))
    top = (agg[agg.positive == 1].orgs.str.split(",").explode().str.strip()
           .value_counts().head(12))
    print("\n   Most common positive organisms:")
    for k, v in top.items():
        if str(k) and str(k) != "":
            print("      %-42s %d" % (str(k)[:42], v))
    agg.to_csv(os.path.join(OUT, "nicu_culture.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 104)
    print("(4) Weight (used for adjustment; MIMIC-III neonates have no numerical gestational age, so weight is the only proxy)")
    print("=" * 104)
    w = pd.read_sql("""
        select ce.hadm_id, ce.charttime, ce.itemid, ce.valuenum, di.label
        from mimiciii.chartevents ce
        join mimiciii.d_items di on ce.itemid=di.itemid
        where ce.hadm_id in (%s) and ce.valuenum is not null
          and (lower(di.label) like '%%weight%%' or lower(di.label) like '%%birth%%')
    """ % ids_sql, c)
    print("   weight-related items %d rows; labels:" % len(w))
    for l, n in w.label.value_counts().head(12).items():
        print("      %-46s %d" % (str(l)[:46], n))
    w.to_csv(os.path.join(OUT, "nicu_weight.csv"), index=False, encoding="utf-8-sig")
    print("\nDone -> %s" % OUT)


if __name__ == "__main__":
    main()
