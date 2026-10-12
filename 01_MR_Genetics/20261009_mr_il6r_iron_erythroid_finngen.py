# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""* The last piece of the mechanism layer: IL-6 pathway -> iron metabolism / erythroid indices (hepcidin axis)

Scientific question:
  The preferred explanation of 12.5(1) / 14.4a is -- **IL-6 signaling up -> hepatic HAMP up -> hepcidin restricts iron utilization ->
  erythropoiesis substrates limited -> heme turnover down -> bilirubin production down**.
  This step occurs at the erythroid/systemic level and cannot be measured at all in hepatocyte monoculture (13.7 limitation 3).
  This section tests it directly **in humans** using **genetics**.

Data (FinnGen R12 * Kanta laboratory quantitative traits, all inverse-rank-normalized -> SD units):
  * All outcomes come from `lab_values/summary_stats/finngen_R12_<OMOPID>.gz`, **with .tbi index** -> fetch only regions, no full-file download
  * Columns of that file: `#chrom pos ref alt pval mlogp beta sebeta af_alt`

Design:
  * Main analysis   IL6R cis-pQTL (4 independent pQTL datasets) -> each iron/erythroid index (Wald ratio)
  * Positive control **TMPRSS6 rs855791**: the well-established A allele -> hepcidin down -> iron utilization up -> MCH/Hb up
             (if this pipeline cannot reproduce that known direction, the allele harmonization or data retrieval is wrong, and all results are void)
  * Locus scan: regional min p in HAMP (chr19:35.1-35.4M) and HFE (chr6:26.05-26.15M)
  * Reverse control: UGT1A1 region -> total bilirubin (a huge peak should appear)
"""
import os, sys, time, json
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\workbuddy_workspace\02-scripts")
import numpy as np
import pandas as pd
from scipy import stats
import finngen_tabix as ft

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
CACHE = os.path.join(DIR, "tabix_cache")
BASE = "https://storage.googleapis.com/finngen-public-data-r12/lab_values/summary_stats/"
os.makedirs(CACHE, exist_ok=True)

# ---- Outcomes: iron / erythroid / bilirubin (OMOPIDs from Kanta_labs_GWAS_results_v2_summary) ----
TRAITS = [
    ("3000963", "Hemoglobin", "erythroid"),
    ("3024731", "MCV", "erythroid"),
    ("3035941", "MCH", "erythroid"),
    ("3026361", "Red blood cell count", "erythroid"),
    ("3009542", "Hematocrit Hct", "erythroid"),
    ("3002385", "RDW", "erythroid"),
    ("3027945", "Reticulocyte %", "erythroid"),
    ("3023520", "Reticulocyte count", "erythroid"),
    ("3024172", "Erythropoietin EPO", "erythroid"),
    ("3001122", "Ferritin", "iron"),
    ("3022948", "Serum iron", "iron"),
    ("3004789", "Transferrin", "iron"),
    ("40763571", "Transferrin saturation TSAT", "iron"),
    ("3015399", "Soluble transferrin receptor sTfR", "iron"),
    ("3006140", "Total bilirubin (reference)", "bili"),
]

REGIONS = {
    "IL6R":   ("1", 154_300_000, 154_800_000),
    "HAMP":   ("19", 35_100_000, 35_400_000),
    "HFE":    ("6",  26_050_000, 26_150_000),
    "TMPRSS6": ("22", 36_950_000, 37_250_000),
    "UGT1A1": ("2", 233_600_000, 233_900_000),
}


def rows_to_df(rows):
    if not rows:
        return pd.DataFrame(columns=["chrom", "pos", "ref", "alt", "pval", "mlogp",
                                     "beta", "sebeta", "af_alt"])
    d = pd.DataFrame([r.split("\t") for r in rows],
                     columns=["chrom", "pos", "ref", "alt", "pval", "mlogp",
                              "beta", "sebeta", "af_alt"][:len(rows[0].split("\t"))])
    for c in ["pos", "pval", "mlogp", "beta", "sebeta", "af_alt"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def wald(b_ex, se_ex, b_out, se_out, ea_ex, nea_ex, ea_out, nea_out):
    if ea_ex == ea_out and nea_ex == nea_out:
        s = 1.0
    elif ea_ex == nea_out and nea_ex == ea_out:
        s = -1.0
    else:
        return None
    bo = s * b_out
    if b_ex == 0:
        return None
    ratio = bo / b_ex
    se = abs(ratio) * np.sqrt((se_out / bo) ** 2 + (se_ex / b_ex) ** 2) if bo != 0 else np.nan
    return dict(ratio=ratio, se=se, lo=ratio - 1.96 * se, hi=ratio + 1.96 * se,
                p=2 * stats.norm.sf(abs(ratio / se)) if se and se > 0 else np.nan,
                F=(b_ex / se_ex) ** 2, align=s)


def load_g38_positions():
    """* Key: the OpenGWAS position is in **GRCh37**, while FinnGen R12 is in **GRCh38** (the IL6R region differs by ~27.5 kb)
       -> never match iv["position"] directly against FinnGen by position (everything would miss).
       Here we resolve the GRCh38 coordinates via the Ensembl REST API and cache them (including the two control SNPs with known directions).
    """
    fp = os.path.join(DIR, "_il6r_snp_coords.json")
    if os.path.exists(fp):
        return json.load(open(fp, encoding="utf-8"))
    import urllib.request
    def eg(rs):
        r = urllib.request.Request("https://rest.ensembl.org/variation/human/%s?content-type=application/json" % rs,
                                   headers={"User-Agent": "Mozilla/5.0"})
        d = json.loads(urllib.request.urlopen(r, timeout=60).read().decode())
        for m in d.get("mappings", []):
            if m.get("assembly_name") == "GRCh38":
                return dict(chrom=m["seq_region_name"], pos=m["start"], alleles=m.get("allele_string"))
        return None
    out = {}
    for rs in ["rs4129267", "rs4601580", "rs7514452", "rs6690230",
               "rs855791", "rs1799945", "rs1800562", "rs10421768"]:
        try:
            out[rs] = eg(rs)
            print("  %-12s GRCh38 %s:%s  %s" % (rs, out[rs]["chrom"], out[rs]["pos"], out[rs]["alleles"]))
        except Exception as e:
            print("  %-12s resolution failed %s" % (rs, str(e)[:60]))
            out[rs] = None
        time.sleep(0.4)
    json.dump(out, open(fp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return out


def find_by_pos(d, pos, tol=0):
    if not len(d):
        return None
    hit = d[d["pos"] == pos]
    if len(hit):
        return hit.iloc[0]
    for k in range(1, tol + 1):
        hit = d[(d["pos"] >= pos - k) & (d["pos"] <= pos + k)]
        if len(hit):
            return hit.iloc[0]
    return None


def main():
    inst = json.load(open(os.path.join(DIR, "MR_il6_saa_instruments.json"), encoding="utf-8"))["instruments"]
    il6r = [iv for k, v in inst.items() if k.startswith("IL6R") for iv in v]
    G38 = load_g38_positions()
    print("IL6R cis instruments (%d): %s" % (len(il6r), [iv["rsid"] for iv in il6r]))
    for iv in il6r:
        g = G38.get(iv["rsid"]) or {}
        iv["g38_chrom"], iv["g38_pos"] = g.get("chrom"), g.get("pos")
        print("   %-12s GRCh37 pos=%s -> GRCh38 pos=%s" % (iv["rsid"], iv["position"], iv["g38_pos"]))

    results, locus, dirchk = [], [], []
    t_start = time.time()
    for omop, name, grp in TRAITS:
        url = BASE + "finngen_R12_%s.gz" % omop
        tbi_fp = os.path.join(CACHE, "finngen_%s.tbi" % omop)
        try:
            tbi = ft.read_tbi(url + ".tbi", tbi_fp)
        except Exception as e:
            print("\n!! %s(%s) tbi failed: %s" % (omop, name, str(e)[:70])); continue
        print("\n" + "=" * 104)
        print("%s  (%s, OMOPID %s)   elapsed %.0f s" % (name, grp, omop, time.time() - t_start))
        print("=" * 104)

        # --- Fetch data for each region ---
        reg = {}
        for rn, (ch, lo, hi) in REGIONS.items():
            reg[rn] = rows_to_df(ft.fetch_region(url, tbi, ch, lo, hi))
        for rn in REGIONS:
            d = reg[rn]
            if len(d):
                print("  region %-8s %5d variants, min p = %.3g" % (rn, len(d), d.pval.min()))
            else:
                print("  region %-8s no data" % rn)
            if len(d) and rn in ("HAMP", "HFE", "UGT1A1"):
                top = d.nsmallest(2, "pval")
                for _, r in top.iterrows():
                    locus.append(dict(trait=name, omop=omop, region=rn, pos=int(r["pos"]),
                                      ref=r["ref"], alt=r["alt"], beta=r["beta"],
                                      se=r["sebeta"], p=r["pval"], af=r["af_alt"]))

        # --- Main analysis: IL6R cis-pQTL -> this index (matched by GRCh38 position) ---
        d = reg["IL6R"]
        if not len(d):
            continue
        print("\n  -- IL6R cis-pQTL -> %s --" % name)
        for iv in il6r:
            if iv.get("g38_pos") is None:
                continue
            r = find_by_pos(d, int(iv["g38_pos"]), tol=2)
            if r is None:
                print("    %-26s %-12s GRCh38 pos %s is missing in FinnGen" % (iv["trait"][:26], iv["rsid"], iv["g38_pos"]))
                continue
            res = wald(iv["beta"], iv["se"], r["beta"], r["sebeta"],
                       iv["ea"], iv["nea"], r["alt"], r["ref"])
            if res is None:
                print("    %-26s %-12s alleles cannot be harmonized (exposure %s>%s, outcome %s>%s)"
                      % (iv["trait"][:26], iv["rsid"], iv["ea"], iv["nea"], r["ref"], r["alt"])); continue
            print("    %-24s %-12s outcome beta=%+.4f(p=%.3g)  Wald %+.4f (95%%CI %+.4f~%+.4f) p=%.3g  F=%.0f"
                  % (iv["trait"][:24], iv["rsid"], r["beta"], r["pval"], res["ratio"],
                     res["lo"], res["hi"], res["p"], res["F"]))
            results.append(dict(trait=name, group=grp, omop=omop, exposure=iv["trait"], rsid=iv["rsid"],
                                b_ex=iv["beta"], se_ex=iv["se"], F=res["F"],
                                b_out=r["beta"], se_out=r["sebeta"], p_out=r["pval"],
                                ratio=res["ratio"], lo=res["lo"], hi=res["hi"], p_mr=res["p"]))

        # --- * Direction controls (iron-axis loci with known literature directions): only check whether the sign matches the literature ---
        # TMPRSS6 rs855791-A -> iron utilization up -> MCH/Hb up; HFE rs1799945 (H63D) -> iron up -> MCH/Hb up
        DIRCHK = [("rs855791", "TMPRSS6", "A"), ("rs1799945", "HFE", "G"), ("rs1800562", "HFE", "A")]
        for rs, ln, expect_alt in DIRCHK:
            g = G38.get(rs)
            if not g:
                continue
            regname = "TMPRSS6" if ln == "TMPRSS6" else "HFE"
            dd = reg[regname]
            rr = find_by_pos(dd, int(g["pos"]), tol=2)
            if rr is None:
                print("    [direction control] %-10s(%s) position missing" % (rs, ln)); continue
            eff = rr["alt"] if rr["alt"] in str(g["alleles"]) else rr["ref"]
            print("    [direction control] %-10s(%s) chr%s:%d %s>%s beta=%+.4f p=%.3g  literature expectation: effective allele -> target index up"
                  % (rs, ln, g["chrom"], int(rr["pos"]), rr["ref"], rr["alt"], rr["beta"], rr["pval"]))
            # also persist the direction controls: the manuscript quotes them, so they must be traceable
            z = abs(float(rr["beta"]) / float(rr["sebeta"])) if float(rr["sebeta"]) else float("nan")
            dirchk.append(dict(trait=name, omop=omop, region=regname, rsid=rs, chrom=g["chrom"],
                               pos=int(rr["pos"]), ref=rr["ref"], alt=rr["alt"], eaf=rr["af_alt"],
                               beta=rr["beta"], se=rr["sebeta"], z=z, pval=rr["pval"],
                               literature_effect_allele=expect_alt))

    if results:
        p = os.path.join(DIR, "MR_il6r_iron_erythroid_results.csv")
        pd.DataFrame(results).to_csv(p, index=False, encoding="utf-8-sig")
        print("\nMain results saved: %s (%d rows)" % (p, len(results)))
    if locus:
        p = os.path.join(DIR, "finngen_lab_locus_scan.csv")
        pd.DataFrame(locus).to_csv(p, index=False, encoding="utf-8-sig")
        print("Locus scan saved: %s (%d rows)" % (p, len(locus)))
    if dirchk:
        p = os.path.join(DIR, "finngen_iron_direction_controls.csv")
        pd.DataFrame(dirchk).to_csv(p, index=False, encoding="utf-8-sig")
        print("Direction controls saved: %s (%d rows)" % (p, len(dirchk)))


if __name__ == "__main__":
    main()
