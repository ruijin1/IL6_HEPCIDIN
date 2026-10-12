# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""IL6R / UGT1A1 / HAMP / TMPRSS6 / HFE -> **FinnGen R12 neonatal phenotypes** (full-file version)

Why re-run:
  The previous round used a **partially downloaded** jaundice file (truncated at chr9), so **HAMP on chr19 and TMPRSS6 on chr22 were not read at all**.
  Now all 5 files are fully downloaded (~750 MB each), so two key pieces of evidence can be added:
    * HAMP (hepcidin) locus -> neonatal jaundice -- an **independent** test of the "hepcidin/iron-erythroid axis" that does not rely on IL6R instruments
    * TMPRSS6 rs855791 -> neonatal jaundice -- a **direction control** with a known direction (A allele -> hepcidin down -> iron up -> heme turnover up -> bilirubin up)

Also added: **IL6R -> neonatal bacterial sepsis** (ncase 203) -- this is the genetic test that best matches the MIMIC observational side
(culture-positive -> bilirubin -0.355) **in terms of phenotype**.

Direction logic (read it backwards and every conclusion flips):
  IL6R cis-pQTL = each SD higher **sIL6R protein** = less membrane-bound IL6R = **weaker classical IL-6 signaling**
  (positive control in this project: Wald ratio for IL6R -> CRP is -0.115, i.e. sIL6R up -> CRP down)
  => If "IL-6 up -> HAMP up -> iron restriction -> heme turnover down -> bilirubin down" holds,
     then **sIL6R up should -> bilirubin up, jaundice risk up** (log-OR > 0)
  => If "IL-6 up -> hepatic UGT1A1 down -> bilirubin up" holds on its own, the direction is also **sIL6R up -> jaundice up** (same direction)
  => The two models make the **same predicted direction** for neonatal jaundice (both positive), so direction cannot distinguish the mechanisms;
     it can only test "whether the IL-6 pathway is genetically associated with neonatal jaundice".

Data: GRCh38; columns #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt
"""
import os, sys, gzip, json, math
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"

PHENO = [
    ("P16_NEONTAL_JAUND_OTH_UNSP_CAUSES", "neonatal jaundice", 499, 499688),
    ("P16_BACTERIAL_SEPSIS_NEWBO", "neonatal bacterial sepsis", 203, 499933),
    ("P16_RESPI_DISTR_NEWBO", "neonatal respiratory distress", 498, 499476),
    ("P16_CHRONIC_RESPI_DISEA_ORIGINA_PERINA_PERIOD", "chronic respiratory disease of perinatal origin (BPD)", 69, 499476),
    ("P16_DISORD_RELATED_GESTATION_LOW_BIRTHWGHTT_NECIFIED", "preterm birth/low birth weight", 782, 499424),
]

# GRCh38 regions
REGIONS = [
    ("IL6R",    "1",  154_300_000, 154_800_000),
    ("UGT1A1",  "2",  233_600_000, 233_900_000),
    ("HFE",     "6",   26_080_000,  26_105_000),
    ("HAMP",    "19",  35_050_000,  35_450_000),
    ("TMPRSS6", "22",  36_950_000,  37_200_000),
]

WANT_RS = {
    # IL6R cis-pQTL
    "rs4601580", "rs7514452", "rs6690230", "rs4129267",
    # UGT1A1 locus (positive control)
    "rs6754859", "rs887829", "rs4148323", "rs6742078", "rs12995326",
    # iron-metabolism direction control
    "rs855791",   # TMPRSS6
    "rs1800562",  # HFE C282Y
    # negative control
    "rs11024589", "rs10430893", "rs3116654",
}

NEEDED_CHROMS = {ch for _, ch, _, _ in REGIONS}


def scan(fp, regions, want_rs, verbose=True):
    """Single streaming scan. The file is sorted by chromosome + position -> we can stop early once the highest target chromosome has been passed."""
    cols = None
    reg = {n: [] for n, _, _, _ in regions}
    byrs = {}
    last = ("?", 0)
    max_needed = max(int(c) for c in NEEDED_CHROMS if c.isdigit())
    n_line = 0
    with gzip.open(fp, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#chrom") or line.startswith("chrom"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                idx = {c: i for i, c in enumerate(cols)}
                continue
            fld = line.rstrip("\n").split("\t")
            if len(fld) < 8:
                continue
            n_line += 1
            ch = fld[idx["chrom"]]
            # Early stop: the highest target chromosome has been passed
            if ch.isdigit():
                if int(ch) > max_needed:
                    break
            elif ch in ("X", "Y", "MT"):
                break
            try:
                ipos = int(fld[idx["pos"]])
            except ValueError:
                continue
            last = (ch, ipos)
            rec = dict(chrom=ch, pos=ipos, ref=fld[idx["ref"]], alt=fld[idx["alt"]],
                       rsids=fld[idx["rsids"]], nearest=fld[idx.get("nearest_genes", 5)],
                       p=float(fld[idx["pval"]]),
                       beta=float(fld[idx["beta"]]) if fld[idx["beta"]] not in ("NA", "") else np.nan,
                       se=float(fld[idx["sebeta"]]) if fld[idx["sebeta"]] not in ("NA", "") else np.nan,
                       af=float(fld[idx["af_alt"]]) if fld[idx["af_alt"]] not in ("NA", "") else np.nan)
            for nm, rch, lo, hi in regions:
                if ch == rch and lo <= ipos <= hi:
                    reg[nm].append(rec)
            for rs in str(rec["rsids"]).split(","):
                rs = rs.strip()
                if rs in want_rs:
                    byrs.setdefault(rs, []).append(rec)
    if verbose:
        print("    Scanned %d lines; endpoint chr%s:%s (stopped early)" % (n_line, last[0], last[1]))
    return cols, reg, byrs


def align(b_ex, se_ex, b_out, se_out, ea_ex, nea_ex, ea_out, nea_out):
    """Allele harmonization + Wald ratio. Returns a dict or None."""
    if ea_ex == ea_out and nea_ex == nea_out:
        s = 1
    elif ea_ex == nea_out and nea_ex == ea_out:
        s = -1
    else:
        return None
    bo = s * b_out
    ratio = bo / b_ex
    se = abs(ratio) * math.sqrt((se_out / bo) ** 2 + (se_ex / b_ex) ** 2) if bo != 0 else np.nan
    p = 2 * stats.norm.sf(abs(ratio / se)) if se and se > 0 else np.nan
    return dict(ratio=ratio, se=se, lo=ratio - 1.96 * se, hi=ratio + 1.96 * se, p=p,
                F=(b_ex / se_ex) ** 2, align=s, b_out_aligned=bo)


def expected_logor(b_sd_per_allele, ncase, ncontrol):
    """Threshold model: convert the "per-allele SD effect" into an expected log-OR (K is very small -> i = φ(z)/K)."""
    K = ncase / (ncase + ncontrol)
    z = stats.norm.isf(K)
    i = stats.norm.pdf(z) / K
    return b_sd_per_allele * i, K, i


def power_at(exp_lor, se, alpha=0.05):
    """Given the expected log-OR and SE, compute the two-sided power at the 5% level."""
    if not (se and se > 0) or not np.isfinite(exp_lor):
        return np.nan
    lam = abs(exp_lor) / se
    return float(stats.norm.cdf(lam - 1.96) + stats.norm.cdf(-lam - 1.96))


def mde(se, power=0.80, alpha=0.05):
    return (1.96 + stats.norm.isf(1 - power)) * se


def run_phenotype(code, label, ncase, ncontrol):
    fp = os.path.join(DIR, "finngen_R12_%s.gz" % code)
    print("\n" + "#" * 116)
    print("# %s  (ncase=%d, ncontrol=%d, K=%.4f%%)" % (label, ncase, ncontrol, 100 * ncase / (ncase + ncontrol)))
    print("#" * 116)
    if not os.path.exists(fp):
        print("  !! missing file", fp); return [], None
    print("  file %.1f MB" % (os.path.getsize(fp) / 2 ** 20))
    cols, reg, byrs = scan(fp, REGIONS, WANT_RS)
    print("  columns: %s" % cols)
    for nm, ch, lo, hi in REGIONS:
        print("    region %-8s chr%s:%-11d-%-11d  %5d variants" % (nm, ch, lo, hi, len(reg[nm])))

    rows = []
    inst = json.load(open(os.path.join(DIR, "MR_il6_saa_instruments.json"), encoding="utf-8"))["instruments"]

    # ---------- (1) Locus signal characteristics ----------
    print("\n  (1) Locus signals (regional min p; see pheno_n for the number of GW-significant variants in this FinnGen phenotype)")
    locus = []
    for nm, ch, lo, hi in REGIONS:
        d = pd.DataFrame(reg[nm])
        if d.empty:
            print("     %-8s no data in region" % nm); continue
        d = d.dropna(subset=["p"])
        best = d.nsmallest(1, "p").iloc[0]
        gw = "*GW" if best["p"] < 5e-8 else ("+" if best["p"] < 1e-5 else "")
        print("     %-8s min p=%.3g %-4s  chr%s:%-11d %s>%s  nearest gene=%s" % (
            nm, best["p"], gw, best["chrom"], best["pos"], best["ref"], best["alt"], str(best["nearest"])[:22]))
        locus.append(dict(phenotype=label, region=nm, min_p=best["p"], pos=int(best["pos"]),
                          ref=best["ref"], alt=best["alt"], beta=best["beta"], nearest=best["nearest"]))

    # ---------- (2) Positive / direction controls ----------
    print("\n  (2) Positive / direction controls (must reproduce the known direction)")
    ctrl = [("rs4148323", "UGT1A1*6 loss-of-function", "2", +1, "-> bilirubin up, jaundice risk up"),
            ("rs887829",  "UGT1A1",         "2", None, ""),
            ("rs855791",  "TMPRSS6 A",      "22", +1, "-> hepcidin down -> iron up -> heme turnover up -> bilirubin up"),
            ("rs1800562", "HFE C282Y",      "6", None, "-> iron overload (control)")]
    for rs, gene, ch, expsign, note in ctrl:
        recs = [r for r in byrs.get(rs, []) if r["chrom"] == ch]
        if not recs:
            print("     %-12s %-18s not found" % (rs, gene)); continue
        r = recs[0]
        flag = ""
        if expsign is not None:
            flag = "OK" if np.sign(r["beta"]) == expsign else "X direction opposite"
        print("     %-12s %-18s chr%s:%-11d %s>%s beta=%+.4f se=%.4f p=%.3g af=%.4f %s %s" % (
            rs, gene, r["chrom"], r["pos"], r["ref"], r["alt"], r["beta"], r["se"], r["p"], r["af"], flag, note))
        rows.append(dict(phenotype=label, block="control", rsid=rs, gene=gene, chrom=r["chrom"], pos=r["pos"],
                         ref=r["ref"], alt=r["alt"], beta=r["beta"], se=r["se"], p=r["p"], af=r["af"]))

    # ---------- (3) Main analysis IL6R -> phenotype ----------
    print("\n  (3) IL6R cis-pQTL -> %s (Wald ratio, alleles harmonized to the sIL6R-raising allele)" % label)
    print("     Expected direction: sIL6R up (= IL-6 signaling down) -> if the hepcidin-iron-erythroid axis holds, **jaundice/yellowing risk down**;")
    print("               if only the hepatic UGT1A1 arm holds, **risk up** (the two models predict opposite directions, so this test can distinguish the mechanisms!)")
    print("     %-14s %-16s %9s %8s %8s %10s %10s %10s %8s %6s" % (
        "rsID", "source dataset", "b_ex", "se_ex", "b_out", "Wald", "95%CI low", "95%CI high", "p", "F"))
    for key, ivs in inst.items():
        if not key.startswith("IL6R"):
            continue
        for iv in ivs:
            rs = iv["rsid"]
            recs = [r for r in byrs.get(rs, []) if r["chrom"] == iv["chr"]]
            if not recs:
                print("     %-14s %-16s not found" % (rs, key[:16])); continue
            r = recs[0]
            res = align(iv["beta"], iv["se"], r["beta"], r["se"], iv["ea"], iv["nea"], r["alt"], r["ref"])
            if res is None:
                print("     %-14s %-16s alleles cannot be harmonized (exposure %s>%s, outcome %s>%s)" % (
                    rs, key[:16], iv["ea"], iv["nea"], r["alt"], r["ref"])); continue
            print("     %-14s %-16s %+9.4f %8.4f %+8.4f %+10.4f %+10.4f %+10.4f %8.3g %6.0f" % (
                rs, key[:16], iv["beta"], iv["se"], r["beta"], res["ratio"],
                res["lo"], res["hi"], res["p"], res["F"]))
            rows.append(dict(phenotype=label, block="main_IL6R", rsid=rs, dataset=key,
                             b_ex=iv["beta"], se_ex=iv["se"], ea=iv["ea"], nea=iv["nea"],
                             b_out=r["beta"], se_out=r["se"], p_out=r["p"], F=res["F"],
                             wald=res["ratio"], wald_lo=res["lo"], wald_hi=res["hi"], p_mr=res["p"]))

    # ---------- (4) Power conversion ----------
    print("\n  (4) Power conversion (threshold model; expected value taken as the adult MR IL6R -> bilirubin effect of +0.0192 SD/allele)")
    for iv in inst.get("IL6R_pQTL_GCST90012005", inst.get("IL6R", [])):
        rs = iv["rsid"]
        recs = [r for r in byrs.get(rs, []) if r["chrom"] == iv["chr"]]
        if not recs:
            continue
        r = recs[0]
        exp_lor, K, i = expected_logor(0.0192, ncase, ncontrol)
        se = r["se"]
        pw = power_at(exp_lor, se)
        md = mde(se)
        print("     %-12s expected log-OR=%+.4f (OR=%.3f)  this SNP SE=%.4f -> log-OR detectable at 80%% power=%.4f (OR=%.3f)" % (
            rs, exp_lor, math.exp(exp_lor), se, md, math.exp(md)))
        print("               -> power = %.1f%%   expected/detectable ratio = %.2f" % (100 * pw, exp_lor / md))
        rows.append(dict(phenotype=label, block="power", rsid=rs, K=K, i=i,
                         exp_logOR=exp_lor, se=se, power=pw, mde=md))
        break

    return rows, locus


def main():
    # Positive control: the known effect from the adult bilirubin GWAS at UGT1A1 (used to illustrate the "ruler")
    print("=" * 116)
    print("IL6R / UGT1A1 / HAMP / TMPRSS6 -> FinnGen R12 neonatal phenotypes (full-file version, GRCh38)")
    print("=" * 116)

    all_rows, all_locus = [], []
    # Run all 5 phenotypes (the earliest 2 were run separately; re-run here to obtain a single consistent results file)
    for code, label, ncase, ncontrol in PHENO:
        rows, locus = run_phenotype(code, label, ncase, ncontrol)
        all_rows += rows
        if locus:
            all_locus += locus

    if all_rows:
        pd.DataFrame(all_rows).to_csv(os.path.join(DIR, "MR_finngen_neonatal_full_results.csv"),
                                      index=False, encoding="utf-8-sig")
        print("\nResults saved: MR_finngen_neonatal_full_results.csv (%d rows)" % len(all_rows))
    if all_locus:
        pd.DataFrame(all_locus).to_csv(os.path.join(DIR, "finngen_neonatal_locus_summary.csv"),
                                       index=False, encoding="utf-8-sig")
        print("Locus summary saved: finngen_neonatal_locus_summary.csv")


if __name__ == "__main__":
    main()
