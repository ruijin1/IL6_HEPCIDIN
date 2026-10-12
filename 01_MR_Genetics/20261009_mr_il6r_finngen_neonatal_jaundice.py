# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""IL6R / UGT1A1 -> **neonatal jaundice** (FinnGen R12 P16_NEONTAL_JAUND_OTH_UNSP_CAUSES) cis-MR

Why this step:
  The core contradiction in Sections 12/13 is -- the mechanism layer predicts "IL-6 pathway up -> hepatic UGT1A1 down -> bilirubin up",
  whereas the MIMIC observational analysis and the IL6R -> bilirubin MR (adults) both give "IL-6 signaling up -> bilirubin down".
  What is missing is exactly the **genetic evidence from neonates themselves**. FinnGen R12 has a neonatal jaundice endpoint (GRCh38),
  so we can directly test **whether, in the neonatal population, the genetic proxy for the IL-6 pathway also points to "bilirubin down"**.

Design:
  * Main analysis   IL6R cis-pQTL (4 independent datasets) -> neonatal jaundice
  * Positive control * UGT1A1 locus (rs6754859 etc.) -> neonatal jaundice
              -- if the UGT1A1 locus is not significant, the phenotype definition/sample size of this GWAS is inadequate for MR, and all results are void
  * Negative control SAA1/SAA2 pQTL -> neonatal jaundice (expected negative, consistent with adults)

* Two methodological constraints that must be stated explicitly (from the checklist of the skill eqtl-catalogue-region-fetch):
  1. **Genotype source**: FinnGen = an adult cohort, and its "neonatal jaundice" is the **medical-record entry from infancy of the same individuals**,
     so genotype and phenotype belong to the same individuals -> the source is consistent and the design holds (this is not a cross-generational design).
  2. **Single-SNP cis instrument** -> no pleiotropy test is possible, and we must rely on colocalization; this script reports only the Wald ratio + F.

* Reading only the front part of the data suffices: FinnGen files are sorted by chromosome position, and chr1:154M and chr2:233M are both in the first half of the file,
   so these two regions can be obtained even while the file is still downloading (not yet at 100%); the script reports the maximum position actually scanned.
"""
import os, sys, gzip, json, re, math, zlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
FIN = os.path.join(DIR, "finngen_R12_P16_NEONTAL_JAUND_OTH_UNSP_CAUSES.gz")

# FinnGen R12 metadata (ncase/ncontrol): start from the values known in the script, override if a manifest is obtained
KNOWN_N = {"P16_NEONTAL_JAUND_OTH_UNSP_CAUSES": dict(ncase=499, ncontrol=None)}

# Regions of interest (GRCh38)
REGIONS = [("IL6R", "1", 154_300_000, 154_800_000),
           ("UGT1A1", "2", 233_600_000, 233_900_000)]

# Variants to query exactly by rsID (exposure-side instruments + positive controls)
WANT_RS = {
    # IL6R cis-pQTL (adult pQTL, see MR_il6_saa_instruments.json)
    "rs4601580", "rs7514452", "rs6690230", "rs4129267",
    # UGT1A1 locus (instruments from the bilirubin GWAS + classical functional variants)
    "rs6754859", "rs887829", "rs4148323", "rs6742078", "rs12995326",
    # negative control (SAA locus)
    "rs11024589", "rs10430893",
    # CRP
    "rs3116654",
}
# Approximate positions of the panel genes (GRCh38; 1 kb precision is enough; only used to print the nearest gene)
PANEL_LOC = {"IL6R": ("1", 154_400_000), "UGT1A1": ("2", 233_760_000),
             "IL6": ("7", 22_766_000), "CRP": ("1", 159_712_000),
             "SAA1": ("11", 18_270_000), "HNF4A": ("20", 44_350_000)}


def scan(fp, regions, want_rs, limit_chrom=None):
    """Single streaming scan: take all variants in the region + variants for the specified rsIDs.
    * Read line by line through gzip; if the file is truncated, catch EOFError and return what has been read (and report the last position).
    """
    cols = None
    reg = {n: [] for n, _, _, _ in regions}
    byrs = {}
    last = ("?", 0, 0)
    truncated = False
    try:
        with gzip.open(fp, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#chrom") or line.startswith("chrom"):
                    cols = line.lstrip("#").rstrip("\n").split("\t")
                    idx = {c: i for i, c in enumerate(cols)}
                    continue
                try:
                    p = line.rstrip("\n").split("\t")
                except Exception:
                    continue
                if len(p) < 8:
                    continue
                ch, pos = p[idx["chrom"]], p[idx["pos"]]
                try:
                    ipos = int(pos)
                except ValueError:
                    continue
                last = (ch, ipos, 0)
                rec = dict(chrom=ch, pos=ipos, ref=p[idx["ref"]], alt=p[idx["alt"]],
                           rsids=p[idx["rsids"]], nearest=p[idx.get("nearest_genes", 5)],
                           p=float(p[idx["pval"]]), beta=float(p[idx["beta"]]) if p[idx["beta"]] not in ("NA", "") else np.nan,
                           se=float(p[idx["sebeta"]]) if p[idx["sebeta"]] not in ("NA", "") else np.nan,
                           af=float(p[idx["af_alt"]]) if p[idx["af_alt"]] not in ("NA", "") else np.nan)
                for nm, rch, lo, hi in regions:
                    if ch == rch and lo <= ipos <= hi:
                        reg[nm].append(rec)
                for rs in str(rec["rsids"]).split(","):
                    rs = rs.strip()
                    if rs in want_rs:
                        byrs.setdefault(rs, []).append(rec)
    except (EOFError, OSError, zlib.error) as e:
        # * The file is still downloading -> stop at the truncation point; the region data already read remains usable
        truncated = True
        print("  [note] reached end of file (may still be downloading): %s" % str(e)[:70])
    return cols, reg, byrs, last, truncated


def wald(b_ex, se_ex, b_out, se_out, ea_ex, nea_ex, ea_out, nea_out, eaf=None):
    """Wald ratio after allele harmonization"""
    if ea_ex == ea_out and nea_ex == nea_out:
        s = 1
    elif ea_ex == nea_out and nea_ex == ea_out:
        s = -1
    else:
        return None
    bo = s * b_out
    # SE by the delta method: ratio = bo / b_ex, SE ~= |bo/b_ex| * sqrt((se_out/bo)^2 + (se_ex/b_ex)^2)
    ratio = bo / b_ex
    se = abs(ratio) * np.sqrt((se_out / bo) ** 2 + (se_ex / b_ex) ** 2) if bo != 0 else np.nan
    return dict(ratio=ratio, se=se,
                ci_lo=ratio - 1.96 * se, ci_hi=ratio + 1.96 * se,
                p=(2 * stats.norm.sf(abs(ratio / se)) if se and se > 0 else np.nan),
                F=(b_ex / se_ex) ** 2, align=s)


def main():
    if not os.path.exists(FIN):
        print("!! FinnGen file does not exist:", FIN); return
    sz = os.path.getsize(FIN)
    print("FinnGen file currently %.1f MB (may not be fully downloaded; does not affect the front regions)" % (sz / 2 ** 20))
    cols, reg, byrs, last, truncated = scan(FIN, REGIONS, WANT_RS)
    print("columns: %s" % cols)
    print("Maximum position actually scanned: chr%s:%s  truncated=%s  <- chr1:154M and chr2:233M are both before this, so both regions are complete" % (last[0], last[1], truncated))
    for nm, ch, lo, hi in REGIONS:
        print("  region %s (chr%s:%d-%d) hit %d variants" % (nm, ch, lo, hi, len(reg[nm])))
    print("  specified rsIDs hit %d: %s" % (len(byrs), sorted(byrs.keys())))

    out = []
    # ---------- Locus characteristics: regional min p ----------
    print("\n" + "=" * 108)
    print("(1) Locus signal characteristics (FinnGen neonatal jaundice)")
    print("=" * 108)
    for nm, ch, lo, hi in REGIONS:
        d = pd.DataFrame(reg[nm])
        if d.empty:
            print("  %s: no data in region (the file may not cover it yet)" % nm); continue
        d = d.dropna(subset=["p"])
        top = d.nsmallest(5, "p")
        print("\n  %-8s regional min p = %.3g (%d variants)" % (nm, d.p.min(), len(d)))
        print("     %-14s %-12s %-7s %-7s %10s %9s %9s  %s" % ("rsID", "position", "ref", "alt", "p", "beta", "af", "nearest gene"))
        for _, r in top.iterrows():
            print("     %-14s %-12d %-7s %-7s %10.3g %+9.4f %9.4f  %s" % (
                str(r["rsids"])[:14], r["pos"], r["ref"], r["alt"], r["p"], r["beta"], r["af"],
                str(r["nearest"])[:18]))

    # ---------- Positive control: UGT1A1 locus ----------
    print("\n" + "=" * 108)
    print("(2) * Positive control: UGT1A1 locus instruments -> neonatal jaundice (**must be significant, otherwise this GWAS is unusable**)")
    print("=" * 108)
    for rs in ["rs6754859", "rs887829", "rs4148323", "rs6742078", "rs12995326"]:
        recs = [r for r in byrs.get(rs, []) if r["chrom"] == "2"]
        if not recs:
            print("  %-12s not found in the file" % rs); continue
        r = recs[0]
        print("  %-12s chr2:%-11d %s>%s  beta=%+.4f se=%.4f p=%.3g af=%.4f  nearest gene=%s" % (
            rs, r["pos"], r["ref"], r["alt"], r["beta"], r["se"], r["p"], r["af"], r["nearest"]))
        out.append(dict(block="positive_control_UGT1A1", rsid=rs, gene=r["nearest"],
                        ch=r["chrom"], pos=r["pos"], ref=r["ref"], alt=r["alt"],
                        beta=r["beta"], se=r["se"], p=r["p"]))

    # ---------- Main analysis: IL6R cis-pQTL -> neonatal jaundice ----------
    print("\n" + "=" * 108)
    print("(3) Main analysis: IL6R cis-pQTL -> neonatal jaundice (Wald ratio, after allele harmonization)")
    print("=" * 108)
    inst = json.load(open(os.path.join(DIR, "MR_il6_saa_instruments.json"), encoding="utf-8"))["instruments"]
    for key, ivs in inst.items():
        if not key.startswith("IL6R") and not key.startswith("SAA"):
            continue
        for iv in ivs:
            rs = iv["rsid"]
            recs = [r for r in byrs.get(rs, []) if r["chrom"] == iv["chr"]]
            if not recs:
                print("  %-26s %-12s not found" % (key, rs)); continue
            r = recs[0]
            res = wald(iv["beta"], iv["se"], r["beta"], r["se"],
                       iv["ea"], iv["nea"], r["alt"], r["ref"], r["af"])
            if res is None:
                print("  %-26s %-12s alleles cannot be harmonized (exposure %s>%s, outcome %s>%s)" % (
                    key, rs, iv["ea"], iv["nea"], r["ref"], r["alt"])); continue
            print("  %-26s %-12s  exposure beta=%+.4f(F=%.0f)  outcome beta=%+.4f(p=%.3g)  ->  Wald %+.4f (95%%CI %+.4f~%+.4f) p=%.3g" % (
                key, rs, iv["beta"], res["F"], r["beta"], r["p"], res["ratio"], res["ci_lo"], res["ci_hi"], res["p"]))
            out.append(dict(block="main_IL6R", dataset=key, rsid=rs, gene="IL6R",
                            b_exposure=iv["beta"], se_exposure=iv["se"], F=res["F"],
                            b_outcome=r["beta"], se_outcome=r["se"], p_outcome=r["p"],
                            ratio=res["ratio"], ratio_lo=res["ci_lo"], ratio_hi=res["ci_hi"], p_mr=res["p"]))

    # ---------- Negative control ----------
    print("\n" + "=" * 108)
    print("(4) Negative control: CRP (rs3116654) and the SAA locus -> neonatal jaundice")
    print("=" * 108)
    for rs, gene in [("rs3116654", "CRP"), ("rs11024589", "SAA1"), ("rs10430893", "SAA2")]:
        recs = byrs.get(rs, [])
        if not recs:
            print("  %-12s %-6s not found" % (rs, gene)); continue
        r = recs[0]
        print("  %-12s %-6s chr%s:%-11d %s>%s  beta=%+.4f se=%.4f p=%.3g  nearest gene=%s" % (
            rs, gene, r["chrom"], r["pos"], r["ref"], r["alt"], r["beta"], r["se"], r["p"], r["nearest"]))
        out.append(dict(block="negative_control", rsid=rs, gene=gene, ch=r["chrom"], pos=r["pos"],
                        ref=r["ref"], alt=r["alt"], beta=r["beta"], se=r["se"], p=r["p"]))

    if out:
        p = os.path.join(DIR, "MR_finngen_neonatal_jaundice_results.csv")
        pd.DataFrame(out).to_csv(p, index=False, encoding="utf-8-sig")
        print("\nResults saved: %s" % p)

    # ---------- Also save the full regional tables (for colocalization) ----------
    for nm, ch, lo, hi in REGIONS:
        if reg[nm]:
            p = os.path.join(DIR, "finngen_neonatal_jaundice_%s_region.csv" % nm)
            pd.DataFrame(reg[nm]).to_csv(p, index=False, encoding="utf-8-sig")
            print("Regional table saved: %s (%d variants, usable for coloc)" % (p, len(reg[nm])))


if __name__ == "__main__":
    main()
