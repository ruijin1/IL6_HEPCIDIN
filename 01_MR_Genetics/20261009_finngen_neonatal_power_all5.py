# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""Power comparison across the 5 FinnGen R12 neonatal phenotypes + summary of the strongest regional signal

Why this is done separately: the power section in `20261009_mr_il6r_finngen_full_neonatal.py` did not print because of an exposure key name mismatch
(the script used `IL6R_pQTL_GCST90012005`, whereas the actual key is `IL6R:ebi-a-GCST90012005`).
This script recomputes directly from the saved results CSV and compares the 5 phenotypes in a single table.

Threshold model (neonatal phenotypes have a very small K -> i = φ(z)/K):
  expected log-OR = b_SD_per_allele x i      (this project uses the adult MR estimate for IL6R -> bilirubin, +0.0192 SD/allele)
  log-OR detectable at 80% power = (1.96 + 0.8416) x SE
  power = Φ(|exp|/SE - 1.96) + Φ(-|exp|/SE - 1.96)
"""
import os, sys, math
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
B_SD = 0.0192          # adult IL6R -> bilirubin (SD/allele)
PHENO = [("neonatal jaundice", 499, 499688, 0), ("neonatal bacterial sepsis", 203, 499933, 1),
         ("neonatal respiratory distress", 498, 499476, 0), ("chronic respiratory disease of perinatal origin (BPD)", 69, 499476, 0),
         ("preterm birth/low birth weight", 782, 499424, 3)]

# The saved result files predate the English release and store the phenotype names in the
# original form. Map them to the canonical English names so the comparisons above work
# against either the old files or freshly regenerated ones.
PHENO_ALIAS = {
    "新生儿黄疸": "neonatal jaundice",
    "新生儿细菌性败血症": "neonatal bacterial sepsis",
    "新生儿呼吸窘迫": "neonatal respiratory distress",
    "围产期慢性呼吸疾病(BPD)": "chronic respiratory disease of perinatal origin (BPD)",
    "早产/低出生体重": "preterm birth/low birth weight",
}


def main():
    R = pd.read_csv(os.path.join(DIR, "MR_finngen_neonatal_full_results.csv"), encoding="utf-8-sig")
    R["phenotype"] = R.phenotype.replace(PHENO_ALIAS)
    M = R[R.block == "main_IL6R"].copy()
    L = pd.read_csv(os.path.join(DIR, "finngen_neonatal_locus_summary.csv"), encoding="utf-8-sig")
    if "phenotype" in L.columns:
        L["phenotype"] = L.phenotype.replace(PHENO_ALIAS)

    print("=" * 116)
    print("(1) IL6R cis-pQTL -> phenotype Wald ratio for each phenotype (alleles aligned to the sIL6R-raising allele)")
    print("=" * 116)
    print("%-22s %-12s %9s %9s %22s %9s %8s" % ("phenotype", "rsID", "Wald", "se", "95%CI", "p", "F"))
    print("-" * 116)
    for lab in M.phenotype.unique():
        s = M[M.phenotype == lab]
        for _, r in s.iterrows():
            se = (r.wald_hi - r.wald_lo) / (2 * 1.96)
            print("%-22s %-12s %+9.4f %9.4f  (%+.4f,%+.4f) %9.3g %8.0f"
                  % (lab[:22], r.rsid, r.wald, se, r.wald_lo, r.wald_hi, r.p_mr, r.F))
        print()

    print("=" * 116)
    print("(2) Power comparison (expecting the adult IL6R -> bilirubin effect of +0.0192 SD/allele; to show that 'it is not that there is no effect, but that it cannot be detected')")
    print("=" * 116)
    rows = []
    for lab, ncase, nctrl, ngw in PHENO:
        K = ncase / (ncase + nctrl)
        z = stats.norm.isf(K); i = stats.norm.pdf(z) / K
        exp_lor = B_SD * i
        s = M[M.phenotype == lab]
        if not len(s):
            continue
        se = float(((s.wald_hi - s.wald_lo) / (2 * 1.96)).median())
        mde = (1.96 + stats.norm.isf(0.20)) * se
        pw = float(stats.norm.cdf(abs(exp_lor) / se - 1.96) + stats.norm.cdf(-abs(exp_lor) / se - 1.96))
        rows.append(dict(phenotype=lab, ncase=ncase, K_pct=100 * K, inflation_i=i,
                         exp_logOR=exp_lor, OR_exp=math.exp(exp_lor),
                         median_SE=se, MDE80=mde, OR_MDE80=math.exp(mde),
                         power_pct=100 * pw, finngen_gw_hits=ngw))
    P = pd.DataFrame(rows)
    print(P.to_string(index=False, float_format=lambda x: "%.4g" % x))
    print("\n* Reading: the smaller K is -> the larger the inflation factor i of the threshold model -> a binary phenotype **inflates** the per-allele log-OR;")
    print("  but the neonatal phenotypes have a too-small ncase -> too-large SE -> power is still extremely low. **All phenotypes except preterm birth/low birth weight are underpowered.**")
    P.to_csv(os.path.join(DIR, "finngen_neonatal_power_all5.csv"), index=False, encoding="utf-8-sig")

    print("\n" + "=" * 116)
    print("(3) Regional minimum p for the 5 regions (note: **the regional minimum p is not the p of the target gene**; check the nearest gene)")
    print("=" * 116)
    piv = L.pivot_table(index="region", columns="phenotype", values="min_p", aggfunc="min")
    piv = piv.reindex(["IL6R", "UGT1A1", "HFE", "HAMP", "TMPRSS6"])
    print(piv.to_string(float_format=lambda x: "%.2g" % x))
    print("\n   Locus and nearest gene of the minimum p in each region (only those with p<1e-4):")
    for _, r in L[L.min_p < 1e-4].sort_values("min_p").iterrows():
        print("     %-8s %-22s chr%s:%-11d p=%.2g  nearest gene=%s"
              % (r.region, r.phenotype[:22], "?", int(r.pos), r.min_p, str(r.nearest)[:24]))
    print("\n* Number of FinnGen genome-wide significant loci for each phenotype itself: jaundice 0 / sepsis 1 / RDS 0 / BPD 0 / preterm-low-birth-weight 3")
    print("  => all the regional signals above (minimum p ~2e-5) **do not reach genome-wide significance (5e-8)** and can only serve as leads to be validated.")


if __name__ == "__main__":
    main()
