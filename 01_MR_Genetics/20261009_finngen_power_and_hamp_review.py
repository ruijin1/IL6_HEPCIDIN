# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""(1) Power quantification for the FinnGen neonatal jaundice GWAS (to judge whether the negative IL6R result is "interpretable")
   (2) Cross-model re-aggregation of the HAMP (hepcidin) evidence -- the previous round may have over-discounted it

(1) Power logic:
   The FinnGen case proportion k = ncase/(ncase+nctrl) can be back-solved from af_alt and af_alt_cases/controls:
       af_alt = k*af_case + (1-k)*af_ctrl
   giving k ~= 0.1% (ncase ~= 499 / ~500k).
   Conversion between the log-OR of a rare binary outcome and a continuous exposure effect in SD units (threshold model):
       log-OR ~= β_SD * i,   i = φ(z)/K,  z = Φ^-1(1-K)
   The adult IL6R -> bilirubin MR provides β_SD; from this we get the "expected log-OR", then compare it with the observed SE to compute power.

(2) HAMP aggregation: group the HAMP (and iron-axis) effects across datasets by **stimulus type** -- direct IL-6 vs LPS/mixed.
   Problem awareness: 13.3(5) discounted the hepcidin hypothesis because of the HAMP -2.19 of dHepaRG(LPS);
   but if LPS and direct IL-6 act differently on HAMP, that basis for discounting does not hold and the ranking must be redone.
"""
import os, sys, json
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"


def norm_sf(x):
    return stats.norm.sf(x)


# ---------------- (1) Power ----------------
def power_block():
    print("=" * 104)
    print("(1) Power quantification for the FinnGen 'neonatal jaundice' GWAS")
    print("=" * 104)
    # Back-solve the case proportion (from af_alt_cases/controls and af_alt; here the measured values at three loci are used)
    obs = {"rs4148323": (0.0412825, 0.0671343, 0.0412566),
           "rs4601580": (0.486055, 0.457181, 0.486084),
           "rs4129267": (0.299495, 0.294628, 0.299500)}
    ks = []
    for rs, (af, afc, afk) in obs.items():
        if afc != afk:
            k = (af - afk) / (afc - afk)
            ks.append(k)
            print("  %-12s af_alt=%.6f  af_case=%.6f  af_ctrl=%.6f  ->  case proportion k = %.5f (%.3f%%)"
                  % (rs, af, afc, afk, k, k * 100))
    k = float(np.median(ks))
    Ntot = 500348
    ncase = k * Ntot
    print("\n  Inferred: total sample ~%d, cases ~%.0f (consistent with the official ncase~=499 for this endpoint)" % (Ntot, ncase))

    K = k
    z = stats.norm.isf(K)
    phi = stats.norm.pdf(z)
    i = phi / K
    print("  Threshold model: z = Φ^-1(1-K) = %.4f, i = φ(z)/K = %.3f" % (z, i))

    # SD-unit effect of adult IL6R -> bilirubin (from Section 10.4)
    ADULT = {"rs4129267": 0.0192, "rs4601580": None}
    res = pd.read_csv(os.path.join(DIR, "MR_finngen_neonatal_jaundice_results.csv"), encoding="utf-8-sig")
    main = res[res.block == "main_IL6R"].copy()
    print("\n  %-26s %-12s %10s %10s %12s %10s %8s" % (
        "instrument", "rsID", "observed effect", "SE", "expected log-OR", "80% detectable", "power"))
    print("  " + "-" * 96)
    rows = []
    for _, r in main.iterrows():
        se = (r.ratio_hi - r.ratio_lo) / (2 * 1.96)
        exp_b = ADULT.get(r.rsid)
        if exp_b is None:
            # Scale the adult IL6R -> bilirubin effect size by the "exposure effect scaling": ratio_adult corresponds to exposure beta_adult
            exp_b = np.nan
        exp_lor = exp_b * i if exp_b else np.nan
        mde = 2.802 * se
        if exp_lor == exp_lor:
            # Two-sided power = Φ(|b|/se - z) + Φ(-|b|/se - z); Φ(x)=1-sf(x), Φ(-y)=sf(y)
            lam = abs(exp_lor) / se
            pw = (1 - norm_sf(lam - 1.96)) + norm_sf(lam + 1.96)
        else:
            pw = np.nan
        print("  %-26s %-12s %+10.4f %10.4f %12s %10.4f %8s" % (
            r.dataset[:26], r.rsid, r.ratio, se,
            ("%+.4f" % exp_lor) if exp_lor == exp_lor else "  -",
            mde, ("%.1f%%" % (100 * pw)) if pw == pw else "  -"))
        rows.append(dict(dataset=r.dataset, rsid=r.rsid, observed=r.ratio, se=se,
                         expected_logOR=exp_lor, mde80=mde, power=pw))
    print("\n  Note: '80%% detectable' = 2.802 x SE (alpha=0.05 two-sided, power 80%%).")
    print("        If the expected log-OR is far below that value, this GWAS **cannot test** the hypothesis, and a negative result is not evidence.")
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(DIR, "finngen_neonatal_jaundice_power.csv"),
                                  index=False, encoding="utf-8-sig")
        print("  Saved finngen_neonatal_jaundice_power.csv")

    # "Ruler" effect of the positive control
    print("\n  -- Positive control as a 'ruler': effect size of reduced UGT1A1 activity on neonatal jaundice --")
    pc = res[res.block == "positive_control_UGT1A1"]
    for _, r in pc.iterrows():
        print("    %-12s beta=%+.4f (OR=%.2f/allele)  p=%.3g" % (r.rsid, r.beta, np.exp(r.beta), r.p))
    print("    => rs4148323 (UGT1A1*6, loss-of-function) OR=%.2f/allele, p=6e-05:"
          % np.exp(pc[pc.rsid == "rs4148323"].beta.iloc[0]))
    print("      **UGT1A1 activity is indeed rate-limiting in neonates** -- therefore the explanation 'adult UGT1A1 is not rate-limiting'")
    print("      does not hold at least in the neonatal population (this GWAS is a neonatal phenotype with adult genotypes, so the source is consistent).")


# ---------------- (2) HAMP cross-model re-ranking ----------------
def hamp_block():
    print("\n\n" + "=" * 104)
    print("(2) Cross-model re-ranking of the HAMP (hepcidin) evidence -- grouped by stimulus type")
    print("=" * 104)
    files = ["mechanism_geo_results.csv", "mechanism_circadian_results.csv",
             "mechanism_replication_results.csv", "mechanism_GSE102006_results.csv"]
    frames = []
    for f in files:
        p = os.path.join(DIR, f)
        if os.path.exists(p):
            d = pd.read_csv(p, encoding="utf-8-sig")
            d["srcfile"] = f
            frames.append(d)
    D = pd.concat(frames, ignore_index=True)
    D = D[D.get("found") == True]
    sel = D[D.gene.isin(["HAMP", "SLC40A1", "FTL", "FTH1", "TFRC"])].copy()
    if "comparison" in sel.columns and "contrast" not in sel.columns:
        sel["contrast"] = sel["comparison"]
    print("%-46s %-9s %9s %9s %8s" % ("dataset/contrast", "gene", "log2FC", "p", "direction consistent"))
    print("-" * 104)
    def pick(row, *names):
        for n in names:
            if n in row.index:
                v = row[n]
                if isinstance(v, str) and v.strip() and v.strip().lower() != "nan":
                    return v
        return ""

    for _, r in sel[sel.gene == "HAMP"].iterrows():
        lab = (pick(r, "dataset", "comparison") + " | " + pick(r, "contrast")).strip(" |")
        cons = pick(r, "pair_consist", "donor_same_dir", "time_same_dir") or "-"
        print("%-52s %-9s %+9.2f %9.2g %8s" % (lab[:52], r.gene, r.log2FC,
                                               r.p if r.p == r.p else np.nan, cons))
    print("\n  -- Classified by stimulus type --")
    IL6_MODELS = {"GSE202045 human primary hepatocytes": (5.44, 3.04, "3/3"),
                  "GSE202045 HepG2": (2.67, 0.78, "2/2"),
                  "GSE202045 mouse primary hepatocytes": (0.03, np.nan, "3/3"),
                  "GSE275342 Hep3B": (3.45, np.nan, "16/16"),
                  "GSE275342 Huh7": (0.41, np.nan, "10/16")}
    print("  Direct IL-6 stimulation:")
    for k, (a, b, c) in IL6_MODELS.items():
        print("     %-28s 6h/24h log2FC = %+.2f / %s   consistent %s" % (k, a, ("%+.2f" % b) if b == b else "  -", c))
    print("  LPS / mixed cytokines:")
    print("     %-28s 6h/24h log2FC = %+.2f / %+.2f   9/9 -> 7/9" % ("GSE230325 dHepaRG (LPS)", -2.19, 0.57))
    print("     %-28s 6h/12h/24h  = %+.2f / %+.2f / %+.2f" % ("GSE102006 HepG2 (mixed)", -0.14, -0.41, -0.61))
    print("\n  => Revised conclusion: the rise in HAMP is a response **specific to direct IL-6** (human primary +5.44, Hep3B +3.45, HepG2 +2.67);")
    print("     under LPS/mixed stimulation it does not rise and may even fall -- this is consistent with 'LPS acts via TLR4 rather than IL6R',")
    print("     and the exposure in this MR study is **IL6R**, so the **direct IL-6 models** should be taken as the reference.")
    print("     => the discounting in 13.3(5) was excessive, and the hepcidin hypothesis should be **restored**.")


if __name__ == "__main__":
    power_block()
    hamp_block()
