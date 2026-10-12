# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""(4) IVW meta-analysis of **correlated instruments** using the official FinnGen LD matrix (the two SNPs with r^2=0.44 cannot be treated as independent)

Why this must be done: the main table of Section 15 uses the **per-SNP Wald ratio** of the 4 instruments;
but in practice some pairs reach r^2 = 0.44 (rs4129267-rs6690230), so they are **not independent instruments**;
a simple average would under/over-estimate precision, and "4/4 concordant" cannot be taken as 4 independent validations.

Approach (the same as the standard approach when individual-level data are unavailable):
  Let g = the vector of per-SNP effects on the exposure, b = the vector of effects on the outcome, Σ = the LD **correlation matrix** among instruments (using r, not r^2)
  β_IVW = (g^T Σ^-1 b) / (g^T Σ^-1 g)
  Under the approximation that "the outcome is a standardized phenotype",  se(β) ~= 1/sqrt(g^T Σ^-1 g)
  (a more rigorous treatment requires the residual variance of the outcome; here the standardized-phenotype approximation is used and noted in the report)
"""
import os, sys, json
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"


def main():
    coords = json.load(open(os.path.join(DIR, "_il6r_snp_coords.json"), encoding="utf-8"))
    inst = json.load(open(os.path.join(DIR, "MR_il6_saa_instruments.json"), encoding="utf-8"))["instruments"]
    tools = []
    for k, v in inst.items():
        if k.startswith("IL6R"):
            for iv in v:
                c = coords.get(iv["rsid"])
                if c:
                    tools.append(dict(rsid=iv["rsid"], pos=c["pos"], b_ex=iv["beta"], se_ex=iv["se"],
                                      ea=iv["ea"], nea=iv["nea"], dataset=k))
    names = [t["rsid"] for t in tools]
    P = [t["pos"] for t in tools]
    print("Instruments (%d):" % len(names))
    for t in tools:
        print("  %-12s chr1:%-10d b_ex=%+.4f (F=%.0f)  %s" % (t["rsid"], t["pos"], t["b_ex"],
                                                             (t["b_ex"] / t["se_ex"]) ** 2, t["dataset"]))

    # ---- Build the correlation matrix Σ from the LD table ----
    ld = pd.read_csv(os.path.join(DIR, "ld_il6r_tools_vs_region.csv"), encoding="utf-8-sig")
    # This table is p1/p2/r2 (only the pairs involving the instruments are kept); switched to rebuilding from the full table
    n = len(P)
    S = np.eye(n)
    print("\nLD correlation matrix among instruments (using the official FinnGen R12 r, 520,210 samples)")
    # * Build the signed correlation matrix from ld_il6r_tools.csv (long table pos1/pos2/r/r2)
    ldf = pd.read_csv(os.path.join(DIR, "ld_il6r_tools.csv"), encoding="utf-8-sig")
    for i, pi in enumerate(P):
        for j, pj in enumerate(P):
            if i == j:
                continue
            h = ldf[((ldf.pos1 == pi) & (ldf.pos2 == pj))]
            if len(h):
                S[i, j] = float(h.iloc[0]["r"])
    print("  " + " ".join("%11s" % x for x in names))
    for i, nm in enumerate(names):
        print("  %-10s" % nm + " ".join("%11.3f" % S[i, j] for j in range(n)))

    try:
        Si = np.linalg.inv(S)
    except np.linalg.LinAlgError:
        print("  !! correlation matrix is not invertible"); return
    g = np.array([t["b_ex"] for t in tools])
    w = Si @ np.ones(n)
    denom = float(g @ Si @ g)
    print("\n  g^T Σ^-1 g = %.4f (information; equals approximately Σ g^2/se^2 when independent)" % denom)

    res = pd.read_csv(os.path.join(DIR, "MR_il6r_iron_erythroid_results.csv"), encoding="utf-8-sig")
    print("\n" + "=" * 104)
    print("IVW meta-analysis of correlated instruments (vs. the simple average of the 4 single-SNP Wald ratios)")
    print("=" * 104)
    print("%-24s %10s %10s %22s %10s %9s %8s %9s" % ("trait", "β_IVW", "se", "95%CI", "p", "mean single SNP", "Q(df=3)", "Q p"))
    out = []
    for tr in res.trait.unique():
        s = res[res.trait == tr].set_index("rsid")
        s = s.reindex(names).dropna(subset=["b_out"])
        if len(s) < 3:
            continue
        bs = s["b_out"].values
        se_o = s["se_out"].values
        # Re-fetch g in the same order
        gg = np.array([next(t["b_ex"] for t in tools if t["rsid"] == r) for r in s.index])
        Ss = Si_small(s.index.values, names, S)          # LD correlation matrix for this subset
        Ssi = np.linalg.inv(Ss)
        num = float(gg @ Ssi @ bs)
        dn = float(gg @ Ssi @ gg)
        beta = num / dn
        # * Standard error of the correlated-instrument IVW: Var(β) = (g^T Σ^-1 Ψ Σ^-1 g) / (g^T Σ^-1 g)^2,
        #   where Ψ_ij = ρ_ij * se_i * se_j is the covariance matrix of the outcome effect estimates
        Psi = Ss * np.outer(se_o, se_o)
        var = float(gg @ Ssi @ Psi @ Ssi @ gg) / (dn ** 2)
        se = np.sqrt(var) if var > 0 else np.nan
        p = 2 * stats.norm.sf(abs(beta / se)) if se and se > 0 else np.nan
        # * Cochran's Q for correlated instruments (multivariate form): Q = (b-βg)^T Ψ^-1 (b-βg), df = n-1
        try:
            resid = bs - beta * gg
            Q = float(resid @ np.linalg.inv(Psi) @ resid)
            Qp = float(stats.chi2.sf(Q, len(s) - 1))
        except Exception:
            Q, Qp = np.nan, np.nan
        mean_single = float(np.mean(s["ratio"].values))
        print("%-24s %+10.4f %10.4f (%+.4f,%+.4f) %10.2g %+10.4f %8.1f %9.2g"
              % (tr[:24], beta, se, beta - 1.96 * se, beta + 1.96 * se, p, mean_single, Q, Qp))
        out.append(dict(trait=tr, beta_ivw=beta, se_ivw=se, lo=beta - 1.96 * se, hi=beta + 1.96 * se,
                        p_ivw=p, mean_single=mean_single, n_snp=len(s), Q=Q, Q_p=Qp))
    if out:
        pd.DataFrame(out).to_csv(os.path.join(DIR, "MR_il6r_iron_IVW_correlated.csv"),
                                 index=False, encoding="utf-8-sig")
        print("\nSaved MR_il6r_iron_IVW_correlated.csv")
    print("\n* Note: Ψ_ij = ρ_ij*se_i*se_j (the LD correlation matrix used as the correlation coefficient of the outcome effect estimates), the standard approximation for correlated-instrument IVW.")
    print("  Only 4 instruments -> MR-Egger/MR-PRESSO pleiotropy tests are not possible; heterogeneity uses Cochran's Q (would need an extension).")


def Si_small(idx, names, Si):
    pos = [names.index(x) for x in idx]
    return Si[np.ix_(pos, pos)]


if __name__ == "__main__":
    main()
