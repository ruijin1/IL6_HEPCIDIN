# -*- coding: utf-8 -*-
r"""* Substitute sensitivity analyses for a four-SNP correlated cis-pQTL instrument set

Scientific question:
  The four IL6R cis-pQTLs are in partial linkage disequilibrium (pairwise r2 0.08-0.44) and come
  from four different pQTL studies. The pleiotropy-robust estimators (MR-Egger, weighted median,
  MR-PRESSO) were designed for several approximately independent variants, so in this design they
  cannot carry the weight they carry in a many-SNP analysis. The question this script answers is
  therefore not "is there pleiotropy-robust evidence" but "what else can be reported, from the same
  four variants, that bears on the same threat, and what does it show".

Inputs (written by the earlier scripts of this layer, all in the same directory):
  MR_il6r_iron_erythroid_results.csv  per-SNP harmonised b_ex / se_ex / F / b_out / se_out / ratio
                                      (20261009_mr_il6r_iron_erythroid_finngen.py)
  ld_il6r_tools.csv                   FinnGen R12 LD among the four tools, long table pos1/pos2/r/r2
                                      (20261009_check_il6r_instrument_independence.py)
  _il6r_snp_coords.json               rsid -> hg19 position (written by the same two scripts)

Outputs (one row per trait per method, unless stated):
  sensitivity_substitutes.csv         IVW_corr, single-instrument Wald range, weighted median,
                                      MR-Egger intercept and slope, I2_GX, leave-one-out range, Q
  sensitivity_loo.csv                 the leave-one-instrument-out estimates themselves
  sensitivity_unmeasured_ld.csv       the one LD pair absent from the reference panel, swept over
                                      r = -0.30, -0.15, 0, +0.15, +0.30

Estimators implemented here (no external MR package, so that every step is inspectable):
  IVW_corr   correlation-aware IVW, beta = (g' S^-1 b)/(g' S^-1 g),
             Var = (g' S^-1 Psi S^-1 g)/(g' S^-1 g)^2 with Psi_ij = rho_ij * se_i * se_j
  Wald_i     single-instrument Wald ratio with a delta-method CI (exposure SE neglected, first order)
  WM         weighted median of the ratio estimates (Bowden 2016), weights 1/se_ratio^2,
             standard error from 20,000 bootstrap resamples
  Egger      weighted linear regression with an intercept (Bowden 2015), instruments oriented to a
             positive exposure effect, sigma-adjusted standard errors, plus I2_GX for NOME
  LOO        drop one instrument, recompute IVW_corr on the remaining three

Interpretation limit, stated here because it is easy to lose in a table of numbers:
  with four variants, two of which are in LD at r2 = 0.44, the Egger intercept has two degrees of
  freedom and negligible power, and the weighted median can be dominated by a single instrument.
  These analyses therefore exclude gross directional pleiotropy; they do not exclude a pleiotropic
  contribution small enough to escape an intercept test, and none of them is used as a primary
  estimator.
"""
import os
import json
import unicodedata

import numpy as np
import pandas as pd
from scipy import stats

# Local working directory; no input table is distributed with this repository.
DIR = os.environ.get("MR_DIR", r"D:\il6_hepcidin_work\03-data\MR_finngen")
OUT = os.environ.get("MR_OUT", DIR)
TOOLS = ["rs4601580", "rs7514452", "rs6690230", "rs4129267"]

# The six traits the manuscript interprets. The LD sweep is summarised over these and not over
# all fifteen traits in the results table: several of the others sit on coefficients close to
# zero, so their relative swing is large without saying anything about the estimates that are
# actually reported. Labels are the ones stored in the result table.
INTERPRETED = [
    "MCH", "MCV", "RDW",
    "可溶性转铁蛋白受体 sTfR",
    "转铁蛋白饱和度 TSAT",
    "总胆红素（参照）",
]


def ld_matrix(names):
    """Signed correlation matrix among `names`, from the long LD table plus rsid -> position."""
    ldf = pd.read_csv(os.path.join(DIR, "ld_il6r_tools.csv"), encoding="utf-8-sig")
    coords = json.load(open(os.path.join(DIR, "_il6r_snp_coords.json"), encoding="utf-8"))
    pos = {k: coords[k]["pos"] for k in names}
    S = np.eye(len(names))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            h = ldf[(ldf.pos1 == pos[a]) & (ldf.pos2 == pos[b])]
            if len(h):
                S[i, j] = float(h.iloc[0]["r"])
    return S, pos


def ivw_corr(g, b, se_out, S):
    """Correlation-aware IVW with the one-step covariance, plus the multivariate Cochran Q."""
    Si = np.linalg.inv(S)
    den = float(g @ Si @ g)
    beta = float(g @ Si @ b) / den
    Psi = S * np.outer(se_out, se_out)
    var = float(g @ Si @ Psi @ Si @ g) / den ** 2
    se = np.sqrt(var)
    resid = b - beta * g
    Q = float(resid @ np.linalg.inv(Psi) @ resid)
    return beta, se, Q


def wmedian(ratios, w, nboot=20000, seed=20261112):
    """Weighted median (Bowden 2016) with a non-parametric bootstrap standard error."""
    rng = np.random.default_rng(seed)

    def wm(r, ww):
        o = np.argsort(r)
        r, ww = np.asarray(r)[o], np.asarray(ww)[o]
        c = np.cumsum(ww) - 0.5 * ww
        c /= np.sum(ww)
        below = np.max(np.where(c < 0.5)[0]) if np.any(c < 0.5) else 0
        return r[below] + (r[below + 1] - r[below]) * (0.5 - c[below]) / (c[below + 1] - c[below])

    est = wm(ratios, w)
    bs = [wm(rng.choice(ratios, len(ratios), replace=True),
             rng.choice(w, len(w), replace=True)) for _ in range(nboot)]
    return est, float(np.std(bs, ddof=1))


def egger(g, b, se_out):
    """MR-Egger: intercept (directional pleiotropy), slope, and I2_GX for the NOME assumption."""
    sgn = np.where(g < 0, -1.0, 1.0)          # orient every instrument to a positive exposure effect
    g, b = g * sgn, b * sgn
    X = np.column_stack([np.ones_like(g), g])
    W = np.diag(1.0 / se_out ** 2)
    XtWX = X.T @ W @ X
    coef = np.linalg.solve(XtWX, X.T @ W @ b)
    resid = b - X @ coef
    df = len(g) - 2
    sigma2 = float(resid @ W @ resid) / df if df > 0 else np.nan
    cov = np.linalg.inv(XtWX) * sigma2
    se = np.sqrt(np.diag(cov))
    Qgx = float(np.sum((g - np.average(g, weights=1 / se_out ** 2)) ** 2 / (se_out ** 2)))
    Igx = max(0.0, (Qgx - (len(g) - 1)) / Qgx) if Qgx > 0 else np.nan
    return coef[0], se[0], coef[1], se[1], Igx


def main():
    res = pd.read_csv(os.path.join(DIR, "MR_il6r_iron_erythroid_results.csv"), encoding="utf-8-sig")
    names = [t for t in TOOLS if t in set(res.rsid)]
    S, pos = ld_matrix(names)

    rows, loo_rows, ld_rows = [], [], []
    print(f"{'trait':34s} {'k':>2s} {'IVWcorr':>9s} {'p':>9s} {'WM':>9s} "
          f"{'EggerInt':>9s} {'p_int':>7s} {'sign':>5s}")
    for tr in res.trait.unique():
        s = res[res.trait == tr].set_index("rsid").reindex(names).dropna(subset=["b_out"])
        if len(s) < 3:
            continue
        g = s["b_ex"].values.astype(float)
        b = s["b_out"].values.astype(float)
        seo = s["se_out"].values.astype(float)
        ratio = b / g
        se_ratio = seo / np.abs(g)
        idx = [names.index(x) for x in s.index]
        Ss = S[np.ix_(idx, idx)]

        beta, se, Q = ivw_corr(g, b, seo, Ss)
        p = 2 * stats.norm.sf(abs(beta / se))

        wm, wm_se = wmedian(ratio, 1 / se_ratio ** 2)
        icpt, icpt_se, slope, slope_se, Igx = egger(g, b, seo)
        p_int = 2 * stats.t.sf(abs(icpt / icpt_se), max(len(g) - 2, 1))

        sign_conc = int(np.sum(np.sign(ratio) == np.sign(beta)))
        loo = []
        for drop in s.index:
            keep = [x for x in s.index if x != drop]
            ii = [names.index(x) for x in keep]
            gg = s.loc[keep, "b_ex"].values.astype(float)
            bb = s.loc[keep, "b_out"].values.astype(float)
            ss = s.loc[keep, "se_out"].values.astype(float)
            bb2, ss2, _ = ivw_corr(gg, bb, ss, S[np.ix_(ii, ii)])
            loo.append(bb2)
            loo_rows.append(dict(trait=tr, dropped=drop, beta_loo=bb2, se_loo=ss2,
                                 p_loo=2 * stats.norm.sf(abs(bb2 / ss2))))
        loo = np.array(loo)

        rows.append(dict(trait=tr, k=len(s), beta_ivw=beta, se_ivw=se, p_ivw=p,
                         beta_wm=wm, se_wm=wm_se, p_wm=2 * stats.norm.sf(abs(wm / wm_se)),
                         egger_intercept=icpt, se_intercept=icpt_se, p_intercept=p_int,
                         egger_slope=slope, se_egger_slope=slope_se, I2_GX=Igx,
                         sign_concordant=f"{sign_conc}/{len(s)}",
                         ratio_min=float(ratio.min()), ratio_max=float(ratio.max()),
                         loo_min=float(loo.min()), loo_max=float(loo.max()),
                         Q=Q, Q_p=stats.chi2.sf(Q, len(s) - 1)))
        print(f"{tr:34s} {len(s):2d} {beta:+9.4f} {p:9.1e} {wm:+9.4f} "
              f"{icpt:+9.5f} {p_int:7.3f} {sign_conc}/{len(s)}")

        # the one LD coefficient that had to be assumed: sweep it rather than fix it at zero
        ia, ib = names.index("rs4601580"), names.index("rs6690230")
        for r in (-0.30, -0.15, 0.0, 0.15, 0.30):
            S2 = Ss.copy()
            if ia in idx and ib in idx:
                pa, pb = idx.index(ia), idx.index(ib)
                S2[pa, pb] = S2[pb, pa] = r
            bb2, ss2, _ = ivw_corr(g, b, seo, S2)
            ld_rows.append(dict(trait=tr, r_unmeasured=r, beta=bb2, se=ss2))

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "sensitivity_substitutes.csv"),
                              index=False, encoding="utf-8-sig")
    pd.DataFrame(loo_rows).to_csv(os.path.join(OUT, "sensitivity_loo.csv"),
                                  index=False, encoding="utf-8-sig")
    d = pd.DataFrame(ld_rows)
    d.to_csv(os.path.join(OUT, "sensitivity_unmeasured_ld.csv"), index=False, encoding="utf-8-sig")

    piv = d.pivot(index="trait", columns="r_unmeasured", values="beta")
    print("\nsensitivity to the unmeasured LD pair (rs4601580-rs6690230), beta_IVW:")
    print(piv.to_string(float_format=lambda x: f"{x:+.4f}"))
    # Swing expressed against each trait's own r = 0 estimate, so the numbers below correspond
    # one to one with the percentages quoted in the manuscript.
    rel = (piv.max(axis=1) - piv.min(axis=1)).abs() / piv[0.0].abs()
    rel = rel[[t for t in INTERPRETED if t in rel.index]].sort_values(ascending=False)
    print("\nswing over r = -0.30..+0.30, against each trait's own r = 0 estimate:")
    for t in rel.index:
        # the trait labels are Chinese, so pad on display width rather than character count
        pad = " " * max(0, 30 - sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in t))
        print("  %s%s %6.1f%%   (%+.4f to %+.4f)"
              % (t, pad, 100 * rel[t], piv.loc[t].min(), piv.loc[t].max()))
    print("  largest of the six: %.1f%% (%s)" % (100 * rel.max(), rel.idxmax()))


if __name__ == "__main__":
    main()
