# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""An LD-free colocalization criterion: **effect-ratio constancy**

Why it is needed (the key to this round's conclusion):
  - `coloc.abf` **does not use LD and allows only a single causal variant** -> even the positive control
    IL6R x CRP is called H3~=1 (wrong) -> results unusable
  - `coloc.susie` needs an LD matrix from the **same source** as the exposure GWAS: the OpenGWAS LD endpoint is
    offline; switching to the FinnGen official LD (Finnish, 520,210 samples) makes `susie_rss` **fail to converge**
    and estimate a negative residual variance -> no usable LD is available, so **a formal coloc cannot currently be run**

Alternative criterion (no LD required):
  If the exposure and the outcome share the same causal variant in this region, then for any two variants i and j
  that both measure the exposure well, we should have
    ratio_i = b_out,i / b_exp,i ~= ratio_j
  Conversely, if each has its own independent signal, the ratios will be systematically non-constant.
  Criteria: (1) the **coefficient of variation CV** of the ratios; (2) a **Cochran Q heterogeneity test** on the
  ratios (H0: the ratios are the same).

Input: coloc/ld_variant_map.csv (provides "the id of the same variant in each file") + the per-region CSVs
Output: coloc/effect_ratio_round2.csv and a printed table
"""
import os, sys, math
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
COLOC = os.path.join(DIR, "coloc")

PAIRS = [
    ("exp_IL6R_folkersen2018", "out_crp_said2022",    "IL6R x CRP (positive control, shared causal variant known)", 3394, 575531),
    ("exp_IL6R_folkersen2018", "out_mch_finngen",     "IL6R x MCH (FinnGen quantitative)", 3394, 447004),
    ("exp_IL6R_folkersen2018", "out_mcv_finngen",     "IL6R x MCV (FinnGen quantitative)", 3394, 447040),
    ("exp_IL6R_folkersen2018", "out_bil_finngen_lab", "IL6R x total bilirubin (FinnGen quantitative)", 3394, 190407),
    ("exp_IL6R_folkersen2018", "out_jaund_finngen",   "IL6R x neonatal jaundice (FinnGen P16)", 3394, 500187),
    ("exp_IL6R_folkersen2018", "out_bil_ukb_neale",   "IL6R x total bilirubin (UKB-Neale, replicating the previous round)", 3394, 342829),
]
MIN_ZEXP = 5.0      # the exposure side must be measured well
MIN_MAF = 0.05
WIN_HALF = 15000    # * take only +/-15 kb around the exposure lead (one LD block).
                    #   Using all 124 SNPs of the full region does not hold: the region has **more than one signal**
                    #   (variants such as rs9616 at the IL6R 3' end have clearly different ratios), which distorts CV/Q.


def rd(f):
    d = pd.read_csv(os.path.join(COLOC, f + ".csv"), encoding="utf-8-sig")
    for c in ("beta", "se", "eaf", "p"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["beta", "se"])


def main():
    MAP = pd.read_csv(os.path.join(COLOC, "ld_variant_map.csv"), encoding="utf-8-sig")
    rows, detail = [], []
    print("=" * 118)
    print("Effect-ratio constancy (LD-free colocalization proxy): multiple variants in the IL6R cis-pQTL region -> effect ratio to the same outcome")
    print("=" * 118)
    for expf, outf, tag, Nexp, Nout in PAIRS:
        try:
            e = rd(expf); o = rd(outf)
        except FileNotFoundError:
            print("  (missing file) %s" % tag); continue
        me = MAP[MAP.file == expf][["id", "ld_id"]]
        mo = MAP[MAP.file == outf][["id", "ld_id"]]
        e = e.merge(me, left_on="rsid", right_on="id", how="inner").drop_duplicates("ld_id")
        o = o.merge(mo, left_on="rsid", right_on="id", how="inner").drop_duplicates("ld_id")
        # * pandas suffixes produce "beta.e"/"beta.o" (with a dot), not "beta_e";
        #   rename them uniformly to _e / _o for dotted access.
        e = e.rename(columns={c: c + "_e" for c in ("rsid", "beta", "se", "eaf", "p", "position")})
        o = o.rename(columns={c: c + "_o" for c in ("rsid", "beta", "se", "eaf", "p", "position")})
        m = e.merge(o, on="ld_id")
        m["z_exp"] = m.beta_e / m.se_e
        m["maf"] = np.minimum(m.eaf_e, 1 - m.eaf_e)
        if not len(m):
            print("\n  %-46s no matching variants -> skipped" % tag); continue
        lead = int(m.loc[m.z_exp.abs().idxmax(), "position_e"])
        sel = m[(m.z_exp.abs() > MIN_ZEXP) & (m.maf > MIN_MAF) &
                (m.position_e.sub(lead).abs() <= WIN_HALF)].copy()
        if len(sel) < 3:
            print("\n  %-46s only %d qualifying variants -> skipped" % (tag, len(sel))); continue
        # ratio + delta-method SE
        sel["ratio"] = sel.beta_o / sel.beta_e
        sel["ratio_se"] = (sel.ratio.abs() *
                           np.sqrt((sel.se_o / sel.beta_o) ** 2 + (sel.se_e / sel.beta_e) ** 2))
        sel = sel[sel.ratio_se.notna() & (sel.ratio_se > 0)]
        w = 1.0 / sel.ratio_se ** 2
        rb = float((sel.ratio * w).sum() / w.sum())
        rb_se = float(np.sqrt(1 / w.sum()))
        Q = float((w * (sel.ratio - rb) ** 2).sum())
        dfree = len(sel) - 1
        Qp = float(stats.chi2.sf(Q, dfree)) if dfree > 0 else np.nan
        cv = float(sel.ratio.std(ddof=1) / abs(sel.ratio.mean())) if len(sel) > 1 else np.nan
        # ratio stability after dropping the single most extreme variant (sensitivity)
        s2 = sel.reindex(sel.ratio.sub(rb).abs().sort_values().index)[:-1] if len(sel) > 3 else sel
        cv2 = float(s2.ratio.std(ddof=1) / abs(s2.ratio.mean())) if len(s2) > 1 else np.nan

        print("\n  * %s" % tag)
        print("     window = exposure lead hg19:%d +/-%d kb" % (lead, WIN_HALF // 1000))
        print("     %d qualifying variants (|z_exp|>%.0f and MAF>%.2f); max exposure z %.1f"
              % (len(sel), MIN_ZEXP, MIN_MAF, sel.z_exp.abs().max()))
        print("     %-13s %-12s %9s %9s %10s %10s" % ("rsID", "hg19 pos", "b_IL6R", "b_outcome", "ratio", "ratio SE"))
        for _, r in sel.sort_values("z_exp", key=lambda s: -s.abs()).head(10).iterrows():
            print("     %-13s %-12d %+9.4f %+9.4f %10.4f %10.4f"
                  % (str(r.rsid_e)[:13], int(r.position_e), r.beta_e, r.beta_o, r.ratio, r.ratio_se))
        print("     weighted ratio = %.4f (SE %.4f)   coefficient of variation CV = %.1f%%   CV after dropping the single most extreme = %.1f%%"
              % (rb, rb_se, 100 * cv, 100 * cv2))
        print("     Cochran Q = %.1f (df=%d)  p = %.3g   -> %s"
              % (Q, dfree, Qp, ("ratios homogeneous (supports a shared causal variant)" if Qp > 0.05 else
                                "ratios heterogeneous (suggests more than one signal / horizontal pleiotropy present)")))
        rows.append(dict(tag=tag, n_snp=len(sel), ratio_weighted=rb, ratio_se=rb_se,
                         CV_pct=100 * cv, CV_trim1_pct=100 * cv2, CochranQ=Q, Q_df=dfree, Q_p=Qp,
                         max_z_exp=float(sel.z_exp.abs().max())))
        detail.append(sel.assign(tag=tag)[["tag", "rsid_e", "position_e", "ld_id", "beta_e", "se_e",
                                           "beta_o", "se_o", "z_exp", "maf", "ratio", "ratio_se"]])

    if rows:
        R = pd.DataFrame(rows)
        R.to_csv(os.path.join(COLOC, "effect_ratio_round2.csv"), index=False, encoding="utf-8-sig")
        pd.concat(detail).to_csv(os.path.join(COLOC, "effect_ratio_round2_detail.csv"),
                                 index=False, encoding="utf-8-sig")
        print("\n" + "=" * 118)
        print("Summary")
        print("=" * 118)
        print(R.to_string(index=False, float_format=lambda x: "%.4g" % x))
        print("\nsaved: effect_ratio_round2.csv / effect_ratio_round2_detail.csv")


if __name__ == "__main__":
    main()
