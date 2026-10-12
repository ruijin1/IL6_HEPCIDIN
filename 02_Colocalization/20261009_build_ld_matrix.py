# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""Build the LD correlation matrix for the IL6R region (FinnGen R12 official LD reference panel, 520,210 samples)

Why this is needed:
  `coloc.abf` **does not use LD** and allows only a single causal variant. In this project the exposure side
  (IL6R pQTL) has p-values saturated at 1e-265, z up to 38, and the exposure peak (chr1:154,426,264) is 20 kb
  away from the outcome peak (CRP: 154,395,212)
  -> coloc.abf is bound to call H3~=1, and **even the positive control (IL6R x CRP) does not pass**, so H3~=1
  cannot be read as "no colocalization".
  The correct approach is `coloc.susie`, which needs an LD matrix from the same source.
  The OpenGWAS LD endpoints are down, but the **FinnGen R12 official LD matrix ships with tabix, so it can be
  queried directly by region**.

Data format (as measured):
  `ld_matrix/finngen_r12_chr1_ld.tsv.gz` long table: `#chr  pos  variant1  variant2  r  r2`
  variants look like `chr1_154420198_T_TA`; the file is sorted by variant1 and **contains long-range pairs**
  (can span several Mb).

Outputs (coloc/ directory):
  ld_matrix_il6r.csv    dense correlation matrix (rows/columns = LD panel variant IDs, with header)
  ld_variant_map.csv    variant in each summary-stat file -> LD panel ID + direction sign (+1/-1)
"""
import os, sys, time
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd

sys.path.insert(0, r"D:\workbuddy_workspace\02-scripts")
import finngen_tabix as ft

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
COLOC = os.path.join(DIR, "coloc")
CACHE = os.path.join(DIR, "tabix_cache")
LDURL = "https://storage.googleapis.com/finngen-public-data-r12/ld_matrix/finngen_r12_chr1_ld.tsv.gz"

# susie window: covers the exposure peak (154.426M GRCh37) and the outcome peak (154.395-154.410M GRCh37)
WIN = (154_370_000, 154_480_000)          # GRCh37 coordinates (converted uniformly to GRCh38 below to query the LD panel)
OFF37TO38 = 27_524                        # measured constant offset for this interval (see 20261009_coloc_mch_mcv_jaundice_prep.py)

FILES = ["exp_IL6R_folkersen2018", "exp_IL6R_folkersen2020", "exp_IL6R_pietzner2020",
         "out_mch_finngen", "out_mcv_finngen", "out_bil_finngen_lab",
         "out_jaund_finngen", "out_crp_said2022"]

# * Which files have GRCh37 positions: everything from OpenGWAS is 37;
#   FinnGen sources (out_*_finngen) are already 38. Mixing them **silently fails all matching**
#   (hit this again this round).
GRCH37_FILES = {"exp_IL6R_folkersen2018", "exp_IL6R_folkersen2020", "exp_IL6R_pietzner2020",
                "out_crp_said2022", "out_bil_barton2021", "out_bil_ukb_neale", "out_bil_bbj",
                "out_bil_barton2021_SAAregion", "out_bil_ukb_neale_SAAregion",
                "exp_SAA1_pietzner2020"}


def main():
    print("=" * 108)
    print("(1) Collect the variants of each summary-stat file")
    print("=" * 108)
    recs = {}   # (file, id) -> (pos, ea, nea)
    for f in FILES:
        p = os.path.join(COLOC, f + ".csv")
        if not os.path.exists(p):
            print("   (missing) %s" % f); continue
        d = pd.read_csv(p, encoding="utf-8-sig")
        d = d[pd.to_numeric(d["position"], errors="coerce").notna()]
        for _, r in d.iterrows():
            recs[(f, str(r["rsid"]))] = (int(r["position"]), str(r["ea"]), str(r["nea"]))
        print("   %-26s %5d variants" % (f, len(d)))
    print("   total %d (file, variant) records" % len(recs))

    allpos = sorted({v[0] for v in recs.values()
                     if WIN[0] <= v[0] + (OFF37TO38 if False else 0) <= WIN[1]})
    print("   susie window (GRCh37 %d-%d)" % WIN)

    print("\n" + "=" * 108)
    print("(2) Fetch the FinnGen official LD (long table -> sparse -> dense matrix)")
    print("=" * 108)
    t0 = time.time()
    tbi = ft.read_tbi(LDURL + ".tbi", os.path.join(CACHE, "finngen_chr1_ld.tbi"))
    rows = ft.fetch_region(LDURL, tbi, "1", WIN[0] + OFF37TO38, WIN[1] + OFF37TO38)
    print("   fetched %d rows in %.0f s" % (len(rows), time.time() - t0))

    # build the LD panel variant dictionary: (pos, frozenset{a1,a2}) -> (ldid, a1, a2)
    ld_by_pos = {}
    R = {}
    for line in rows:
        f = line.split("\t")
        if len(f) < 6:
            continue
        try:
            pos1 = int(f[1]); r = float(f[4])
        except ValueError:
            continue
        v1 = f[2]
        # v1 looks like chr1_154420198_T_TA
        parts = v1.split("_")
        if len(parts) < 4:
            continue
        p1 = int(parts[1]); a1a, a1b = parts[2], "_".join(parts[3:])
        ld_by_pos[(p1, frozenset((a1a, a1b)))] = (v1, a1a, a1b)
        v2 = f[3]
        parts2 = v2.split("_")
        if len(parts2) < 4:
            continue
        p2 = int(parts2[1]); a2a, a2b = parts2[2], "_".join(parts2[3:])
        ld_by_pos[(p2, frozenset((a2a, a2b)))] = (v2, a2a, a2b)
        R[(v1, v2)] = r
    print("   LD panel variants (within window) %d, pairs %d" % (len(ld_by_pos), len(R)))

    print("\n" + "=" * 108)
    print("(3) Map summary-stat variants onto the LD panel (position + allele set), record the direction sign")
    print("=" * 108)
    map_rows, miss = [], {}
    for (f, vid), (pos, ea, nea) in recs.items():
        # * convert uniformly to GRCh38 to query the LD panel
        p38 = pos + OFF37TO38 if f in GRCH37_FILES else pos
        if not (WIN[0] + OFF37TO38 <= p38 <= WIN[1] + OFF37TO38):
            continue
        hit = ld_by_pos.get((p38, frozenset((ea, nea))))
        if hit is None:
            miss[f] = miss.get(f, 0) + 1
            continue
        ldid, a1, a2 = hit
        sign = 1 if (ea, nea) == (a1, a2) else (-1 if (ea, nea) == (a2, a1) else 0)
        if sign == 0:
            miss[f] = miss.get(f, 0) + 1
            continue
        map_rows.append(dict(file=f, id=vid, ld_id=ldid, pos37=(pos if f in GRCH37_FILES else pos - OFF37TO38),
                             pos38=p38, ea=ea, nea=nea, sign=sign))
    M = pd.DataFrame(map_rows)
    M.to_csv(os.path.join(COLOC, "ld_variant_map.csv"), index=False, encoding="utf-8-sig")
    print("   mapped successfully %d records" % len(M))
    for f in FILES:
        n = int((M.file == f).sum())
        print("     %-26s mapped %4d" % (f, n))

    print("\n" + "=" * 108)
    print("(4) Assemble the dense correlation matrix")
    print("=" * 108)
    ids = sorted(set(M.ld_id))
    idx = {v: i for i, v in enumerate(ids)}
    n = len(ids)
    S = np.eye(n, dtype=np.float32)
    got = 0
    for (v1, v2), r in R.items():
        i, j = idx.get(v1), idx.get(v2)
        if i is None or j is None:
            continue
        S[i, j] = r; S[j, i] = r
        got += 1
    print("   matrix %d x %d; filled %d off-diagonal elements" % (n, n, got))
    S = (S + S.T) / 2
    np.fill_diagonal(S, 1.0)
    # * FinnGen's LD long table **lists only part of the pairs** (measured 849,059 pairs within this window
    #   << all pairs), so the unlisted pairs are set to 0 -> the matrix is not positive definite and
    #   `susie_rss` will reject it.
    #   Handling: apply the **nearest positive-definite projection** (eigenvalue truncation at eps, then
    #   renormalize the diagonal to 1), and faithfully record the change in the minimum eigenvalue.
    w0 = np.linalg.eigvalsh(S.astype(np.float64))
    before = (w0[0], w0[-1])
    eps = 1e-4
    if w0[0] <= eps:
        w, V = np.linalg.eigh(S.astype(np.float64))
        w = np.clip(w, eps, None)
        S = (V * w) @ V.T
        d = np.sqrt(np.diag(S))
        S = S / np.outer(d, d)                      # turn it back into a correlation matrix
        np.fill_diagonal(S, 1.0)
        w1 = np.linalg.eigvalsh(S)
        print("   not positive definite -> nearest positive-definite projection applied: min eigenvalue %.4f -> %.6f, max %.2f -> %.2f (must be declared in the manuscript)"
              % (before[0], w1[0], before[1], w1[-1]))
    else:
        print("   matrix is positive definite (min eigenvalue %.4f)" % w0[0])
    S = S.astype(np.float32)
    out = pd.DataFrame(S, index=ids, columns=ids)
    fp = os.path.join(COLOC, "ld_matrix_il6r.csv")
    out.to_csv(fp, encoding="utf-8")
    print("   saved %s (%.1f MB)" % (fp, os.path.getsize(fp) / 2 ** 20))


if __name__ == "__main__":
    main()
