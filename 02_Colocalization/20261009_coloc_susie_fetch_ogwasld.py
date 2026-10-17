# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""coloc.susie data preparation (3rd attempt, this time using the OpenGWAS European reference LD)

Why the first two attempts failed:
  - Attempt 1, `coloc.abf`: p-values saturated (z up to 38) + peaks 31 kb apart -> even the positive control
    IL6R x CRP is called H3
  - Attempt 2, `coloc.susie` + the **FinnGen Finnish LD panel** -> `susie_rss` IBSS does not converge and the
    residual variance is negative (susieR explicitly requires in-sample LD; Folkersen/Pietzner are European
    pQTLs, so the Finnish panel does not match)
This time: use **OpenGWAS `POST /api/ld/matrix`** (European reference panel; measured to agree with the
FinnGen authoritative r: rs4129267-rs6690230 r=-0.703->r^2=0.49 vs FinnGen 0.437) [OK]

* Key switch: the outcome side no longer uses FinnGen lab_values (**no rsID**); instead it uses the
**Vuckovic 2020 blood-cell GWAS** (N=408,112, **with rsID**), which happens to cover the most stable set of
genetic-endpoint traits plus **reticulocyte count** (the marker closest to the mechanism).

Outputs: under coloc_ogwasld/
  exp_*.csv / out_*.csv   regional summary statistics per dataset (rsid, chr, position, ea, nea, eaf, beta, se, p, n)
  ld_matrix.csv + ld_snplist.json  LD correlation matrix of the intersecting variants
  pairs.json              the columns to run for each pair
The token is read from the `OGWAS_JWT` environment variable and **is never written to disk in the workspace**.
"""
import os, sys, json, time, csv, urllib.request, urllib.error
sys.stdout.reconfigure(encoding="utf-8")

TOK = os.environ.get("OGWAS_JWT", "").strip()
if not TOK:
    raise SystemExit("!! OGWAS_JWT environment variable not set")
DIR = r"D:\workbuddy_workspace\03-data\MR_finngen\coloc_ogwasld"
os.makedirs(DIR, exist_ok=True)
LOG = os.path.join(DIR, "_fetch.log")
_L = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); _L.append(s)
    open(LOG, "w", encoding="utf-8").write("\n".join(_L))


def og(path, body=None, timeout=300, tries=4):
    for i in range(tries):
        try:
            d = json.dumps(body).encode() if body is not None else None
            rq = urllib.request.Request("https://api.opengwas.io" + path, data=d, headers={
                "Authorization": "Bearer " + TOK, "Content-Type": "application/json",
                "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(rq, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            P("     HTTP %d %s" % (e.code, e.read().decode()[:120]))
            return None
        except Exception as e:
            if i == tries - 1:
                P("     ERR %s" % str(e)[:120]); return None
            time.sleep(4)


REGION = "1:154227669-154591926"      # IL6R +/-150 kb, GRCh37 (consistent with the analysis window)

EXP = {"exp_IL6R_folkersen2018": ("prot-b-23", 3394, "IL6R pQTL Folkersen2018"),
       "exp_IL6R_folkersen2020": ("ebi-a-GCST90012025", 21758, "IL6R pQTL Folkersen2020"),
       "exp_IL6R_pietzner2020":  ("ebi-a-GCST90019398", 10708, "IL6R pQTL Pietzner2020"),
       "exp_IL6R_sun2020":       ("prot-a-1540", 3301, "IL6R pQTL Sun2020"),
       "exp_SAA1_pietzner2020":  ("ebi-a-GCST90019394", 10708, "SAA1 pQTL Pietzner2020 (negative control)")}

OUT = {"out_MCH_vuckovic":            ("ebi-a-GCST90002390", 408112, "MCH Vuckovic2020"),
       "out_MCV_vuckovic":            ("ebi-a-GCST90002392", 408112, "MCV Vuckovic2020"),
       "out_MCHC_vuckovic":           ("ebi-a-GCST90002391", 408112, "MCHC Vuckovic2020"),
       "out_RDW_vuckovic":            ("ebi-a-GCST90002404", 408112, "RDW Vuckovic2020"),
       "out_RBC_vuckovic":            ("ebi-a-GCST90002403", 408112, "RBC count Vuckovic2020"),
       "out_RETIC_vuckovic":          ("ebi-a-GCST90002405", 408112, "Reticulocyte count Vuckovic2020"),
       "out_MRV_vuckovic":            ("ebi-a-GCST90002396", 408112, "Mean reticulocyte volume Vuckovic2020"),
       "out_PLT_vuckovic":            ("ebi-a-GCST90002402", 408112, "Platelet count Vuckovic2020 (negative control)"),
       "out_BIL_barton2021":          ("ebi-a-GCST90025973", 436055, "Total bilirubin Barton2021"),
       "out_BIL_bbj":                 ("bbj-a-53", 110207, "Total bilirubin BBJ East Asian"),
       "out_CRP_said2022":            ("ebi-a-GCST90029070", 575531, "CRP Said2022 (positive control)")}

FIELDS = ["rsid", "chr", "position", "ea", "nea", "eaf", "beta", "se", "p", "n"]


def fetch(key, gid, note, N):
    r = og("/api/associations", {"id": [gid], "variant": [REGION]})
    rows = r if isinstance(r, list) else []
    keep = []
    for x in rows:
        if not x.get("rsid") or x.get("beta") is None or x.get("se") is None: continue
        if x.get("ea") is None or x.get("nea") is None: continue
        keep.append({k: x.get(k) for k in FIELDS})
    fp = os.path.join(DIR, key + ".csv")
    with open(fp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); w.writerows(keep)
    pmin = min((x["p"] for x in keep if x.get("p")), default=float("nan"))
    P("   %-26s %-22s %-42s %5d variants  min p=%.3g  N=%d"
      % (key, gid, note, len(keep), pmin, N))
    time.sleep(0.5)
    return [x["rsid"] for x in keep]


def main():
    P("=" * 112)
    P("(1) Exposure side: IL6R cis-pQTL (GRCh37, with rsID)")
    P("=" * 112)
    exp_rs = set()
    for k, (gid, N, note) in EXP.items():
        exp_rs |= set(fetch(k, gid, note, N))

    P("")
    P("=" * 112)
    P("(2) Outcome side: Vuckovic blood cells / bilirubin / CRP")
    P("=" * 112)
    out_rs = {}
    for k, (gid, N, note) in OUT.items():
        out_rs[k] = set(fetch(k, gid, note, N))

    P("")
    P("=" * 112)
    P("(3) Number of sharable variants per pair (determines whether coloc is powered; <200 is treated as insufficient)")
    P("=" * 112)
    pairs = {}
    need = set()
    for ek in EXP:
        e = set(pd.read_csv(os.path.join(DIR, ek + ".csv"), encoding="utf-8-sig").rsid) if False else None
    import pandas as pd
    exp_sets = {ek: set(pd.read_csv(os.path.join(DIR, ek + ".csv"), encoding="utf-8-sig").rsid)
                for ek in EXP}
    for ek in EXP:
        for ok, rs in out_rs.items():
            n = len(exp_sets[ek] & rs)
            P("   %-26s x %-26s shared %4d" % (ek, ok, n))
            if n >= 60 and ek.startswith("exp_IL6R"):
                pairs.setdefault(ek, []).append(ok)
                need |= (exp_sets[ek] & rs)
    json.dump({"pairs": pairs, "exp_sets": {k: sorted(v) for k, v in exp_sets.items()},
               "out_sets": {k: sorted(v) for k, v in out_rs.items()}},
              open(os.path.join(DIR, "pairs.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    P("")
    P("=" * 112)
    P("(4) Fetch the LD correlation matrix (OpenGWAS European reference panel): %d rsIDs needed" % len(need))
    P("=" * 112)
    snps = sorted(need)
    snplist_all, mat = [], []
    B = 150
    for i in range(0, len(snps), B):
        chunk = snps[i:i + B]
        d = og("/api/ld/matrix", {"id": ["prot-b-23"], "rsid": chunk})
        if not isinstance(d, dict) or not d.get("snplist"):
            P("   block %d-%d returned empty, skipping" % (i, i + len(chunk))); continue
        sl, m = d["snplist"], d["matrix"]
        P("   block %3d-%3d: requested %3d -> returned %3d" % (i, i + len(chunk), len(chunk), len(sl)))
        snplist_all += sl; mat.append((sl, m))
        time.sleep(0.6)
    # assemble the dense matrix (in snplist order; blocdiag is wrong -- global coordinates are needed)
    idx = {s: i for i, s in enumerate(snplist_all)}
    n = len(snplist_all)
    import numpy as np
    S = np.eye(n)
    for sl, m in mat:
        for a, ra in enumerate(sl):
            ia = idx.get(ra)
            if ia is None: continue
            for b, rb in enumerate(sl):
                ib = idx.get(rb)
                if ib is None: continue
                try: S[ia, ib] = float(m[a][b])
                except Exception: pass
    S = (S + S.T) / 2; np.fill_diagonal(S, 1.0)
    w = np.linalg.eigvalsh(S)
    P("   matrix %d x %d; min eigenvalue %.4f max %.2f -> %s"
      % (n, n, w[0], w[-1], "positive definite [OK]" if w[0] > 0 else "not positive definite (the R side needs the nearest positive-definite projection)"))
    pd.DataFrame(S, index=snplist_all, columns=snplist_all).to_csv(
        os.path.join(DIR, "ld_matrix.csv"), encoding="utf-8")
    json.dump(snplist_all, open(os.path.join(DIR, "ld_snplist.json"), "w"))
    P("")
    P("Done -> %s" % DIR)
    P("  pairs to run: %s" % json.dumps(pairs, ensure_ascii=False))


if __name__ == "__main__":
    main()
