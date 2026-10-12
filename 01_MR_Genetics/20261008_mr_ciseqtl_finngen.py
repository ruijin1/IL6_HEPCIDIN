# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Design A, round 1: MR from four-gene cis-eQTL -> FinnGen R12 preterm birth/low birth weight

WARNING: the token is read from the environment variable OGWAS_JWT, and is **never written to this file or any log**.
   Usage: export OGWAS_JWT='...'  then run this script.

Design notes and empirical basis (2026-10-08):
  - Ensembl IDs must be confirmed by reverse lookup from the gene symbol -- writing them from memory is error-prone (already hit in this project: 158321 is actually AUTS2, 110448 is actually CD5)
  - OpenGWAS SERVICES__PIPELINE is in Error -> /api/tophits returns 502
    switched to the **interval form** of /api/associations ("11:118215059-118225876", hg19)
  - FinnGen summary statistics are hg38, but matching by **rsID** avoids the need for liftover
"""
import os, sys, gzip, json, time, argparse, urllib.request, urllib.error
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd

OG = "https://api.opengwas.io"
JWT = os.environ.get("OGWAS_JWT", "")

# -- Verified gene definitions (Ensembl resolved by symbol via MyGene; coordinates are hg19, matching the OpenGWAS dataset build) --
GENES = {
    "CD3G":  dict(ens="ENSG00000160654", hg19=("11", 118215059, 118225876)),
    "PRTN3": dict(ens="ENSG00000196415", hg19=("19", 840963, 848175)),
    "CD6":   dict(ens="ENSG00000013725", hg19=("11", 60739115, 60787849)),
    "THBS1": dict(ens="ENSG00000137801", hg19=("15", 39873280, 39891667)),
}

OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"
FINNGEN = os.path.join(OUT, "finngen_R12_preterm_lowbw.gz")
FINNGEN_EXPECT = 792_903_109   # official Content-Length; gzip raises EOFError when the download is incomplete
P_THRESH = 5e-8          # genome-wide significant
WINDOW = 1_000_000       # distance-pruning window (bp), an approximation of LD clumping


def og(path, body=None, timeout=90):
    if not JWT:
        raise SystemExit("!! OGWAS_JWT environment variable not set")
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(OG + path, data=data, headers={
        "Authorization": "Bearer " + JWT, "Content-Type": "application/json",
        "User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        print("   !! HTTP %s %s" % (e.code, e.read().decode()[:160]))
        return None


def step1_instruments():
    """Fetch cis-region eQTLs from OpenGWAS, filter p<5e-8, and prune by distance"""
    print("=" * 88)
    print("[1] Fetch cis-eQTL instruments (OpenGWAS interval query, bypassing the broken tophits)")
    print("=" * 88)
    ins = {}
    for gene, cfg in GENES.items():
        ch, s, e = cfg["hg19"]
        var = "%s:%d-%d" % (ch, s, e)
        r = og("/api/associations", {"id": ["eqtl-a-" + cfg["ens"]], "variant": [var]})
        rows = r if isinstance(r, list) else []
        sig = [x for x in rows if isinstance(x, dict) and (x.get("p") or 1) < P_THRESH]
        sig.sort(key=lambda x: x["p"])
        kept = []
        for x in sig:
            if all(abs(x["position"] - k["position"]) > WINDOW for k in kept):
                kept.append(x)
        ins[gene] = kept
        print("   %-6s %-18s region returned %-4d rows | p<5e-8: %-3d | instruments after pruning: %d"
              % (gene, cfg["ens"], len(rows), len(sig), len(kept)))
        for k in kept[:4]:
            print("        %-13s chr%s:%-10s p=%.2e beta=%+.4f se=%.4f ea=%s nea=%s"
                  % (k.get("rsid"), k.get("chr"), k.get("position"), k["p"],
                     k.get("beta", 0), k.get("se", 0), k.get("ea"), k.get("nea")))
    print()
    return ins


def step2_outcome(rsids):
    """Extract instrument effects from the local FinnGen summary statistics"""
    print("=" * 88)
    print("[2] Extract outcome from FinnGen R12 preterm birth/low birth weight (rsID matching)")
    print("=" * 88)
    if not os.path.exists(FINNGEN):
        print("   !! %s not found -- run the download first" % FINNGEN)
        return pd.DataFrame()
    got = os.path.getsize(FINNGEN)
    if got < FINNGEN_EXPECT:
        print("   !! Incomplete file: %s / %s bytes (%.1f%%) -- the download has not finished, re-run later"
              % (format(got, ","), format(FINNGEN_EXPECT, ","), 100.0 * got / FINNGEN_EXPECT))
        return pd.DataFrame()
    print("   File complete: %.1f MB" % (got / 2 ** 20))
    want = set(rsids)
    hits = []
    t0 = time.time()
    n = 0
    with gzip.open(FINNGEN, "rt", encoding="utf-8", errors="replace") as f:
        hdr = None
        for line in f:
            if line.startswith("##"):
                continue
            if hdr is None and line.startswith("#"):
                hdr = line.lstrip("#").rstrip("\n").split("\t")
                print("   header:", ",".join(hdr))
                continue
            n += 1
            f0 = line.split("\t", 8)
            if len(f0) < 8:
                continue
            rs = f0[3]
            if rs in want:
                hits.append(f0[:9])
    print("   Scanned %s rows, matched %d instrument SNPs (%.1fs)"
          % (format(n, ","), len(hits), time.time() - t0))
    if not hits:
        return pd.DataFrame()
    df = pd.DataFrame(hits, columns=(hdr[:9] if hdr else
                                     ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
                                      "pval", "beta", "sebeta"]))
    print()
    print(df.to_string(index=False, max_colwidth=26))
    print()
    return df


def step3_mr(ins, out):
    """Compute the MR estimate per gene (single SNP: Wald ratio; multiple SNPs: fixed-effect IVW)"""
    print("=" * 88)
    print("[3] MR estimate (exposure = gene expression, outcome = preterm birth/low birth weight)")
    print("=" * 88)
    if out.empty:
        print("   No outcome data, skipping")
        return None
    omap = {r["alt"] if "alt" in r else r.get("rsids"): r for _, r in out.iterrows()}
    omap = {r.iloc[3]: r for _, r in out.iterrows()}
    res = []
    for gene, kept in ins.items():
        bx, bxse, by, byse, used = [], [], [], [], []
        for k in kept:
            rs = k.get("rsid")
            if rs not in omap:
                continue
            o = omap[rs]
            by.append(float(o.iloc[7])); byse.append(float(o.iloc[8]))
            bx.append(k["beta"]); bxse.append(k["se"]); used.append(rs)
        if not bx:
            res.append(dict(gene=gene, n_snp=0, beta=np.nan, se=np.nan, p=np.nan,
                            note="no matching outcome SNP"))
            continue
        bx = np.array(bx); bxse = np.array(bxse)
        by = np.array(by); byse = np.array(byse)
        if len(bx) == 1:
            b = by[0] / bx[0]
            se = abs(byse[0] / bx[0])
            note = "Wald ratio"
        else:
            w = 1.0 / (byse ** 2 / bx ** 2)
            b = np.sum(w * (by / bx)) / np.sum(w)
            se = np.sqrt(1.0 / np.sum(w))
            note = "IVW (%d SNP)" % len(bx)
        z = b / se if se else np.nan
        from math import erfc
        p = erfc(abs(z) / np.sqrt(2)) if se else np.nan
        res.append(dict(gene=gene, n_snp=len(bx), beta=b, se=se, p=p, note=note))
    rdf = pd.DataFrame(res)
    print(rdf.to_string(index=False, float_format=lambda x: "%.4g" % x))
    print()
    print("   Interpretation: beta>0 means higher expression of this gene is associated with a higher risk of preterm birth/low birth weight;")
    print("                   p<0.05 is nominally significant; the single-SNP Wald ratio is not adjusted for multiple comparisons and must be stated in the main text.")
    p = os.path.join(OUT, "MR_ciseqtl_preterm_result.csv")
    rdf.to_csv(p, index=False, encoding="utf-8-sig")
    print("   Saved:", p)
    return rdf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-download-check", action="store_true")
    a = ap.parse_args()
    print("OpenGWAS cis-eQTL x FinnGen R12 preterm MR  %s" % time.strftime("%Y-%m-%d %H:%M"))
    print("token: %s" % ("read from environment variable (not persisted)" if JWT else "**missing**"))
    print()
    ins = step1_instruments()
    rsids = [k.get("rsid") for v in ins.values() for k in v]
    print("Total instrument SNPs: %d (before dedup)\n" % len(rsids))
    out = step2_outcome(set(rsids))
    step3_mr(ins, out)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "MR_instruments.json"), "w", encoding="utf-8") as f:
        json.dump({g: [{k: v for k, v in x.items()} for x in lst] for g, lst in ins.items()},
                  f, ensure_ascii=False, indent=1)
    print("\nDone [OK]")


if __name__ == "__main__":
    main()
