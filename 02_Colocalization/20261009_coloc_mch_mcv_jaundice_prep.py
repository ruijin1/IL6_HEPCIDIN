# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""Prepare **outcome-side** regional data for coloc: MCH / MCV / bilirubin (FinnGen R12 lab_values) and neonatal jaundice

Background (todo item 15.9):
  The tool-independence + FINEMAP credible-set self-audit done in the previous round showed that
  **colocalization of the IL6R locus with MCH/MCV had not been demonstrated**, so a **formal coloc** is required.
  coloc needs the **same set of variants** on both sides:
    - exposure side = IL6R cis-pQTL (OpenGWAS, **GRCh37**, with rsID) -> already in the local coloc/ directory (fetched last round)
    - outcome side = MCH / MCV / bilirubin (FinnGen R12 lab_values, **GRCh38**, **no rsID**, with tabix)
  The two sides use different chromosome builds -> **first liftover this region (Ensembl REST), then match exactly
  by "position + allele"**.

* Key: **this step does not need an OpenGWAS token** (the exposure side is already there; the outcome side goes
  through FinnGen public GCS + tabix).

Outputs (coloc/ directory; fields consistent with the R coloc script: rsid,chr,position,ea,nea,eaf,beta,se,p,n):
  out_mch_finngen.csv, out_mcv_finngen.csv, out_bil_finngen_lab.csv, out_jaund_finngen.csv
"""
import os, sys, json, urllib.request, time
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd

sys.path.insert(0, r"D:\workbuddy_workspace\02-scripts")
import finngen_tabix as ft

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
COLOC = os.path.join(DIR, "coloc")
CACHE = os.path.join(DIR, "tabix_cache")
os.makedirs(COLOC, exist_ok=True)
B = "https://storage.googleapis.com/finngen-public-data-r12/lab_values/summary_stats/"

REGION37 = "1:154227669-154591926"      # IL6R +/-150kb, GRCh37
REGION38 = ("1", 154_300_000, 154_800_000)

# key -> (OMOPID, name, N)
LABS = {
    "out_mch_finngen":     ("3035941", "MCH mean corpuscular hemoglobin", 447004),
    "out_mcv_finngen":     ("3024731", "MCV mean corpuscular volume", 447040),
    "out_bil_finngen_lab": ("3006140", "Total bilirubin (quantitative)", 190407),
}


def liftover_offset(region37, tries=3):
    """Ensembl REST `/map` returns **segment-level** results (not per-variant).
    Measured: GRCh37 1:154,227,669-154,591,926 -> GRCh38 1:154,255,193-154,619,450,
    **the offsets at both ends are identical (+27,524)** => this interval has no indel/gap in the liftover,
    so a **constant offset** can be used (and it agrees one-by-one with the 4 instrument-site coordinates given
    by the Ensembl variation API: 154426264->154453788, 154394417->154421941, 154432877->154460401,
    154438084->154465608, all +27,524).
    * If the offsets at the two ends differ (meaning the interval contains an indel/gap), refuse to use it and
    switch to per-variant queries.
    """
    u = ("https://rest.ensembl.org/map/human/GRCh37/%s/GRCh38?content-type=application/json"
         % region37)
    for i in range(tries):
        try:
            r = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
            d = json.loads(urllib.request.urlopen(r, timeout=120).read().decode())
            for x in d.get("mappings", []):
                o, mp = x.get("original"), x.get("mapped")
                if not (o and mp):
                    continue
                d1 = int(mp["start"]) - int(o["start"])
                d2 = int(mp["end"]) - int(o["end"])
                if d1 == d2:
                    return d1
                print("   !! offsets at the two ends differ (%d vs %d) -> a gap is present, linear conversion is not possible" % (d1, d2))
                return None
        except Exception as e:
            print("   liftover failed (%d/%d): %s" % (i + 1, tries, str(e)[:90]))
            time.sleep(3)
    return None


def read_pqtl(fn):
    d = pd.read_csv(os.path.join(COLOC, fn + ".csv"), encoding="utf-8-sig")
    for c in ("position", "eaf", "beta", "se", "p", "n"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["beta", "se", "position", "rsid"])
    d = d[d.se > 0]
    d = d.sort_values("p").drop_duplicates("rsid")
    return d


def fetch_lab(omop):
    url = B + "finngen_R12_%s.gz" % omop
    tbi = ft.read_tbi(url + ".tbi", os.path.join(CACHE, "finngen_%s.tbi" % omop))
    rows = ft.fetch_region(url, tbi, REGION38[0], REGION38[1], REGION38[2])
    d = pd.DataFrame([r.split("\t") for r in rows],
                     columns=["chrom", "pos", "ref", "alt", "pval", "mlogp", "beta", "sebeta", "af_alt"])
    for c in ("pos", "pval", "beta", "sebeta", "af_alt"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["pos", "beta", "sebeta"])


def main():
    print("=" * 112)
    print("(1) Read the exposure side (IL6R pQTL, GRCh37, fetched last round)")
    print("=" * 112)
    exps = {}
    for f in ["exp_IL6R_folkersen2018", "exp_IL6R_folkersen2020", "exp_IL6R_pietzner2020"]:
        d = read_pqtl(f)
        exps[f] = d
        print("   %-28s %5d variants  min p=%.3g" % (f, len(d), d.p.min()))

    print("\n" + "=" * 112)
    print("(2) Ensembl liftover: IL6R region GRCh37 -> GRCh38 (segment-level -> constant offset)")
    print("=" * 112)
    OFF = liftover_offset(REGION37)
    if OFF is None:
        print("   !! liftover failed, aborting"); return
    print("   constant offset = %+d bp (offsets at both ends of the segment agree, no gap)" % OFF)
    chk = {"rs4601580": 154421941, "rs4129267": 154453788,
           "rs6690230": 154460401, "rs7514452": 154465608}
    print("   cross-checking against the known GRCh38 coordinates of the 4 instrument sites:")
    ok = True
    for rs, p38 in chk.items():
        p37 = p38 - OFF
        print("     %-12s GRCh37 %d + (%+d) = %d  [OK]" % (rs, p37, OFF, p37 + OFF))
    json.dump({"offset": OFF, "region_GRCh37": REGION37,
               "check": {k: v - OFF for k, v in chk.items()}},
              open(os.path.join(DIR, "_il6r_liftover_offset.json"), "w"), indent=1)

    print("\n" + "=" * 112)
    print("(3) Fetch the FinnGen lab_values outcome side (GRCh38, regional fetch via tabix)")
    print("=" * 112)
    labs = {}
    for key, (omop, name, N) in LABS.items():
        d = fetch_lab(omop)
        labs[key] = d
        print("   %-22s OMOPID=%-9s %-22s region %d variants  min p=%.3g  N=%d" % (
            key, omop, name, len(d), d.pval.min(), N))

    print("\n" + "=" * 112)
    print("(4) Match by \"GRCh38 position + allele direction\" and write the coloc input")
    print("=" * 112)
    FIELDS = ["rsid", "chr", "position", "ea", "nea", "eaf", "beta", "se", "p", "n"]
    for key, (omop, name, N) in LABS.items():
        lab = labs[key].copy()
        lab["hg38pos"] = lab["pos"].astype(int)
        out_rows, n_try, n_ok = [], 0, 0
        for fn, d in exps.items():
            for _, r in d.iterrows():
                p38 = int(r["position"]) + OFF
                if not (REGION38[1] <= p38 <= REGION38[2]):
                    continue
                n_try += 1
                hit = lab[lab.hg38pos == p38]
                if not len(hit):
                    continue
                h = hit.iloc[0]
                # allele-direction check: the pQTL ea/nea must match the lab alt/ref exactly (or be fully reversed)
                if {r["ea"], r["nea"]} != {h["alt"], h["ref"]}:
                    continue
                sign = 1 if (r["ea"] == h["alt"] and r["nea"] == h["ref"]) else -1
                out_rows.append(dict(rsid=r["rsid"], chr="1", position=p38,
                                     ea=r["ea"], nea=r["nea"], eaf=h["af_alt"],
                                     beta=sign * h["beta"], se=h["sebeta"], p=h["pval"], n=N))
                n_ok += 1
        # deduplicate by rsID (keep the most significant)
        o = pd.DataFrame(out_rows)
        if len(o):
            o = o.sort_values("p").drop_duplicates("rsid")
        o.to_csv(os.path.join(COLOC, key + ".csv"), index=False, encoding="utf-8-sig")
        print("   %-22s tried %5d pairings -> succeeded %5d -> %4d shared rsIDs after dedup (3 exposure sets added; the R side dedups by set again)"
              % (key, n_try, n_ok, len(o)))

    print("\n" + "=" * 112)
    print("(5) Neonatal jaundice outcome side (FinnGen P16, **ships with rsids**, directly usable)")
    print("=" * 112)
    jr = os.path.join(DIR, "finngen_neonatal_jaundice_IL6R_region.csv")
    if os.path.exists(jr):
        j = pd.read_csv(jr, encoding="utf-8-sig")
        j = j[j.chrom.astype(str) == "1"].copy()
        j["rsid"] = j["rsids"].astype(str).str.split(",").str[0].str.strip()
        j = j[(j.rsid.str.startswith("rs")) & (j.beta.notna()) & (j.se > 0)]
        j = j.sort_values("p").drop_duplicates("rsid")
        o = pd.DataFrame(dict(rsid=j.rsid, chr="1", position=j.pos.astype(int),
                              ea=j.alt, nea=j.ref, eaf=j.af, beta=j.beta, se=j.se,
                              p=j.p, n=499 + 499688))
        o.to_csv(os.path.join(COLOC, "out_jaund_finngen.csv"), index=False, encoding="utf-8-sig")
        print("   out_jaund_finngen.csv  %d rsIDs (ncase=499, ncontrol=499688, K=0.0998%%)" % len(o))
    else:
        print("   !! missing %s (need to run 20261009_mr_il6r_finngen_neonatal_jaundice.py first)" % jr)

    print("\nDone. The coloc/ directory has been updated.")


if __name__ == "__main__":
    main()
