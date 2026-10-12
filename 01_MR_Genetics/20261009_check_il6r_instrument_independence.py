# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""* Self-check: are the 4 IL6R instruments **independent**? Is the MCH/MCV signal in the IL6R region **the same causal variant as the instrument loci**?

Two things (neither requires an OpenGWAS token):
  (1) **LD check**: the official FinnGen R12 LD matrix (`ld_matrix/finngen_r12_chr1_ld.tsv.gz`, 10 GB but **tabix-indexed**)
     long-table format: `#chr pos variant1 variant2 r r2` (LD measured in 520,210 FinnGen samples)
     -> what is the r^2 between each pair of my 4 instruments? If r^2 is high -> "4/4 concordant" is not 4 independent validations.
  (2) **Credible-set comparison (token-free colocalization)**: `lab_values/finemap/full/finemap/finngen_R12_<OMOPID>.FINEMAP.snp.bgz`
     columns include `prob`(PIP) and `cs`...`cs10` (credible-set index; -1 = not in any set)
     -> which SNPs are in the 95% credible set of MCH/MCV in the IL6R region? **Does it include the 4 instrument loci?**
       - if yes -> the evidence that the instrument loci and the outcome signal **share a causal variant is strengthened**
       - if no (the credible set is elsewhere) -> limitation 2 of section 15.8 is confirmed and the MR conclusion must be substantially weakened

* These two checks are key to the strength of the conclusion in section 15 and must be completed before finalizing.
"""
import os, sys, json
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\workbuddy_workspace\02-scripts")
import numpy as np
import pandas as pd
import finngen_tabix as ft

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
CACHE = os.path.join(DIR, "tabix_cache")
B = "https://storage.googleapis.com/finngen-public-data-r12/"
os.makedirs(CACHE, exist_ok=True)

IL6R_CH, IL6R_LO, IL6R_HI = "1", 154_300_000, 154_800_000
POS = {"rs4601580": 154421941, "rs4129267": 154453788, "rs6690230": 154460401, "rs7514452": 154465608}
TRAITS = [("3035941", "MCH"), ("3024731", "MCV"), ("3000963", "hemoglobin"),
          ("3022948", "serum iron"), ("3015399", "sTfR"), ("3006140", "total bilirubin")]


def ld_block():
    url = B + "ld_matrix/finngen_r12_chr1_ld.tsv.gz"
    print("=" * 110)
    print("(1) Instrument independence check (official FinnGen R12 LD matrix, 520,210 samples)")
    print("=" * 110)
    try:
        tbi = ft.read_tbi(url + ".tbi", os.path.join(CACHE, "finngen_chr1_ld.tbi"))
    except Exception as e:
        print("  failed to fetch tbi: %s" % str(e)[:120]); return None
    rows = ft.fetch_region(url, tbi, IL6R_CH, IL6R_LO, IL6R_HI)
    print("  Retrieved %d variant pairs in the IL6R region (r^2>0.01, window 3 Mb)" % len(rows))
    if not rows:
        print("  !! empty result"); return None
    d = pd.DataFrame([r.split("\t") for r in rows],
                     columns=["chrom", "pos", "v1", "v2", "r", "r2"][:len(rows[0].split("\t"))])
    for c in ["pos", "r", "r2"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    # variant format looks like 1:154421941:T:C
    def key(v):
        # * the format is `chr1_154420198_T_TA` (underscores), not `chr1:154420198:T:TA`
        try:
            return int(str(v).split("_")[1])
        except Exception:
            return None
    d["p1"] = d["v1"].map(key); d["p2"] = d["v2"].map(key)
    # consider only pairs between instrument loci
    P = list(POS.values())
    sub = d[d.p1.isin(P) & d.p2.isin(P)]
    print("\n  Pairwise LD (r^2) between the 4 instrument loci:")
    names = list(POS.keys())
    print("            " + " ".join("%-11s" % n for n in names))
    for n1, p1 in POS.items():
        line = "  %-10s" % n1
        for n2, p2 in POS.items():
            if p1 == p2:
                line += " %-11s" % "1.000"
            else:
                v = sub[((sub.p1 == p1) & (sub.p2 == p2)) | ((sub.p1 == p2) & (sub.p2 == p1))]
                line += " %-11s" % ("%.3f" % v.iloc[0]["r2"] if len(v) else "missing (<0.01)")
        print(line)
    # LD between the instrument loci and the strongest regional signal
    print("\n  LD between instrument loci x regional MCH/related top hits (see the position comparison in the next section)")
    out = sub[["p1", "p2", "r", "r2"]].rename(columns={"p1": "pos1", "p2": "pos2"})
    out.to_csv(os.path.join(DIR, "ld_il6r_tools.csv"), index=False, encoding="utf-8-sig")
    # LD between all regional SNPs and the instruments (to judge whether the regional top hit is in the same LD block as the instrument)
    full = d[["p1", "p2", "r2"]].copy()
    full = full[full.p1.isin(P) | full.p2.isin(P)]
    full.to_csv(os.path.join(DIR, "ld_il6r_tools_vs_region.csv"), index=False, encoding="utf-8-sig")
    print("  Saved ld_il6r_tools.csv / ld_il6r_tools_vs_region.csv (%d rows)" % len(full))
    return d


def finemap_block():
    print("\n\n" + "=" * 110)
    print("(2) Credible set of the IL6R region (FinnGen FINEMAP) -- are the instrument loci inside the 95% credible set of MCH/MCV?")
    print("=" * 110)
    allrows = []
    for omop, name in TRAITS:
        u = B + "lab_values/finemap/full/finemap/finngen_R12_%s.FINEMAP.snp.bgz" % omop
        try:
            tbi = ft.read_tbi(u + ".tbi", os.path.join(CACHE, "efm_%s.tbi" % omop))
        except Exception as e:
            print("\n  %-10s tbi failed %s" % (name, str(e)[:80])); continue
        # * in finemap the chromosome is in column 7 (tbi col_seq=6, col_beg=7); the column positions must be passed explicitly
        rows = ft.fetch_region(u, tbi, "chr1", IL6R_LO, IL6R_HI, cseq=6, cpos=7)
        if not rows:
            print("\n  %-10s no data in region" % name); continue
        cols = ["trait", "region", "v", "index", "rsid", "chromosome", "position", "allele1", "allele2",
                "maf", "beta", "se", "z", "prob", "log10bf", "mean", "sd", "mean_incl", "sd_incl", "p",
                "cs", "cs2", "cs3", "cs4", "cs5", "cs6", "cs7", "cs8", "cs9", "cs10"]
        d = pd.DataFrame([r.split("\t") for r in rows], columns=cols[:len(rows[0].split("\t"))])
        for c in ["position", "prob", "p", "maf", "beta"]:
            if c in d.columns:
                d[c] = pd.to_numeric(d[c], errors="coerce")
        d = d.dropna(subset=["position"])
        print("\n  -- %s --  %d SNPs in region" % (name, len(d)))
        # credible-set members (cs != -1 means it belongs to some credible set)
        cs_cols = [c for c in d.columns if c == "cs" or c.startswith("cs")]
        d["in_cs"] = (d[cs_cols].apply(pd.to_numeric, errors="coerce") != -1).any(axis=1)
        csd = d[d.in_cs]
        print("    Number of SNPs in the 95% credible set = %d" % len(csd))
        if len(csd):
            for _, r in csd.nsmallest(8, "p").iterrows():
                tag = ""
                for n, p in POS.items():
                    if abs(p - r["position"]) <= 2000:
                        tag = " * contains instrument " + n
                print("       chr1:%-11d PIP=%.3f p=%.2g beta=%+.4f maf=%.3f%s"
                      % (r["position"], r["prob"], r["p"], r["beta"], r["maf"], tag))
        # PIP of the instrument loci themselves
        print("    PIP of the instrument loci themselves: ", end="")
        for n, p in POS.items():
            h = d[(d["position"] >= p - 2) & (d["position"] <= p + 2)]
            print("%s=%.3f%s " % (n, h.iloc[0]["prob"] if len(h) else float("nan"),
                                  "(in CS)" if len(h) and h.iloc[0]["in_cs"] else ""), end="")
        print()
        allrows.append(d.assign(trait_name=name, omop=omop))
    if allrows:
        pd.concat(allrows, ignore_index=True).to_csv(
            os.path.join(DIR, "finemap_IL6R_region_all_traits.csv"), index=False, encoding="utf-8-sig")
        print("\n  Saved finemap_IL6R_region_all_traits.csv")


if __name__ == "__main__":
    ld_block()
    finemap_block()
