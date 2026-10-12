# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""Mechanism pillar summary: merge the panel results from each dataset into a single "cross-dataset directional concordance" master table

Input (same directory 03-data\MR_finngen\):
  mechanism_geo_results.csv          GSE202045 (IL-6 x human primary/HepG2/mouse) etc.
  mechanism_circadian_results.csv    GSE218242 (circadian time course) / GSE218243 (ARNTL KD)
  mechanism_replication_results.csv  GSE230325 (LPS dHepaRG) / GSE275342 (IL-6 Hep3B/Huh7)
Output:
  mechanism_summary_panel.csv   log2FC / direction / concordance for each gene x each contrast
  print the core gene matrix to the console (ready to copy straight into the manuscript)

* Lesson: the three scripts use different names for the "contrast" column (comparison / none / contrast), and geo_results has no dataset column,
   so after concat the pivot silently drops entire blocks because dataset=NaN -- we must **explicitly tag by source file**.
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd

OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"

CORE = ["UGT1A1", "UGT1A6", "UGT1A9", "SLCO1B1", "ABCC2", "ABCG2", "HMOX1", "BLVRA",
        "NR0B2", "HNF4A", "HAMP", "SLC40A1", "SOD2", "CAT", "CRP", "SOCS3", "SAA1",
        "CYP3A4", "CYP1A2", "ALB", "ARNTL", "NR1D1", "DBP", "PER2", "PER3"]

SOURCES = [
    ("mechanism_geo_results.csv", "GSE202045", "comparison"),
    ("mechanism_circadian_results.csv", "GSE218243", None),
    ("mechanism_replication_results.csv", "GSE230325", "contrast"),
]


def main():
    frames = []
    for fn, tag, ccol in SOURCES:
        p = os.path.join(OUT, fn)
        if not os.path.exists(p):
            print("  (missing) %s" % fn); continue
        d = pd.read_csv(p, encoding="utf-8-sig")
        # * src prefers the dataset column embedded in the file (one file may contain multiple datasets, e.g. replication contains GSE230325+GSE275342)
        if "dataset" in d.columns:
            d["src"] = d["dataset"].astype(str).str.replace(r"\s+.*$", "", regex=True).replace("", tag)
            d["src"] = d["src"].where(d["src"].astype(str).str.startswith("GSE"), tag)
        else:
            d["src"] = tag
        d["contrast_unified"] = d[ccol].astype(str) if ccol and ccol in d.columns else ""
        frames.append(d)
    if not frames:
        print("no input"); return
    D = pd.concat(frames, ignore_index=True)
    if "found" in D.columns:
        D = D[D["found"] == True]
    D.to_csv(os.path.join(OUT, "mechanism_summary_panel.csv"), index=False, encoding="utf-8-sig")
    print("saved mechanism_summary_panel.csv (%d rows)" % len(D))

    # keep only rows with log2FC and a readable contrast name (the circadian time-course log2FC column is empty, so amplitude is used instead)
    sub = D[D["log2FC"].notna()].copy()
    sub = sub[sub["contrast_unified"].astype(str).str.len() > 0]
    sub["col"] = sub["src"] + " | " + sub["contrast_unified"].str.slice(0, 34)
    piv = sub.pivot_table(index="gene", columns="col", values="log2FC", aggfunc="mean")
    pref = ["GSE202045", "GSE230325", "GSE275342"]
    order = [c for c in piv.columns if any(c.startswith(p) for p in pref)]
    piv = piv[order + [c for c in piv.columns if c not in order]]
    print("\n" + "=" * 118)
    print("Core gene x dataset log2FC direction matrix (%d contrast columns)" % len(piv.columns))
    print("=" * 118)
    for i, c in enumerate(piv.columns):
        print("  [%2d] %s" % (i + 1, c))
    print()
    print("  %-9s %s" % ("gene", "  ".join("[%2d]" % (i + 1) for i in range(len(piv.columns)))))
    for g in CORE:
        if g in piv.index:
            row = piv.loc[g]
            print("  %-9s %s" % (g, "  ".join("%+4.2f" % v if pd.notna(v) else "  -  " for v in row.values)))


if __name__ == "__main__":
    main()

