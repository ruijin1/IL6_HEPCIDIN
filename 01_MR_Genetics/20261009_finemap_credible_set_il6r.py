# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""(3) FinnGen FINEMAP credible set: does the MCH/MCV signal in the IL6R region **share a causal variant** with the 4 pQTL instrument loci?

(token-free "colocalization-lite". True coloc needs region-level data on the pQTL side -> requires an OpenGWAS JWT.)

* Data: `lab_values/finemap/full/finemap/finngen_R12_<OMOPID>.FINEMAP.snp.bgz` (1-3 MB)
  columns: trait region v index rsid chromosome position allele1 allele2 maf beta se z **prob(PIP)** log10bf
      mean sd mean_incl sd_incl p **cs cs2...cs10** (credible-set index; -1 = not in any set)
* Empirically observed: the `.tbi` of these small files is **not fully consistent** with the actual file (bin-based reads under-cover and large-range queries return 0 rows)
  -> download the whole file instead (see `finngen_tabix.fetch_all`).
"""
import os, sys, time, json
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\workbuddy_workspace\02-scripts")
import numpy as np
import pandas as pd
import finngen_tabix as ft

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"
B = "https://storage.googleapis.com/finngen-public-data-r12/lab_values/finemap/full/finemap/"
LO, HI = 154_300_000, 154_800_000
POS = json.load(open(os.path.join(DIR, "_il6r_snp_coords.json"), encoding="utf-8"))
TOOLS = {k: v["pos"] for k, v in POS.items() if k in
         ("rs4601580", "rs4129267", "rs6690230", "rs7514452") and v}
TRAITS = [("3035941", "MCH"), ("3024731", "MCV"), ("3000963", "hemoglobin"),
          ("3022948", "serum iron"), ("3015399", "sTfR"), ("3001122", "ferritin"),
          ("3006140", "total bilirubin")]
CACHE = os.path.join(DIR, "tabix_cache")
os.makedirs(CACHE, exist_ok=True)


def load(omop, name):
    """* Fetch by tabix region (column positions derived automatically from the tbi) -- FINEMAP files can reach 400 MB and cannot be downloaded whole.
       Empirically: `fetch_region`'s cseq/cpos default to being derived from the tbi's col_seq/col_beg (1-based -> -1),
       which is exactly the root cause of the earlier "no data in region"."""
    url = B + "finngen_R12_%s.FINEMAP.snp.bgz" % omop
    tbi = ft.read_tbi(url + ".tbi", os.path.join(CACHE, "efm_%s.tbi" % omop))
    t0 = time.time()
    rows = ft.fetch_region(url, tbi, "chr1", LO, HI)
    print("    %-8s retrieved %d rows in region, %.0f s" % (name, len(rows), time.time() - t0))
    if not rows:
        return None
    cols = rows[0].split("\t") if rows[0].startswith("trait") else \
        ["trait", "region", "v", "index", "rsid", "chromosome", "position", "allele1", "allele2",
         "maf", "beta", "se", "z", "prob", "log10bf", "mean", "sd", "mean_incl", "sd_incl", "p",
         "cs", "cs2", "cs3", "cs4", "cs5", "cs6", "cs7", "cs8", "cs9", "cs10"][:len(rows[0].split("\t"))]
    d = pd.DataFrame([l.split("\t") for l in rows if l.count("\t") >= len(cols) - 1], columns=cols)
    for c in ["position", "prob", "p", "maf", "beta", "se"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["position"])


def main():
    print("=" * 112)
    print("(3) FinnGen FINEMAP credible set (IL6R region chr1:%d-%d)" % (LO, HI))
    print("   instrument loci GRCh38:", ", ".join("%s=%d" % (k, v) for k, v in TOOLS.items()))
    print("=" * 112)
    rows = []
    for omop, name in TRAITS:
        try:
            d = load(omop, name)
        except Exception as e:
            print("\n  %-10s failed %s" % (name, str(e)[:90])); continue
        if d is None:
            print("\n  %-10s no data in region (this segment was not fine-mapped)" % name); continue
        sub = d[(d.position >= LO) & (d.position <= HI)].copy()
        if not len(sub):
            print("\n  %-10s this region is not within its fine-mapped segment (FINEMAP only covers segments with signal)" % name); continue
        cs_cols = [c for c in d.columns if c == "cs" or (c.startswith("cs") and c[2:].isdigit())]
        for c in cs_cols:
            sub[c] = pd.to_numeric(sub[c], errors="coerce")
        sub["in_cs"] = (sub[cs_cols] != -1).any(axis=1)
        top = sub.nsmallest(5, "p")
        print("\n  -- %-8s %d SNPs in region, %d in the 95%% credible set; min p=%.3g --"
              % (name, len(sub), int(sub.in_cs.sum()), sub.p.min()))
        print("     %-14s %8s %10s %9s %8s %6s" % ("position", "PIP", "p", "beta", "maf", "in CS"))
        for _, r in top.iterrows():
            tag = ""
            for k, v in TOOLS.items():
                if abs(v - r["position"]) <= 2000:
                    tag = "(4) " + k
            print("     chr1:%-9d %8.3f %10.2g %+9.4f %8.3f %6s %s"
                  % (r["position"], r["prob"], r["p"], r["beta"], r["maf"], "yes" if r["in_cs"] else "no", tag))
        print("     The 4 instrument loci themselves: ", end="")
        for k, v in TOOLS.items():
            h = sub[(sub.position >= v - 2) & (sub.position <= v + 2)]
            if len(h):
                print("%s(PIP=%.3f%s) " % (k, h.iloc[0]["prob"], " -in CS" if h.iloc[0]["in_cs"] else ""), end="")
            else:
                print("%s(missing) " % k, end="")
        print()
        # whether the instrument loci fall inside the credible set (+/-2 kb)
        hit = [k for k, v in TOOLS.items()
               if len(sub[(sub.position >= v - 2) & (sub.position <= v + 2)]) and
               bool(sub[(sub.position >= v - 2) & (sub.position <= v + 2)].iloc[0]["in_cs"])]
        print("     -> instrument loci inside the 95%% credible set: %s" % (", ".join(hit) if hit else "**none**"))
        sub["trait_name"] = name
        rows.append(sub)
    if rows:
        pd.concat(rows, ignore_index=True).to_csv(
            os.path.join(DIR, "finemap_IL6R_region_all_traits.csv"), index=False, encoding="utf-8-sig")
        print("\nSaved finemap_IL6R_region_all_traits.csv")


if __name__ == "__main__":
    main()
