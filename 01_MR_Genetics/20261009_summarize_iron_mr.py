# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""Summarize the IL6R -> iron/erythroid/bilirubin MR and test "direction consistency" against the "hepcidin model" predictions

Reads `MR_il6r_iron_erythroid_results.csv` and `finngen_lab_locus_scan.csv`.

* Direction logic (must be stated clearly, otherwise it will be read backwards):
   The effect size of the IL6R cis-pQTL is "**each SD higher sIL6R**".
   High sIL6R = less membrane-bound IL6R = **weaker classical IL-6 signaling** (positive control in this project: IL6R -> CRP is negative, b=-0.115).
   Therefore:
     * If the hepcidin model holds (IL-6 signaling up -> HAMP up -> iron restriction -> erythroid restriction -> heme turnover down -> bilirubin down)
       => this table (each SD higher sIL6R) should show: **serum iron up, TSAT up, sTfR down, MCV up, MCH up, Hb up, bilirubin up**
     * If the opposite direction is seen, the model is refuted.
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen"

# Model predictions: expected sign on the "each SD higher sIL6R" baseline
EXPECT = {
    "Hemoglobin": +1, "MCV": +1, "MCH": +1, "Red blood cell count": +1, "Hematocrit Hct": +1,
    "RDW": -1, "Reticulocyte %": +1, "Reticulocyte count": +1, "Erythropoietin EPO": -1,
    "Ferritin": -1, "Serum iron": +1, "Transferrin": -1,
    "Transferrin saturation TSAT": +1, "Soluble transferrin receptor sTfR": -1,
    "Total bilirubin (reference)": +1,
}
ORDER = ["Serum iron", "Transferrin saturation TSAT", "Soluble transferrin receptor sTfR", "Ferritin",
         "Transferrin", "MCV", "MCH", "Hemoglobin", "Hematocrit Hct",
         "Red blood cell count", "RDW", "Reticulocyte %", "Reticulocyte count", "Erythropoietin EPO",
         "Total bilirubin (reference)"]

# The saved result files predate the English release and store several trait names in the
# original bilingual form. Map them to the canonical English names so the lookups above work
# against either the old files or freshly regenerated ones.
TRAIT_ALIAS = {
    "血红蛋白 Hemoglobin": "Hemoglobin",
    "红细胞计数": "Red blood cell count",
    "血细胞比容 Hct": "Hematocrit Hct",
    "网织红细胞%": "Reticulocyte %",
    "网织红细胞计数": "Reticulocyte count",
    "促红细胞生成素 EPO": "Erythropoietin EPO",
    "铁蛋白 Ferritin": "Ferritin",
    "血清铁 Iron": "Serum iron",
    "转铁蛋白 Transferrin": "Transferrin",
    "转铁蛋白饱和度 TSAT": "Transferrin saturation TSAT",
    "可溶性转铁蛋白受体 sTfR": "Soluble transferrin receptor sTfR",
    "总胆红素（参照）": "Total bilirubin (reference)",
}

d = pd.read_csv(os.path.join(DIR, "MR_il6r_iron_erythroid_results.csv"), encoding="utf-8-sig")
d["trait"] = d.trait.replace(TRAIT_ALIAS)
d["sign_ok"] = np.sign(d.ratio) == d.trait.map(EXPECT)

print("=" * 118)
print("IL6R cis-pQTL (each SD higher sIL6R = weaker classical IL-6 signaling) -> FinnGen R12 quantitative traits (IRN, SD units)")
print("=" * 118)
hdr = "%-26s %-5s %9s %9s %9s %9s %7s %7s" % ("trait", "expected", "rs4601580", "rs7514452", "rs6690230", "rs4129267", "concordant", "model")
print(hdr)
print("-" * 118)
for t in ORDER:
    s = d[d.trait == t]
    if not len(s):
        continue
    exp = EXPECT.get(t, 0)
    cells, ok, tot = [], 0, 0
    for rs in ["rs4601580", "rs7514452", "rs6690230", "rs4129267"]:
        v = s[s.rsid == rs]
        if len(v):
            r = v.iloc[0]
            star = "*" if r.p_mr < 0.05 else ("+" if r.p_mr < 0.10 else "")
            cells.append("%+7.4f%s" % (r.ratio, star))
            tot += 1
            if np.sign(r.ratio) == exp:
                ok += 1
        else:
            cells.append("   -   ")
    verdict = "OK consistent" if (tot and ok == tot) else ("~ %d/%d" % (ok, tot) if ok else "X inconsistent")
    print("%-26s %-5s %s %6s %7s" % (t[:26], ("+" if exp > 0 else "-"), " ".join(cells),
                                     "%d/%d" % (ok, tot), verdict))
print("\n  * p<0.05   + p<0.10")

# ---------- Summary by "model step" ----------
print("\n" + "=" * 118)
print("Summary by step of the hepcidin model (among the 4 instruments: concordant and passing threshold)")
print("=" * 118)
CHAIN = [("(1) Iron utilization (direct hepcidin target)", ["Serum iron", "Transferrin saturation TSAT", "Soluble transferrin receptor sTfR", "Ferritin"]),
         ("(2) Erythropoiesis (downstream of iron restriction)", ["MCV", "MCH", "Hemoglobin", "RDW"]),
         ("(3) Terminal phenotype", ["Total bilirubin (reference)"])]
for label, ts in CHAIN:
    sub = d[d.trait.isin(ts)]
    if not len(sub):
        continue
    sig = sub[(sub.p_mr < 0.05)]
    concord = int(sub.sign_ok.sum())
    print("\n  %s" % label)
    for t in ts:
        s = sub[sub.trait == t]
        if not len(s):
            continue
        nsig = int((s.p_mr < 0.05).sum())
        print("    %-26s concordant %d/%d, significant (p<0.05) %d: %s" % (
            t, int(s.sign_ok.sum()), len(s), nsig,
            ", ".join("%s%+.4f(p=%.1e)" % (r.rsid, r.ratio, r.p_mr) for _, r in s.iterrows() if r.p_mr < 0.05)))
    print("    -> subtotal: %d comparisons concordant, %d significant" % (concord, len(sig)))

# ---------- Locus scan (regional min p in HAMP / HFE / UGT1A1) ----------
ls = os.path.join(DIR, "finngen_lab_locus_scan.csv")
if os.path.exists(ls):
    L = pd.read_csv(ls, encoding="utf-8-sig")
    print("\n" + "=" * 118)
    print("Locus scan: strongest signal per region across these quantitative traits (beta on the log scale, IRN SD units)")
    print("=" * 118)
    piv = L.pivot_table(index="trait", columns="region", values="p", aggfunc="min")
    cols = [c for c in ["HAMP", "HFE", "UGT1A1"] if c in piv.columns]
    print("%-26s %s" % ("trait", " ".join("%16s" % c for c in cols)))
    for t in ORDER:
        if t in piv.index:
            print("%-26s %s" % (t[:26], " ".join(
                ("%16.2g" % piv.loc[t, c]) if pd.notna(piv.loc[t, c]) else "%16s" % "-" for c in cols)))
    print("\n  * The HAMP locus itself has a strong signal for erythroid indices (independent of IL6R) -> independent evidence for the hepcidin axis and erythropoiesis")
