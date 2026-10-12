# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""R3 · GSE102006 (HepG2 +/- cytokine mix, 6/12/24h, n=252, Agilent GPL6480)

* Use only the **solvent control arm** (DMSO / TL_DMSO), avoiding the pharmacological effects of the 6 test compounds:
   `HepG2_DMSO_with(out)_cytokines_<t>h_exp<n>`
   `HepG2_TL_DMSO_with(out)_cytokines_<t>h_exp<n>`
   contrast = with cytokines vs without cytokines, run separately at each of 6/12/24h.
* Cytokine mix composition (GEO original text): TNF 100 + IFNγ 100 + IL-1α 20 + **IL-6 20** ng/mL + LPS 10 µg/mL
   -- **not IL-6 alone**, so it can only serve as weak support for "inflammation" and must not be treated as IL-6-specific evidence (already qualified in the report).
* probe -> symbol uses the GENE_SYMBOL column of GPL6480.annot.gz; for multiple probes of the same name, take the one with the **highest mean expression**.
* The values are already on a log scale (Agilent processed signal), so no further log transform.
"""
import os, sys, gzip, re
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"

PANEL = [
    ("UGT1A1", "*bilirubin conjugation"), ("UGT1A6", "UGT1A6"), ("UGT1A9", "UGT1A9"), ("UGT2B7", "UGT2B7"),
    ("HMOX1", "heme oxygenase 1"), ("HMOX2", "heme oxygenase 2"), ("BLVRA", "biliverdin reductase A"), ("BLVRB", "biliverdin reductase B"),
    ("ALAS1", "heme synthesis"), ("FECH", "ferrochelatase"),
    ("SLCO1B1", "hepatic uptake OATP1B1"), ("SLCO1B3", "hepatic uptake OATP1B3"),
    ("ABCC2", "efflux MRP2"), ("ABCC3", "efflux MRP3"), ("ABCG2", "efflux BCRP"),
    ("NR0B2", "SHP (represses UGT1A1)"), ("HNF4A", "HNF4α (UGT1A1 transcription factor)"), ("NR1H4", "FXR"),
    ("HAMP", "*hepcidin"), ("SLC40A1", "FPN1"), ("TFRC", "TFRC"), ("FTL", "ferritin light chain"), ("FTH1", "ferritin heavy chain"),
    ("IL6", "*IL-6"), ("IL6R", "IL-6 receptor"), ("IL6ST", "gp130"), ("SOCS3", "*IL-6/STAT3 target"),
    ("STAT3", "STAT3"), ("CRP", "*CRP"), ("SAA1", "*SAA1"), ("LBP", "*LBP"),
    ("ORM1", "α1-acid glycoprotein"), ("SERPINA1", "α1-antitrypsin"),
    ("FGA", "fibrinogen α"), ("FGB", "fibrinogen β"), ("FGG", "fibrinogen γ"), ("HP", "haptoglobin"), ("C3", "complement C3"),
    ("CYP3A4", "CYP3A4 (liver homeostasis)"), ("CYP1A2", "CYP1A2 (liver homeostasis)"), ("CYP2C9", "CYP2C9"), ("ALB", "albumin"),
    ("APOB", "ApoB"), ("TTR", "TTR"), ("NR1D1", "REV-ERBα"), ("ARNTL", "BMAL1"), ("DBP", "DBP"),
    ("NFKB1", "NF-κB1"), ("RELA", "RELA"), ("TNFAIP3", "A20"), ("IL1B", "IL-1β"),
    ("ACTB", "β-actin (housekeeping)"), ("GAPDH", "GAPDH (housekeeping)"),
]
ACUTE = ["CRP", "SAA1", "LBP", "ORM1", "SERPINA1", "FGA", "FGB", "FGG", "HP", "C3"]


def read_symbols(fp):
    """GPL6480.annot.gz -> {probe: symbol}
    * The header is actually `ID` / `Gene symbol` (lowercase s, with a space), not `GENE_SYMBOL`
      -- testing with `"GENE_SYMBOL" in h.upper()` misses it (`GENE SYMBOL` contains a space).
    """
    m = {}
    ip = isym = None
    with gzip.open(fp, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(("^", "!", "#")):
                continue
            flds = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if ip is None:
                hdr = [h.upper().replace(" ", "_").replace("_REF", "") for h in flds]
                ip = next((i for i, h in enumerate(hdr) if h in ("ID", "ID_REF")), 0)
                isym = next((i for i, h in enumerate(hdr)
                             if h in ("GENE_SYMBOL", "GENE_SYMBOL_S", "SYMBOL")), None)
                if isym is None:
                    print("  !! annotation table has no gene symbol column, header:", flds[:12]); return {}
                continue
            if len(flds) <= max(ip, isym):
                continue
            sym = flds[isym]
            if sym and sym not in ("", "---"):
                m[flds[ip]] = sym.split("///")[0].strip()
    return m


def main():
    smx = os.path.join(DIR, "GSE102006_series_matrix.txt.gz")
    ann = os.path.join(DIR, "GPL6480.annot.gz")
    if not os.path.exists(smx) or not os.path.exists(ann):
        print("series matrix or GPL6480 annotation missing, skipping"); return
    pmap = read_symbols(ann)
    print("GPL6480 probe->symbol mappings: %d" % len(pmap))

    # ---- read sample titles + expression table ----
    titles, accs = [], []
    data_started = False
    rows = {}
    with gzip.open(smx, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                data_started = True
                next(f)  # header row
                continue
            elif data_started:
                if line.startswith("!series_matrix_table_end"):
                    break
                flds = line.rstrip("\n").split("\t")
                rows[flds[0].strip('"')] = flds[1:]
    print("samples %d; probes %d" % (len(accs), len(rows)))
    if not accs or not rows:
        print("!! parsing failed"); return
    # * must use orient="index": otherwise pandas treats the dict keys as **column names**, producing a transposed 252x28796 matrix
    if len(accs) != len(next(iter(rows.values()))):
        print("!! sample count (%d) and data column count (%d) differ, aborting" % (len(accs), len(next(iter(rows.values()))))); return
    X = pd.DataFrame.from_dict(rows, orient="index")
    X.columns = accs
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X[(X.notna().sum(axis=1) > 0)]
    X = X.fillna(X.mean(axis=1))
    print("expression matrix %d probes x %d samples" % X.shape)
    # probe -> symbol (for the same name, take the highest mean expression)
    mean_expr = X.mean(axis=1)
    best = {}
    for p in X.index:
        s = pmap.get(p)
        if not s:
            continue
        if s not in best or mean_expr[p] > mean_expr[best[s]]:
            best[s] = p
    print("usable gene symbols: %d" % len(best))
    symX = X.loc[list(best.values())]
    symX.index = list(best.keys())

    # ---- grouping: use only the solvent control arm ----
    meta = pd.DataFrame({"acc": accs, "title": titles})
    meta["cyto"] = meta.title.str.contains("with_cytokines", case=False)
    meta["nocyto"] = meta.title.str.contains("without_cytokines", case=False)
    meta["solvent"] = meta.title.str.contains(r"_DMSO_|_TL_DMSO_", case=False)
    meta["time"] = meta.title.str.extract(r"_(\d+)h_")[0]
    sub = meta[meta.solvent & (meta.cyto | meta.nocyto)].copy()
    print("\nsolvent control arm samples %d (of %d)" % (len(sub), len(meta)))
    print(sub.groupby(["time", "cyto"]).size().to_string())

    rows_out = []
    for t in sorted(sub.time.dropna().unique(), key=lambda x: int(x)):
        case = list(sub[(sub.time == t) & sub.cyto].acc)
        ctrl = list(sub[(sub.time == t) & sub.nocyto].acc)
        if len(case) < 2 or len(ctrl) < 2:
            print("\n[%sh] too few samples (%d/%d), skipping" % (t, len(case), len(ctrl))); continue
        npair = "%d/%d" % (len(case) * len(ctrl), len(case) * len(ctrl))
        print("\n" + "=" * 108)
        print("GSE102006 HepG2 cytokine mix vs solvent   %sh   (n=%d vs %d)" % (t, len(case), len(ctrl)))
        print("=" * 108)
        recs = []
        for sym, note in PANEL:
            if sym not in symX.index:
                rows_out.append(dict(dataset="GSE102006 HepG2 cytokines", contrast="%sh" % t,
                                     gene=sym, note=note, found=False)); continue
            a = symX.loc[sym, case].astype(float).values
            b = symX.loc[sym, ctrl].astype(float).values
            fc = float(np.mean(a) - np.mean(b))
            pairs = [(x - y) for x in a for y in b]
            cons = sum(1 for v in pairs if (v > 0) == (fc > 0)) if fc != 0 else 0
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) > 1 and len(b) > 1 else np.nan
            rec = dict(dataset="GSE102006 HepG2 cytokines", contrast="%sh" % t, gene=sym, note=note, found=True,
                       mean_case=float(np.mean(a)), mean_ctrl=float(np.mean(b)), log2FC=fc,
                       pair_consist="%d/%d" % (cons, len(pairs)), p=p)
            recs.append(rec); rows_out.append(rec)
        recs.sort(key=lambda r: -abs(r["log2FC"]))
        print("%-9s %8s %8s %8s %8s %9s  %s" % ("gene", "cytokine", "solvent", "log2FC", "dir_cons", "p", "note"))
        print("-" * 108)
        for x in recs:
            print("%-9s %8.2f %8.2f %+8.2f %8s %9.2g  %s" % (x["gene"], x["mean_case"], x["mean_ctrl"],
                                                             x["log2FC"], x["pair_consist"], x["p"], str(x["note"])[:20]))
        ac = [x for x in recs if x["gene"] in ACUTE]
        up = [x for x in ac if x["log2FC"] > 0.3 and x["pair_consist"] == npair]
        print("\n  -- hard validity gate (acute-phase genes with %s concordance and log2FC>0.3) -- %d/%d passed: %s" % (
            npair, len(up), len(ac), ", ".join("%s%+.2f" % (x["gene"], x["log2FC"]) for x in up)))
        print("     verdict: %s" % ("✓ model active, results are trustworthy" if len(up) >= 3 else "✗ model not active -- results at this time point must not be trusted"))
        key = {x["gene"]: x for x in recs}
        for g in ["UGT1A1", "UGT1A6", "UGT1A9", "HMOX1", "HAMP", "CYP3A4", "CYP1A2", "ALB", "IL6", "SOCS3"]:
            if g in key:
                print("        %-8s log2FC=%+.2f  %s" % (g, key[g]["log2FC"], key[g]["pair_consist"]))

    if rows_out:
        d = pd.DataFrame(rows_out)
        p = os.path.join(OUT, "mechanism_GSE102006_results.csv")
        d.to_csv(p, index=False, encoding="utf-8-sig")
        print("\nresults saved: %s (%d rows)" % (p, len(d)))


if __name__ == "__main__":
    main()
