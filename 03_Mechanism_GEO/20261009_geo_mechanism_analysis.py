# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Mechanism pillar: effect of IL-6 / infection on bilirubin metabolism pathway gene expression

Datasets (sample grouping is always parsed from GEO metadata / sample names, never guessed):
  T1  GSE202045  IL-6-stimulated hepatocytes (treatment: IL6 stimulated vs No treatment; time 0/6/24h; 3 replicates)
      · human primary hepatocytes (hPH, r1-r3) = primary model
      · HepG2 (r1-r2) = replication model
      · mouse primary hepatocytes = cross-species replication
  T2  GSE25504   neonatal sepsis whole blood mRNA (case vs control)
  T3  GSE218242/218243  primary hepatocytes: circadian rhythm x inflammation x drug metabolism (pending annotation parsing)

Methods: CPM -> log2 (counts) or log2(TPM); Mann-Whitney U + per-donor directional concordance
* Mandatory positive control: are the known IL-6 / acute-phase target genes (SOCS3/CRP/SAA1/LBP/ORM1/FGA/FGB/SERPINA1) induced?
  -- if they are negative, the stimulation model or the grouping parsing is wrong, and none of the other results can be used
"""
import os, gzip, re, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"

PANEL = [
    ("HMOX1", "heme oxygenase 1 (produces bilirubin)"), ("HMOX2", "heme oxygenase 2"),
    ("BLVRA", "biliverdin reductase A"), ("BLVRB", "biliverdin reductase B"),
    ("ALAS1", "rate-limiting enzyme of heme synthesis"), ("FECH", "ferrochelatase"),
    ("UGT1A1", "*bilirubin conjugation (determines serum bilirubin)"), ("UGT1A6", "UGT1A6"), ("UGT1A9", "UGT1A9"),
    ("SLCO1B1", "hepatic uptake"), ("SLCO1B3", "hepatic uptake"), ("ABCC2", "bilirubin efflux MRP2"),
    ("ABCC3", "efflux MRP3"), ("ABCG2", "efflux BCRP"),
    ("G6PD", "hemolysis pathway"), ("PKLR", "pyruvate kinase"), ("SPTA1", "spectrin α"), ("HP", "haptoglobin"),
    ("IL6", "IL-6"), ("IL6R", "*IL-6 receptor"), ("IL6ST", "gp130"), ("SOCS3", "*IL-6/STAT3 target gene"),
    ("STAT3", "STAT3"), ("CRP", "*CRP"), ("SAA1", "*SAA1"), ("LBP", "LBP"),
    ("ORM1", "α1-acid glycoprotein"), ("SERPINA1", "α1-antitrypsin"), ("FGA", "fibrinogen α"), ("FGB", "fibrinogen β"),
    ("MPO", "myeloperoxidase (oxidizes bilirubin)"), ("NCF1", "NOX2 subunit"), ("CYBB", "NOX2"),
    ("SOD1", "SOD1"), ("SOD2", "SOD2"), ("CAT", "catalase"), ("GPX1", "GPX1"),
    ("NQO1", "NQO1"), ("GCLC", "glutamate-cysteine ligase"),
    ("NR0B2", "SHP (represses UGT1A1)"), ("HNF4A", "HNF4α (UGT1A1 transcription factor)"),
    # hepcidin / iron metabolism axis -- used to test the hypothesis "IL-6 -> hepcidin -> iron restriction -> erythropoiesis down -> heme turnover down -> bilirubin production down"
    ("HAMP", "*hepcidin (classically induced by IL-6, restricts iron utilization)"), ("HFE", "HFE"), ("SLC40A1", "ferroportin FPN1"),
    ("TFRC", "transferrin receptor 1"), ("FTL", "ferritin light chain"), ("FTH1", "ferritin heavy chain"),
    ("BMP6", "BMP6 (upstream of hepcidin)"), ("EPO", "erythropoietin"), ("EPOR", "EPO receptor"),
    ("HBA1", "hemoglobin α (erythroid source marker)"), ("HBB", "hemoglobin β"),
]
ACUTE_PHASE = ["SOCS3", "CRP", "SAA1", "LBP", "ORM1", "FGA", "FGB", "SERPINA1"]


def load(fp):
    with (gzip.open if fp.endswith(".gz") else open)(fp, "rt", encoding="utf-8", errors="replace") as f:
        first = f.readline().rstrip("\n")
    sep = "\t" if first.count("\t") >= first.count(",") else ","
    return pd.read_csv(fp, sep=sep, index_col=0, compression="infer"), sep


def find_gene(df, sym):
    idx = [str(x) for x in df.index]
    pat = re.compile(r"^(?:.*[|;\s])?%s(?:[|;\s].*)?$" % re.escape(sym), re.I)
    for i in idx:
        if pat.match(i.strip()):
            return i
    for i in idx:
        if sym.upper() == i.strip().upper():
            return i
    return None


def norm(df, counts=True):
    x = df.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0)
    if counts:
        tot = x.sum(axis=0).replace(0, np.nan)
        return np.log2(x.div(tot, axis=1) * 1e6 + 1)
    return np.log2(x + 1)


def groups_from_names(cols, case_pat, ctrl_pat):
    case = [c for c in cols if re.search(case_pat, str(c), re.I)]
    ctrl = [c for c in cols if re.search(ctrl_pat, str(c), re.I)]
    return case, ctrl


def analyze(logx, case, ctrl, title, donors=None):
    print()
    print("=" * 112)
    print(title)
    print("  case n=%d %s | control n=%d %s" % (len(case), case[:6], len(ctrl), ctrl[:6]))
    print("=" * 112)
    rows = []
    for sym, note in PANEL:
        gi = find_gene(logx, sym)
        if gi is None:
            rows.append(dict(gene=sym, note=note, found=False)); continue
        a = logx.loc[gi, case].astype(float).values
        b = logx.loc[gi, ctrl].astype(float).values
        rec = dict(gene=sym, note=note, found=True, n_case=len(a), n_ctrl=len(b),
                   mean_case=float(np.mean(a)), mean_ctrl=float(np.mean(b)),
                   log2FC=float(np.mean(a) - np.mean(b)))
        rec["p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) > 1 and len(b) > 1 else np.nan
        # per-donor directional concordance: pair by the r index
        if donors:
            cons = []
            for d in donors:
                ca = [c for c in case if d in str(c) and c in logx.columns]
                cb = [c for c in ctrl if d in str(c) and c in logx.columns]
                if ca and cb:
                    cons.append(np.mean(logx.loc[gi, ca].astype(float)) - np.mean(logx.loc[gi, cb].astype(float)))
            if cons:
                rec["donor_same_dir"] = "%d/%d" % (sum(1 for v in cons if (v > 0) == (rec["log2FC"] > 0)), len(cons))
        rows.append(rec)
    r = pd.DataFrame(rows)
    rr = r[r["found"] == True].dropna(subset=["log2FC"]).sort_values("p")
    print("%-9s %8s %8s %9s %10s %8s  %s" % ("gene", "mean_case", "mean_ctrl", "log2FC", "p", "donor_cons", "note"))
    print("-" * 112)
    for _, x in rr.iterrows():
        star = "***" if x.p < 0.001 else ("**" if x.p < 0.01 else ("*" if x.p < 0.05 else ""))
        print("%-9s %8.2f %8.2f %+9.2f %10.2g %8s  %s %s" % (
            x.gene, x.mean_case, x.mean_ctrl, x.log2FC, x.p, x.get("donor_same_dir", "-"),
            str(x.note)[:22], star))
    fl = rr
    if len(fl):
        up = fl[fl.log2FC > 0]; dn = fl[fl.log2FC < 0]
        print("\n  up %d / down %d (measurable %d)" % (len(up), len(dn), len(fl)))
        print("  *positive control (IL-6 acute-phase target genes): %s" % ", ".join(
            "%s%+.2f(p=%.2g)" % (g, fl[fl.gene == g].log2FC.iloc[0], fl[fl.gene == g].p.iloc[0])
            for g in ACUTE_PHASE if g in set(fl.gene)))
    r["comparison"] = title
    return r


def main():
    print("Mechanism pillar analysis   %s" % DIR)
    allres = []

    # ---------- T1: GSE202045 ----------
    MODELS = [("humanPrimaryHep", "GSE202045_humanPrimaryHep_mRNA_counts.csv.gz", "human primary hepatocytes", ["r1", "r2", "r3"]),
              ("HepG2", "GSE202045_HepG2_mRNA_counts.csv.gz", "HepG2", ["r1", "r2"]),
              ("mousePrimaryHep", "GSE202045_mousePrimaryHep_mRNA_counts.csv.gz", "mouse primary hepatocytes", ["r1", "r2", "r3"])]
    for tag, fn, label, donors in MODELS:
        fp = os.path.join(DIR, fn)
        if not os.path.exists(fp):
            print("\n[%s] file missing, skipping" % label); continue
        df, _ = load(fp)
        cols = [str(c) for c in df.columns]
        logx = norm(df, counts=True)
        print("\n\n########## T1 %s: %d genes x %d samples ##########" % (label, *df.shape))
        for hh in ["0", "6", "24"]:
            case = [c for c in cols if re.search(r"il6pos", c, re.I) and re.search(r"[.\-]%s[.\-]" % hh, c)]
            ctrl = [c for c in cols if re.search(r"il6neg", c, re.I) and re.search(r"[.\-]%s[.\-]" % hh, c)]
            if len(case) < 2 or len(ctrl) < 2:
                print("  [%sh] too few samples (%d/%d), skipping" % (hh, len(case), len(ctrl))); continue
            allres.append(analyze(logx, case, ctrl,
                                  "GSE202045 %s: IL-6 stimulated vs unstimulated, %sh" % (label, hh), donors))

    # ---------- T2: GSE25504 ----------
    fp = os.path.join(DIR, "GSE25504_fixed.txt.gz")
    if os.path.exists(fp):
        raw, sep = load(fp)
        cols = [str(c) for c in raw.columns]
        print("\n\n########## T2 GSE25504 whole blood: %d x %d ##########" % raw.shape)
        print("  first 24 sample names: %s" % cols[:24])
        case, ctrl = [], []
        for c in cols:
            lc = c.lower()
            if re.search(r"sepsis|septic|infect|case|diseased|patient", lc):
                case.append(c)
            elif re.search(r"control|healthy|normal|donor|ctl", lc):
                ctrl.append(c)
        print("  auto grouping: case=%d ctrl=%d" % (len(case), len(ctrl)))
        if len(case) >= 3 and len(ctrl) >= 3:
            allres.append(analyze(norm(raw, counts=False), case, ctrl,
                                  "GSE25504 neonatal sepsis whole blood: case vs control"))
        else:
            print("  !! grouping parsing failed, sample annotation needs manual review")

    if allres:
        d = pd.concat(allres, ignore_index=True)
        p = os.path.join(OUT, "mechanism_geo_results.csv")
        d.to_csv(p, index=False, encoding="utf-8-sig")
        print("\n\nresults saved: %s  (%d rows)" % (p, len(d)))
    else:
        print("\n!! no results")


if __name__ == "__main__":
    main()
