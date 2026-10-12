# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Mechanism pillar · **independent replication**: inflammatory stimulation -> UGT1A1 / bilirubin-iron axis in human hepatocytes

Aim: GSE202045 (human primary hepatocytes + IL-6) already gave UGT1A1 down (24h, 3/3 donors concordant).
     This script replicates with **additional, methodologically independent human hepatocyte systems**, avoiding "isolated single-dataset evidence".

Datasets (grouping is always parsed from GEO metadata, never guessed):
  R1  GSE230325  dHepaRG (differentiated human liver progenitor cell line, with primary hepatocyte features) +/- LPS, 0/6/24h, 3 replicates each
                 -- the original abstract states: LPS induces acute-phase genes via **autocrine IL-6 signaling** while downregulating liver homeostatic function
  R2  GSE275342  Hep3B + Huh7 (human liver-cancer-derived hepatocyte lines) +/- IL-6, RNA-seq
  R3  GSE102006  HepG2 + cytokine mix (TNF+IFNγ+IL-1α+IL-6+LPS), 6/12/24h (array)

* The positive control is a hard gate: if the acute-phase genes (CRP/SAA1/LBP/ORM1/SERPINA1/SOCS3) are not induced,
  the stimulation model / grouping parsing is wrong, and none of the remaining results from that dataset may be trusted (a verdict is printed directly).
"""
import os, sys, gzip, re, json, tarfile, glob
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"

PANEL = [
    ("UGT1A1", "*bilirubin conjugation (determines serum bilirubin)"), ("UGT1A6", "UGT1A6"), ("UGT1A9", "UGT1A9"),
    ("UGT2B7", "UGT2B7"),
    ("HMOX1", "heme oxygenase 1 (bilirubin production)"), ("HMOX2", "heme oxygenase 2"),
    ("BLVRA", "biliverdin reductase A"), ("BLVRB", "biliverdin reductase B"),
    ("ALAS1", "rate-limiting enzyme of heme synthesis"), ("FECH", "ferrochelatase"),
    ("SLCO1B1", "hepatic uptake OATP1B1"), ("SLCO1B3", "hepatic uptake OATP1B3"),
    ("ABCC2", "bilirubin efflux MRP2"), ("ABCC3", "efflux MRP3"), ("ABCG2", "efflux BCRP"),
    ("NR0B2", "SHP (represses UGT1A1)"), ("HNF4A", "HNF4α (UGT1A1 transcription factor)"), ("NR1H4", "FXR"),
    ("HAMP", "*hepcidin (classically induced by IL-6)"), ("HFE", "HFE"), ("SLC40A1", "FPN1"),
    ("TFRC", "transferrin receptor 1"), ("FTL", "ferritin light chain"), ("FTH1", "ferritin heavy chain"),
    ("BMP6", "BMP6"), ("EPO", "erythropoietin"), ("EPOR", "EPO receptor"),
    ("HBA1", "hemoglobin α"), ("HBB", "hemoglobin β"),
    ("IL6", "*IL-6"), ("IL6R", "*IL-6 receptor"), ("IL6ST", "gp130"),
    ("SOCS3", "*IL-6/STAT3 target"), ("STAT3", "STAT3"), ("JAK1", "JAK1"),
    ("CRP", "*CRP"), ("SAA1", "*SAA1"), ("SAA2", "SAA2"), ("LBP", "*LBP"),
    ("ORM1", "*α1-acid glycoprotein"), ("ORM2", "ORM2"), ("SERPINA1", "*α1-antitrypsin"),
    ("FGA", "fibrinogen α"), ("FGB", "fibrinogen β"), ("FGG", "fibrinogen γ"),
    ("HP", "haptoglobin"), ("C3", "complement C3"), ("APOA1", "ApoA1"), ("TTR", "transthyretin"),
    ("MPO", "myeloperoxidase"), ("SOD1", "SOD1"), ("SOD2", "SOD2"), ("CAT", "catalase"),
    ("GPX1", "GPX1"), ("NQO1", "NQO1"), ("GCLC", "GCLC"),
    ("CYP1A2", "CYP1A2 (liver homeostatic drug-metabolizing enzyme control; inflammation should downregulate it)"), ("CYP3A4", "CYP3A4 (as above)"),
    ("CYP2C9", "CYP2C9"), ("ALB", "albumin (liver homeostasis)"), ("APOB", "ApoB"),
    ("ARNTL", "BMAL1 (clock gene; inflammation should downregulate it)"), ("NR1D1", "REV-ERBα"), ("PER2", "PER2"), ("DBP", "DBP"),
    ("NFKB1", "NF-κB1"), ("RELA", "RELA(p65)"), ("IL1B", "IL-1β"), ("TNFAIP3", "A20"),
]
ACUTE = ["CRP", "SAA1", "SAA2", "LBP", "ORM1", "ORM2", "SERPINA1", "FGA", "FGB", "FGG", "HP", "C3"]
STAT3T = ["SOCS3", "IL6", "STAT3"]

DEPRECATED = {"ARNTL": "BMAL1"}


def _ensg_of(h):
    ens = h.get("ensembl")
    if isinstance(ens, list):
        ens = ens[0] if ens else {}
    return (ens or {}).get("gene") if isinstance(ens, dict) else None


def mygene_map(symbols, cache=os.path.join(OUT, "_mygene_symbol2ensg2.json")):
    import urllib.request, urllib.parse, time
    d = {}
    if os.path.exists(cache):
        try:
            d = json.load(open(cache, encoding="utf-8"))
        except Exception:
            d = {}
    todo = [s for s in symbols if not d.get(s)]
    if todo:
        print("  MyGene lookup for %d symbols..." % len(todo))
        for i in range(0, len(todo), 40):
            chunk = todo[i:i + 40]
            q = "symbol:(" + " OR ".join(chunk) + ") AND taxid:9606"
            u = "https://mygene.info/v3/query?q=%s&fields=symbol,ensembl.gene&species=human&size=200" % urllib.parse.quote(q)
            try:
                req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
                hits = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())["hits"]
                got = {}
                for h in hits:
                    sym = (h.get("symbol") or "").upper()
                    g = _ensg_of(h)
                    if g and sym and sym not in got:
                        got[sym] = g
                for s in chunk:
                    d[s] = got.get(s.upper())
            except Exception as e:
                print("    batch failed:", str(e)[:70])
            time.sleep(0.4)
    miss = [s for s in symbols if not d.get(s)]
    for s in miss:
        for term in ["symbol:%s" % s, s, "symbol:%s" % DEPRECATED.get(s, s)]:
            if not term:
                continue
            try:
                u = ("https://mygene.info/v3/query?q=%s&fields=symbol,ensembl.gene&species=human&size=10"
                     % urllib.parse.quote(term))
                req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
                hits = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())["hits"]
                for h in hits:
                    g = _ensg_of(h)
                    if g and (h.get("symbol", "").upper() in (s.upper(), DEPRECATED.get(s, "").upper())):
                        d[s] = g
                        break
                if d.get(s):
                    break
            except Exception:
                pass
            time.sleep(0.3)
        if not d.get(s):
            print("    WARNING: unresolved:", s)
    json.dump(d, open(cache, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return d


def analyze(M, contrasts, ctrl, title, mapd=None):
    """M: log2 expression matrix; contrasts: [(name, list of column names)]; ctrl: list of control column names
    * with n=3 vs 3, the two-sided Mann-Whitney minimum p = 2/C(6,3) = 0.1, so the p-value bottoms out by nature;
      therefore the **primary criterion is changed to log2FC + pairwise directional concordance**, and the MWU p is only a descriptive reference.
    """
    print("\n" + "=" * 112)
    print(title)
    print("  control n=%d %s" % (len(ctrl), ctrl))
    print("=" * 112)
    idx = set(str(x) for x in M.index)

    def inv(sym):
        for cand in [sym, mapd.get(sym) if mapd else None]:
            if cand and cand in idx:
                return cand
        return None

    rows = []
    for gn, cols in contrasts:
        print("\n  - %s  (n=%d %s)" % (gn, len(cols), cols[:6]))
        print("  %-9s %8s %8s %9s %7s %10s  %s" % ("gene", "stim", "ctrl", "log2FC", "dir_cons", "p_MWU", "note"))
        print("  " + "-" * 106)
        sub = []
        for sym, nt in PANEL:
            k = inv(sym)
            if k is None:
                rows.append(dict(dataset=title, contrast=gn, gene=sym, note=nt, found=False)); continue
            a = M.loc[k, cols].astype(float).values
            b = M.loc[k, ctrl].astype(float).values
            fc = float(np.mean(a) - np.mean(b))
            # pairwise directional concordance (3x3=9 pairs)
            pairs = [(x - y) for x in a for y in b]
            cons = sum(1 for v in pairs if (v > 0) == (fc > 0)) if fc != 0 else 0
            p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) > 1 and len(b) > 1 else np.nan
            rec = dict(dataset=title, contrast=gn, gene=sym, note=nt, found=True,
                       mean_case=float(np.mean(a)), mean_ctrl=float(np.mean(b)),
                       log2FC=fc, p=p, pair_consist=("%d/%d" % (cons, len(pairs))))
            rows.append(rec); sub.append(rec)
        sub = sorted(sub, key=lambda r: -abs(r["log2FC"]))
        for x in sub:
            mark = "**" if x["pair_consist"] == "9/9" and abs(x["log2FC"]) >= 1 else ("*" if x["pair_consist"] == "9/9" else "")
            print("  %-9s %8.2f %8.2f %+9.2f %7s %10.2g  %s %s" % (
                x["gene"], x["mean_case"], x["mean_ctrl"], x["log2FC"], x["pair_consist"],
                x["p"], str(x["note"])[:26], mark))
        # hard validity gate
        npair = "%d/%d" % (len(cols) * len(ctrl), len(cols) * len(ctrl))
        ac = [x for x in sub if x["gene"] in ACUTE]
        up = [x for x in ac if x["log2FC"] > 0.5 and x["pair_consist"] == npair]
        print("    -- hard validity gate (acute-phase genes require %s concordant upregulation and log2FC>0.5) -- %d/%d passed: %s" % (
            npair, len(up), len(ac), ", ".join("%s%+.2f" % (x["gene"], x["log2FC"]) for x in up)))
        ok = len(up) >= 3
        print("       verdict: %s" % ("✓ stimulation model active, results are trustworthy" if ok else "✗ model not active -- results from this dataset must not be trusted"))
        key = {x["gene"]: x for x in sub}
        for g in ["UGT1A1", "UGT1A6", "UGT1A9", "HMOX1", "HAMP", "CYP3A4", "CYP1A2", "ALB", "ARNTL", "NR1D1", "IL6", "SOCS3"]:
            if g in key:
                print("        %-8s log2FC=%+.2f  %s" % (g, key[g]["log2FC"], key[g]["pair_consist"]))
        for r in rows:
            if r.get("contrast") == gn:
                r["model_ok"] = ok
    return rows


# ---------------- R1 GSE230325 ----------------
def r1():
    fp = os.path.join(DIR, "GSE230325_normcounts.csv.gz")
    if not os.path.exists(fp):
        print("!! GSE230325 missing"); return []
    df = pd.read_csv(fp, index_col=0, compression="infer")
    df.index = [str(i).strip() for i in df.index]
    df = df[~df.index.duplicated(keep="first")].apply(pd.to_numeric, errors="coerce").fillna(0)
    cols = [str(c) for c in df.columns]
    lg = np.log2(df + 1)
    mp = mygene_map([s for s, _ in PANEL])
    c = [x for x in cols if re.search(r"control", x, re.I)]
    l6 = [x for x in cols if re.search(r"6h", x, re.I)]
    l24 = [x for x in cols if re.search(r"24h", x, re.I)]
    print("\n\n########## R1  GSE230325  dHepaRG +/- LPS (human hepatocyte line)  %d x %d ##########" % df.shape)
    print("  WARNING: n=3 vs 3: two-sided Mann-Whitney minimum p = 0.1 (bottoms out), so the primary criterion is log2FC + 9/9 pairwise directional concordance")
    return analyze(lg, [("LPS 6h vs 0h", l6), ("LPS 24h vs 0h", l24)], c,
                   "GSE230325 dHepaRG LPS 0/6/24h", mapd=mp)


# ---------------- R2 GSE275342 ----------------
def r2():
    tar = os.path.join(DIR, "GSE275342_RAW.tar")
    if not os.path.exists(tar) or os.path.getsize(tar) == 0:
        print("\n\n########## R2  GSE275342  skipped: RAW.tar has not finished downloading ##########")
        return []
    exdir = os.path.join(DIR, "GSE275342_raw")
    os.makedirs(exdir, exist_ok=True)
    if not glob.glob(os.path.join(exdir, "*")):
        with tarfile.open(tar) as t:
            t.extractall(exdir)
    files = sorted(glob.glob(os.path.join(exdir, "*")))
    print("\n\n########## R2  GSE275342  Hep3B/Huh7 +/- IL-6 (100 ng/mL, 24h) ##########")
    print("  %d files inside RAW" % len(files))
    if not files:
        print("  !! no files inside tar"); return []
    # * the filename carries the grouping: GSM..._Hep3BnIL6_rep1.txt.gz -- no guessing, parse the filename directly
    mat = {}
    for fp in files:
        base = os.path.basename(fp)
        # * GSE275342's suppl is StringTie/RSEM style:
        #     index = `SYMBOL_SYMBOL` (e.g. A1BG_A1BG), columns = gene_id, transcript_id(s), length,
        #     effective_length, expected_count, TPM, FPKM
        #   -> do not take the "last numeric column" (that is FPKM); **prefer TPM**; the index must drop the duplicated `_SYMBOL` suffix
        op = gzip.open if base.endswith(".gz") else open
        try:
            with op(fp, "rt", encoding="utf-8", errors="replace") as fh:
                d = pd.read_csv(fh, sep="\t", comment="#", index_col=0)
        except Exception as e:
            print("    read failed %s : %s" % (base, str(e)[:60])); continue
        col = None
        for cand in ["TPM", "tpm", "expected_count", "counts"]:
            if cand in d.columns:
                col = cand; break
        if col is None:
            num = d.select_dtypes(include=[np.number])
            col = num.columns[-1] if len(num.columns) else d.columns[-1]
        s = pd.to_numeric(d[col], errors="coerce").fillna(0)
        s.index = [str(x).strip().split("_")[0] for x in s.index]   # A1BG_A1BG -> A1BG
        s = s[~s.index.duplicated(keep="first")]
        mat[base] = s
    if not mat:
        print("  !! no expression values recognized"); return []
    X = pd.DataFrame(mat).fillna(0)
    lg = np.log2(X + 1)
    print("  assembled matrix %d genes x %d samples; examples:" % lg.shape)
    for c in list(lg.columns)[:4]:
        print("    ", c)
    mp = mygene_map([s for s, _ in PANEL])
    rows = []
    for cl in ["Hep3B", "Huh7"]:
        case = [c for c in lg.columns if cl in c and re.search(r"pIL6", c)]
        ctrl = [c for c in lg.columns if cl in c and re.search(r"nIL6", c)]
        if len(case) >= 2 and len(ctrl) >= 2:
            rows += analyze(lg, [("IL-6 vs vehicle (%s)" % cl, case)], ctrl,
                            "GSE275342 %s human liver cancer cell line IL-6 24h" % cl, mapd=mp)
    # * no "two-line merge": the expression baselines of Hep3B and Huh7 differ greatly (UGT1A1 baseline differs by ~0.6 log2,
    #   CYP3A4 by ~4 orders of magnitude), and merging would only introduce between-line variance without adding power (removed in the previous version).
    print("\n  WARNING: no Hep3B+Huh7 merged analysis: the two lines differ too much in metabolic phenotype baseline, so merging has no statistical meaning")
    return rows


def main():
    print("Mechanism pillar · independent replication (inflammatory stimulation of human hepatocytes)   %s" % DIR)
    rows = r1() + r2()
    if rows:
        d = pd.DataFrame(rows)
        p = os.path.join(OUT, "mechanism_replication_results.csv")
        d.to_csv(p, index=False, encoding="utf-8-sig")
        print("\nresults saved: %s  (%d rows)" % (p, len(d)))


if __name__ == "__main__":
    main()
