# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Mechanism pillar · circadian-rhythm dataset analysis (GSE218242 / GSE218243) -- honest positioning: orthogonal mechanism probing, not IL-6 replication

* Prerequisite facts (verified, no hedging allowed):
  The superseries GSE218244 (PMID 38657074, Sci Adv 2024) contains only 2 subseries, **neither with an IL-6/inflammation stimulation arm**:
    GSE218242  primary hepatocyte circadian time course (24-72h, every 3h, 2-3 replicates) -- single factor: time
    GSE218243  primary hepatocyte ARNTL(BMAL1) siRNA vs non-targeting control (36/42/48/54h, 2 replicates each)
  -> Therefore this script answers **two answerable questions** and never masquerades as IL-6 replication:
    Q1 [GSE218243] Does knockdown of the core clock gene ARNTL in hepatocytes change bilirubin metabolism / heme / iron-hepcidin axis genes?
                   (rationale: inflammation and IL-6 are known to suppress hepatic clock genes; if ARNTL down also lowers UGT1A1,
                     then IL-6 -> UGT1A1 down has a candidate mediation chain "via clock genes")
    Q2 [GSE218242] Do the bilirubin pathway genes (UGT1A1/HMOX1/HAMP/LPIN1...) themselves show a 24h rhythm?
                   (cosinor fitting; the in-panel positive control = core clock genes must have a strong amplitude)

* Gene IDs are always looked up and verified live at MyGene.info, never written from memory (once recorded CD3G as AUTS2)
"""
import os, sys, gzip, json, re, urllib.request, time
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
from scipy import stats

DIR = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"

PANEL = [
    ("HMOX1", "heme oxygenase 1 (rate-limiting for bilirubin production)"), ("HMOX2", "heme oxygenase 2"),
    ("BLVRA", "biliverdin reductase A"), ("BLVRB", "biliverdin reductase B"),
    ("ALAS1", "rate-limiting enzyme of heme synthesis"), ("FECH", "ferrochelatase"),
    ("UGT1A1", "*bilirubin conjugation (determines serum bilirubin)"), ("UGT1A6", "UGT1A6"), ("UGT1A9", "UGT1A9"),
    ("UGT2B7", "UGT2B7"), ("SLCO1B1", "hepatic uptake OATP1B1"), ("SLCO1B3", "hepatic uptake OATP1B3"),
    ("ABCC2", "bilirubin efflux MRP2"), ("ABCC3", "efflux MRP3"), ("ABCG2", "efflux BCRP"),
    ("NR0B2", "SHP (represses UGT1A1)"), ("HNF4A", "HNF4α (UGT1A1 transcription factor)"), ("NR1H4", "FXR"),
    ("HAMP", "*hepcidin (classically induced by IL-6)"), ("HFE", "HFE"), ("SLC40A1", "FPN1 ferroportin"),
    ("TFRC", "transferrin receptor 1"), ("FTL", "ferritin light chain"), ("FTH1", "ferritin heavy chain"),
    ("BMP6", "BMP6"), ("EPO", "erythropoietin"), ("EPOR", "EPO receptor"),
    ("HBA1", "hemoglobin α"), ("HBB", "hemoglobin β"),
    ("IL6", "IL-6"), ("IL6R", "*IL-6 receptor"), ("IL6ST", "gp130"), ("SOCS3", "*IL-6/STAT3 target"),
    ("STAT3", "STAT3"), ("CRP", "*CRP"), ("SAA1", "*SAA1"), ("LBP", "LBP"),
    ("ORM1", "α1-acid glycoprotein"), ("SERPINA1", "α1-antitrypsin"), ("FGA", "fibrinogen α"), ("FGB", "fibrinogen β"),
    # core clock genes (positive control)
    ("ARNTL", "*BMAL1 (the KD target of this study)"), ("CLOCK", "CLOCK"),
    ("PER1", "PER1"), ("PER2", "*PER2"), ("PER3", "PER3"),
    ("CRY1", "CRY1"), ("CRY2", "CRY2"),
    ("NR1D1", "*REV-ERBα"), ("NR1D2", "REV-ERBβ"),
    ("DBP", "*DBP"), ("TEF", "TEF"), ("HLF", "HLF"),
    ("NFIL3", "E4BP4"), ("RORA", "RORα"), ("BHLHE40", "DEC1"), ("BHLHE41", "DEC2"),
    ("PPARA", "PPARα"), ("PPARG", "PPARγ"), ("NR1H3", "LXRα"), ("CYP2E1", "CYP2E1 (clock-regulated drug-metabolizing enzyme control)"),
    # hepatocyte identity / housekeeping -- used to rule out the alternative explanation that "ARNTL KD is merely global hepatocyte dedifferentiation" (specificity check)
    ("ALB", "albumin (liver identity)"), ("APOB", "ApoB (liver identity)"), ("TTR", "TTR (liver identity)"),
    ("CYP3A4", "CYP3A4 (liver identity drug-metabolizing enzyme)"), ("CYP1A2", "CYP1A2 (liver identity drug-metabolizing enzyme)"),
    ("ACTB", "β-actin (housekeeping)"), ("GAPDH", "GAPDH (housekeeping)"), ("B2M", "B2M (housekeeping)"),
]
ACUTE = ["SOCS3", "CRP", "SAA1", "LBP", "ORM1", "FGA", "FGB", "SERPINA1"]
CLOCKGENES = ["ARNTL", "CLOCK", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2",
              "DBP", "TEF", "HLF", "NFIL3", "RORA", "BHLHE40", "BHLHE41"]


# ---------------- live lookup and verification: symbol -> ENSG ----------------
# * deprecated symbol mapping (HGNC renaming; MyGene's symbol: query can no longer resolve the old name, real tests return 0 hits)
#   ARNTL -> BMAL1 (tested 2026-10-09: `symbol:ARNTL` 0 hits, `symbol:BMAL1` -> ENSG00000133794)
DEPRECATED = {"ARNTL": "BMAL1"}


def _ensg_of(h):
    ens = h.get("ensembl")
    if isinstance(ens, list):
        ens = ens[0] if ens else {}
    return (ens or {}).get("gene") if isinstance(ens, dict) else None


def mygene_map(symbols, cache=os.path.join(OUT, "_mygene_symbol2ensg.json")):
    d = {}
    if os.path.exists(cache):
        try:
            with open(cache, encoding="utf-8") as f:
                d = json.load(f)
        except Exception:
            d = {}
    # * the cache-hit test must use "value is non-empty", not `s in d` -- otherwise None is treated as resolved and the follow-up lookup is skipped (hit this before)
    todo = [s for s in symbols if not d.get(s)]
    if todo:
        print("  MyGene lookup for %d symbols..." % len(todo))
        for i in range(0, len(todo), 40):
            chunk = todo[i:i + 40]
            q = "symbol:(" + " OR ".join(chunk) + ") AND taxid:9606"
            u = ("https://mygene.info/v3/query?q=%s&fields=symbol,ensembl.gene&species=human&size=200"
                 % urllib.parse.quote(q))
            try:
                req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
                hits = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())["hits"]
                got = {}
                for h in hits:
                    sym = (h.get("symbol") or "").upper()
                    gid = _ensg_of(h)
                    if gid and sym and sym not in got:
                        got[sym] = gid
                for s in chunk:
                    d[s] = got.get(s.upper())
            except Exception as e:
                print("    MyGene batch failed:", str(e)[:70])
            time.sleep(0.4)
    # * fallback: batch OR queries miss entries (in real tests ARNTL=BMAL1 was missed -- missing the KD target ruins the analysis)
    miss = [s for s in symbols if not d.get(s)]
    if miss:
        print("  fallback single lookups for %d: %s" % (len(miss), miss))
        for s in miss:
            for term in ["symbol:%s" % s, s, "symbol:%s" % DEPRECATED.get(s, s)]:
                if not term:
                    continue
                try:
                    u = ("https://mygene.info/v3/query?q=%s&fields=symbol,ensembl.gene,alias"
                         "&species=human&size=10" % urllib.parse.quote(term))
                    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
                    hits = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())["hits"]
                    for h in hits:
                        gid = _ensg_of(h)
                        if gid and (h.get("symbol", "").upper() in (s.upper(), DEPRECATED.get(s, "").upper())):
                            d[s] = gid
                            break
                    if d.get(s):
                        break
                except Exception as e:
                    print("    %-9s (%s) lookup failed %s" % (s, term, str(e)[:40]))
                time.sleep(0.3)
            print("    %-9s -> %s" % (s, d.get(s)))
    with open(cache, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1)
    return d


def load_l2(fp):
    """GEO l2tpm / normcounts: first column is gene ID, the remaining columns are samples."""
    df = pd.read_csv(fp, sep="\t" if fp.endswith(".txt.gz") else ",", index_col=0,
                     compression="infer", low_memory=False)
    df.index = [str(i).strip() for i in df.index]
    df = df[~df.index.duplicated(keep="first")]
    return df.apply(pd.to_numeric, errors="coerce")


# ---------------- Q1: GSE218243 ARNTL KD ----------------
def q1():
    fp = os.path.join(DIR, "GSE218243_l2tpm.txt.gz")
    if not os.path.exists(fp):
        print("!! GSE218243 missing"); return []
    df = load_l2(fp)
    cols = [str(c) for c in df.columns]
    print("\n" + "#" * 108)
    print("# Q1  GSE218243  primary hepatocytes ARNTL(BMAL1) siRNA vs non-targeting control   %d genes x %d samples" % df.shape)
    print("#" * 108)
    print("  sample columns:", cols)
    # parsing: KD group = title contains BMAL; control group = contains NT
    kd = [c for c in cols if re.search(r"BMAL", c, re.I)]
    ct = [c for c in cols if re.search(r"\bNT", c, re.I) or re.search(r"NT[ab]", c, re.I)]
    if not ct:
        ct = [c for c in cols if not re.search(r"BMAL", c, re.I)]
    print("  KD n=%d %s" % (len(kd), kd))
    print("  control n=%d %s" % (len(ct), ct))
    if len(kd) < 2 or len(ct) < 2:
        print("  !! grouping failed"); return []

    m = mygene_map([s for s, _ in PANEL])
    global _GENEMAP
    _GENEMAP = m
    rows = []
    for sym, note in PANEL:
        gid = m.get(sym)
        if not gid or gid not in df.index:
            rows.append(dict(dataset="GSE218243_ARNTL_KD", gene=sym, note=note, found=False))
            continue
        a = df.loc[gid, kd].astype(float).values
        b = df.loc[gid, ct].astype(float).values
        fc = float(np.mean(a) - np.mean(b))
        p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) > 1 and len(b) > 1 else np.nan
        # per-timepoint directional concordance: KD vs NT at 36/42/48/54h respectively
        cons, ntot = 0, 0
        for hh in ["36", "42", "48", "54"]:
            ka = [c for c in kd if hh + "h" in c]
            kb = [c for c in ct if hh + "h" in c]
            if ka and kb:
                ntot += 1
                if (np.mean(df.loc[gid, ka]) - np.mean(df.loc[gid, kb]) > 0) == (fc > 0):
                    cons += 1
        rows.append(dict(dataset="GSE218243_ARNTL_KD", gene=sym, note=note, found=True,
                         n_kd=len(a), n_ctrl=len(b), mean_kd=float(np.mean(a)), mean_ctrl=float(np.mean(b)),
                         log2FC=fc, p=p, time_same_dir=("%d/%d" % (cons, ntot)) if ntot else "-"))
    r = pd.DataFrame(rows)
    rr = r[r.found == True].sort_values("p")
    print("\n%-9s %8s %8s %9s %10s %9s  %s" % ("gene", "mean_KD", "mean_ctrl", "log2FC", "p", "time_cons", "note"))
    print("-" * 108)
    for _, x in rr.iterrows():
        star = "***" if x.p < 0.001 else ("**" if x.p < 0.01 else ("*" if x.p < 0.05 else ""))
        print("%-9s %8.2f %8.2f %+9.2f %10.2g %9s  %s %s" % (
            x.gene, x.mean_kd, x.mean_ctrl, x.log2FC, x.p, x.time_same_dir, str(x.note)[:24], star))
    # positive control assessment
    print("\n  -- positive control (did the KD actually work) --")
    for g in ["ARNTL", "CLOCK", "PER2", "NR1D1", "DBP", "BHLHE40", "BHLHE41"]:
        s = rr[rr.gene == g]
        if len(s):
            print("    %-8s log2FC=%+.2f  p=%.2g  %s" % (g, s.log2FC.iloc[0], s.p.iloc[0],
                  "✓ matches KD expectation" if (s.log2FC.iloc[0] < 0 or s.p.iloc[0] < 0.05) else "✗ does not match"))
    # direction summary
    for blk, names, title in [("bilirubin/heme", ["HMOX1", "HMOX2", "BLVRA", "BLVRB", "ALAS1", "FECH", "UGT1A1",
                                            "UGT1A6", "UGT1A9", "SLCO1B1", "SLCO1B3", "ABCC2", "ABCC3", "ABCG2"], "bilirubin pathway"),
                              ("iron/hepcidin axis", ["HAMP", "HFE", "SLC40A1", "TFRC", "FTL", "FTH1", "BMP6"], "iron-hepcidin axis"),
                              ("acute phase", ACUTE, "acute-phase genes")]:
        sub = rr[rr.gene.isin(names)]
        if len(sub):
            print("  %s: %d measured, up %d / down %d; significant (p<0.05) %d: %s" % (
                title, len(sub), int((sub.log2FC > 0).sum()), int((sub.log2FC < 0).sum()),
                int((sub.p < 0.05).sum()),
                ", ".join("%s%+.2f" % (g, f) for g, f in zip(sub[sub.p < 0.05].gene, sub[sub.p < 0.05].log2FC))))
    genomewide_shift(df, kd, ct, "GSE218243 ARNTL KD vs NT")
    return rows


def genomewide_shift(df, case, ctrl, title):
    """* key specificity check: genome-wide median log2FC ~= 0 -> panel gene changes are not a normalization/composition artifact;
       also reports the genome-wide percentile of key genes (with n=8 the p-value bottoms out, so the percentile is more convincing)."""
    exp = df[df.sum(axis=1) > 0]
    fc = exp[case].mean(axis=1) - exp[ctrl].mean(axis=1)
    print("\n  -- genome-wide shift (specificity check): %s --" % title)
    print("     measurable genes %d; median log2FC %+.3f, mean %+.3f, SD %.3f" % (
        len(fc), fc.median(), fc.mean(), fc.std()))
    print("     median ~= 0 -> no global shift, panel gene changes cannot be explained by normalization/composition artifacts")
    for g in ["ARNTL", "UGT1A1", "UGT1A6", "UGT1A9", "NR1D1", "DBP", "ORM1", "SERPINA1",
              "ALB", "APOB", "TTR", "GAPDH", "ACTB", "B2M"]:
        gid = _GENEMAP.get(g)
        if gid and gid in fc.index:
            v = fc[gid]
            pct = float((fc < v).mean() * 100)
            print("       %-9s log2FC=%+.2f   |FC| percentile ~= %.2f%% (smaller = larger change)" % (
                g, v, pct if v < 0 else 100 - pct))


_GENEMAP = {}


# ---------------- Q2: GSE218242 circadian rhythm ----------------
def cosinor(v, t):
    """Fit log2 expression ~ a + b·cos(2πt/24) + c·sin(2πt/24); return amplitude and F-test p"""
    X = np.column_stack([np.ones_like(t, dtype=float),
                         np.cos(2 * np.pi * t / 24.0), np.sin(2 * np.pi * t / 24.0)])
    y = np.asarray(v, dtype=float)
    ok = np.isfinite(y)
    if ok.sum() < 8:
        return np.nan, np.nan, np.nan
    Xo, yo = X[ok], y[ok]
    beta, *_ = np.linalg.lstsq(Xo, yo, rcond=None)
    resid = yo - Xo @ beta
    rss = float(resid @ resid)
    tss = float(((yo - yo.mean()) ** 2).sum())
    k = 3
    n = ok.sum()
    if n - k <= 0 or tss <= 0:
        return np.nan, np.nan, np.nan
    r2 = 1 - rss / tss
    F = (r2 / (k - 1)) / ((1 - r2) / (n - k)) if r2 < 1 else np.inf
    p = float(stats.f.sf(F, k - 1, n - k))
    amp = float(np.hypot(beta[1], beta[2]))
    return amp, r2, p


def q2():
    fp = os.path.join(DIR, "GSE218242_l2tpm.txt.gz")
    if not os.path.exists(fp):
        print("!! GSE218242 missing"); return []
    df = load_l2(fp)
    cols = [str(c) for c in df.columns]
    t = np.array([float(re.match(r"(\d+)", c).group(1)) for c in cols])
    step = sorted(set(np.diff(np.unique(t)).astype(int).tolist()))
    print("\n" + "#" * 108)
    print("# Q2  GSE218242  primary hepatocyte circadian time course  %d genes x %d samples (%s-%s h, step %s)" % (
        df.shape[0], df.shape[1], int(t.min()), int(t.max()), step or "?"))
    print("#" * 108)
    print("  WARNING: this is only a screening conclusion about 'presence/absence of a 24h amplitude'; samples are not fully balanced (24h 3 replicates, 27h 2 replicates), "
          "and the 72h batch number differs (d25 vs d33), so a batch effect exists and the amplitude values should not be quoted on their own.")
    m = mygene_map([s for s, _ in PANEL])
    rows = []
    for sym, note in PANEL:
        gid = m.get(sym)
        if not gid or gid not in df.index:
            rows.append(dict(dataset="GSE218242_circadian", gene=sym, note=note, found=False)); continue
        amp, r2, p = cosinor(df.loc[gid].values, t)
        rows.append(dict(dataset="GSE218242_circadian", gene=sym, note=note, found=True,
                         amplitude=amp, r2=r2, p=p, mean_log2=float(np.nanmean(df.loc[gid].values))))
    r = pd.DataFrame(rows)
    rr = r[r.found == True].sort_values("p")
    # BH correction
    from statsmodels.stats.multitest import multipletests
    try:
        rr = rr.assign(q=multipletests(rr.p.fillna(1).values, method="fdr_bh")[1])
    except Exception:
        rr = rr.assign(q=np.nan)
    print("\n%-9s %9s %8s %10s %8s  %s" % ("gene", "amplitude", "R²", "p", "q(BH)", "note"))
    print("-" * 108)
    for _, x in rr.iterrows():
        star = "***" if x.p < 0.001 else ("**" if x.p < 0.01 else ("*" if x.p < 0.05 else ""))
        print("%-9s %9.3f %8.3f %10.2g %8.3f  %s %s" % (x.gene, x.amplitude, x.r2, x.p, x.q, str(x.note)[:26], star))
    cg = rr[rr.gene.isin(CLOCKGENES)]
    print("\n  -- positive control: core clock genes (expected to have a strong 24h amplitude) --  significant (q<0.05) %d/%d" % (
        int((cg.q < 0.05).sum()), len(cg)))
    print("    " + ", ".join("%s(amp=%.2f,q=%.2g)" % (g, a, q) for g, a, q in zip(cg.gene, cg.amplitude, cg.q)))
    bili = rr[rr.gene.isin(["HMOX1", "HMOX2", "BLVRA", "BLVRB", "ALAS1", "FECH", "UGT1A1", "UGT1A6",
                            "UGT1A9", "SLCO1B1", "SLCO1B3", "ABCC2", "ABCC3", "ABCG2"])]
    print("  bilirubin pathway significant (q<0.05): %s" % (", ".join("%s(amp=%.2f)" % (g, a)
          for g, a, q in zip(bili.gene, bili.amplitude, bili.q) if q < 0.05) or "none"))
    return rows


def main():
    print("Mechanism pillar · circadian-rhythm datasets (GSE218242 / GSE218243)   %s" % DIR)
    rows = q1() + q2()
    if rows:
        d = pd.DataFrame(rows)
        p = os.path.join(OUT, "mechanism_circadian_results.csv")
        d.to_csv(p, index=False, encoding="utf-8-sig")
        print("\nresults saved: %s  (%d rows)" % (p, len(d)))


if __name__ == "__main__":
    main()
