# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""IL-6 pathway (IL6R pQTL) and SAA (SAA1/SAA2 pQTL) -> serum bilirubin  bidirectional cis-MR

Background (empirical constraints as of 2026-10-08, all reported in the Methods):
  - **IL-6 protein itself cannot serve as an instrument in OpenGWAS**: the only two IL-6 level datasets
    (ebi-a-GCST90012005 Folkersen2020 n=21,758; prot-b-2 Folkersen2018) have, within +/-50kb of the IL6 gene,
    their strongest association at only p~=3e-3 / 1.8e-3, **failing to reach genome-wide significance**.
  - This script therefore uses two workable exposures instead:
      (1) **IL6R (IL-6 signalling pathway)**: a very strong cis-pQTL (prot-b-23 p=2.2e-265, F~=1457);
         this is the standard approach for MR of the IL-6 pathway (IL-6R blockade = the tocilizumab target)
      (2) **SAA1/SAA2 (acute-phase proteins)**: ebi-a-GCST90019394 / 90019411 (Pietzner 2020 Olink)
  - Positive control: IL6R -> CRP. It is well established that IL6R variants lower CRP; failure to detect this in our pipeline would indicate an error in the pipeline.

WARNING: the token is read from the environment variable OGWAS_JWT and is never written to this file.
Method notes: cis region (gene +/-50kb) + distance pruning (500kb approximating LD clumping, to be declared in the Methods);
all coordinates are resolved dynamically via MyGene (this project has already hit the "writing Ensembl IDs from memory" pitfall twice).
"""
import os, sys, json, math, time, csv, urllib.request, urllib.error, urllib.parse
sys.stdout.reconfigure(encoding="utf-8")

JWT = os.environ.get("OGWAS_JWT", "")
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"
LOG = os.path.join(OUT, "_mr_il6_saa.log")
P_CIS, PAD, PRUNE = 5e-8, 50_000, 500_000
_L = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    _L.append(s)
    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(_L))


def og(path, body=None, timeout=180, tries=4):
    if not JWT:
        raise SystemExit("!! OGWAS_JWT not set")
    for i in range(tries):
        try:
            data = json.dumps(body).encode() if body is not None else None
            req = urllib.request.Request("https://api.opengwas.io" + path, data=data, headers={
                "Authorization": "Bearer " + JWT, "Content-Type": "application/json",
                "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i == tries - 1:
                P("   !! API %s" % str(e)[:100]); return None
            time.sleep(3)


def web(u, timeout=30):
    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def gene_coord(symbol):
    d = web("https://mygene.info/v3/query?" + urllib.parse.urlencode(
        {"q": "symbol:%s AND taxid:9606" % symbol,
         "fields": "ensembl.gene,symbol,name,genomic_pos_hg19", "size": 5}))
    for h in d.get("hits", []):
        if h.get("symbol") != symbol:
            continue
        # * MyGene's ensembl may be a dict ({"gene": "ENSG..."}) or a list (multiple transcripts/multiple entries)
        eg = h.get("ensembl")
        if isinstance(eg, list):
            eg = eg[0] if eg else None
        if isinstance(eg, dict):
            eg = eg.get("gene")
        if isinstance(eg, list):
            eg = eg[0] if eg else None
        g19 = h.get("genomic_pos_hg19")
        if isinstance(g19, list):
            g19 = g19[0] if g19 else None
        if eg and isinstance(g19, dict) and g19.get("start"):
            return eg, g19.get("chr"), int(g19["start"]), int(g19["end"]), h.get("name")
    return None, None, None, None, None


def get_cis(gid, ch, s, e, label, quiet=False):
    var = "%s:%d-%d" % (ch, s - PAD, e + PAD)
    r = og("/api/associations", {"id": [gid], "variant": [var]})
    rows = r if isinstance(r, list) else []
    sig = [x for x in rows if (x.get("p") or 1) < P_CIS]
    sig.sort(key=lambda x: x["p"])
    kept = []
    for x in sig:
        if all(abs(x["position"] - k["position"]) > PRUNE for k in kept):
            kept.append(x)
    if not quiet:
        best = ("%.3e" % rows[0]["p"]) if False else (min((x.get("p", 1) for x in rows), default=float("nan")))
        P("   %-30s %-24s returned %-5d p<5e-8: %-5d after pruning: %d" % (label, var, len(rows), len(sig), len(kept)))
        for k in kept[:4]:
            P("        %-13s chr%s:%-10s p=%.3e beta=%+.4f se=%.4f ea=%s nea=%s eaf=%s F=%.0f"
              % (k.get("rsid"), k.get("chr"), k.get("position"), k["p"], k["beta"], k["se"],
                 k.get("ea"), k.get("nea"), k.get("eaf"), (k["beta"] / k["se"]) ** 2))
        if not kept:
            P("        WARNING: no genome-wide significant SNP (this exposure cannot support cis-MR)")
    return kept


def fetch(gid, snps):
    if not snps:
        return {}
    r = og("/api/associations", {"id": [gid], "variant": list(snps)})
    return {x["rsid"]: x for x in (r if isinstance(r, list) else []) if x.get("rsid")}


def align(a_ea, a_nea, b_ea, b_nea, bx, bxse, verbose=True, tag=""):
    if a_ea is None or b_ea is None:
        if verbose:
            P("        %-13s WARNING: allele information missing, skipping" % tag)
        return None, None, "noallele"
    if a_ea == b_ea and a_nea == b_nea:
        return bx, bxse, "same"
    if a_ea == b_nea and a_nea == b_ea:
        return -bx, bxse, "flip"
    if verbose:
        P("        %-13s WARNING: allele mismatch exposure %s/%s vs outcome %s/%s" % (tag, a_ea, a_nea, b_ea, b_nea))
    return None, None, "mismatch"


def run_mr(expo, outmap, label, verbose=True):
    bx, bxse, by, byse, used = [], [], [], [], []
    for e in expo:
        rs = e.get("rsid")
        o = outmap.get(rs)
        if not o:
            if verbose:
                P("        %-13s WARNING: missing in outcome" % rs)
            continue
        nx, nsx, st = align(e["ea"], e["nea"], o.get("ea"), o.get("nea"), e["beta"], e["se"], verbose, rs)
        if nx is None:
            continue
        bx.append(nx); bxse.append(nsx); by.append(o["beta"]); byse.append(o["se"]); used.append(rs)
        if verbose:
            P("        SNP %-13s exposure b=%+.4f/se=%.4f(%s/%s) | outcome b=%+.4f/se=%.4f(%s/%s) | %s"
              % (rs, e["beta"], e["se"], e["ea"], e["nea"], o["beta"], o["se"],
                 o.get("ea"), o.get("nea"), st))
    if not bx:
        P("   -> %-46s no usable SNP" % label)
        return None
    bx = [float(v) for v in bx]; bxse = [float(v) for v in bxse]
    by = [float(v) for v in by]; byse = [float(v) for v in byse]
    if len(bx) == 1:
        b = by[0] / bx[0]; se = abs(byse[0] / bx[0]); meth = "Wald (1 SNP)"
    else:
        w = [1.0 / ((s_ / b_) ** 2) for b_, s_ in zip(bx, bxse)]
        b = sum(wi * (byi / bxi) for wi, bxi, byi in zip(w, bx, by)) / sum(w)
        se = math.sqrt(1.0 / sum(w)); meth = "IVW (%d SNP)" % len(bx)
    z = b / se if se else float("nan")
    p = math.erfc(abs(z) / math.sqrt(2)) if se else float("nan")
    P("   -> %-46s %-14s b=%+.4f SE=%.4f p=%.4g%s"
      % (label, meth, b, se, p, "  **" if p < 0.05 else ""))
    return dict(label=label, method=meth, nsnp=len(bx), b=b, se=se, p=p, snps=",".join(used))


def main():
    t0 = time.time()
    P("=" * 100)
    P("IL-6 pathway (IL6R) and SAA -> serum bilirubin  bidirectional cis-MR      %s" % time.strftime("%Y-%m-%d %H:%M"))
    P("=" * 100)
    P("")
    P("[0] Dynamically resolve gene coordinates (not from memory)")
    co = {}
    for sym in ("IL6", "IL6R", "SAA1", "SAA2", "UGT1A1", "CRP"):
        ens, ch, s, e, nm = gene_coord(sym)
        co[sym] = (ens, ch, s, e)
        P("   %-7s -> %-18s %-38s hg19:chr%s:%s-%s" % (sym, ens, str(nm)[:38], ch, s, e))
    P("")

    P("=" * 100)
    P("[1] Exposure and outcome dataset metadata (empirically verified, written into the Methods)")
    P("=" * 100)
    IDS = ["ebi-a-GCST90012025", "ebi-a-GCST90019398", "prot-a-1540", "prot-b-23",
           "ebi-a-GCST90019394", "ebi-a-GCST90019411", "ebi-a-GCST90012005",
           "ebi-a-GCST90029070", "ieu-b-35",
           "ebi-a-GCST90025973", "ebi-a-GCST90025983", "bbj-a-53", "ukb-d-30840_irnt",
           "finn-b-P16_NEONTAL_JAUND_OTH_UNSP_CAUSES", "finn-b-P16_BACTERIAL_SEPSIS_NEWBO",
           "finn-b-P16_DISORD_RELATED_GESTATION_LOW_BIRTHWGHTT_NECIFIED",
           "finn-b-P16_RESPI_DISTR_NEWBO"]
    meta = {}
    r = og("/api/gwasinfo", {"id": IDS})
    for v in (r if isinstance(r, list) else []):
        meta[v.get("id")] = v
        P("   %-46s %-40s N=%-8s ncase=%-6s %s %s"
          % (v.get("id"), str(v.get("trait"))[:40], v.get("sample_size"), v.get("ncase"),
             v.get("author"), v.get("year")))
    for c in IDS:
        if c not in meta:
            P("   %-46s [FAIL] does not exist (this outcome is unavailable)" % c)
    P("")

    P("=" * 100)
    P("[2] cis instrument extraction")
    P("=" * 100)
    ex = {}
    P(" > IL6R (IL-6 signalling pathway)")
    for gid in ["ebi-a-GCST90012025", "ebi-a-GCST90019398", "prot-a-1540", "prot-b-23"]:
        if gid in meta:
            ex["IL6R:" + gid] = get_cis(gid, co["IL6R"][1], co["IL6R"][2], co["IL6R"][3], "IL6R %s" % gid)
    P("")
    P(" > SAA1 / SAA2 (acute-phase proteins)")
    for gid in ["ebi-a-GCST90019394", "ebi-a-GCST90019411"]:
        if gid in meta:
            # SAA1 and SAA2 are in the same cluster (hg19 chr11:18.25-18.29Mb); the window takes the union of both genes
            s1 = min(co["SAA1"][2] or 9e9, co["SAA2"][2] or 9e9)
            e1 = max(co["SAA1"][3] or 0, co["SAA2"][3] or 0)
            ch = co["SAA1"][1] or co["SAA2"][1]
            ex["SAA:" + gid] = get_cis(gid, ch, s1, e1, "SAA %s" % gid)
    P("")
    P(" > for reverse MR: bilirubin instrument (UGT1A1) + for the positive control: CRP instrument")
    bil_ins = get_cis("ebi-a-GCST90025973", co["UGT1A1"][1], co["UGT1A1"][2], co["UGT1A1"][3], "Bil(Barton21)")
    crp_ins = get_cis("ebi-a-GCST90029070", co["CRP"][1], co["CRP"][2], co["CRP"][3], "CRP(Said22)")
    P("")

    results = []
    P("=" * 100)
    P("[3] * Positive control: IL6R -> CRP (IL6R variants are known to lower CRP; failure to detect this indicates a pipeline error)")
    P("=" * 100)
    for gid in ["ebi-a-GCST90029070", "ieu-b-35"]:
        om = fetch(gid, [x["rsid"] for x in ex.get("IL6R:ebi-a-GCST90012025", [])])
        rr = run_mr(ex.get("IL6R:ebi-a-GCST90012025", []), om, "IL6R -> CRP [%s]" % gid)
        if rr:
            results.append(rr)
    om = fetch("ebi-a-GCST90029070", [x["rsid"] for x in ex.get("IL6R:prot-b-23", [])])
    rr = run_mr(ex.get("IL6R:prot-b-23", []), om, "IL6R(prot-b-23) -> CRP [Said22]")
    if rr:
        results.append(rr)
    P("")

    P("=" * 100)
    P("[4] Main analysis: IL6R / SAA -> bilirubin")
    P("=" * 100)
    OUTS_BIL = [("ebi-a-GCST90025973", "total bilirubin Barton21"), ("ukb-d-30840_irnt", "total bilirubin UKB-Neale"),
                ("bbj-a-53", "bilirubin BBJ (East Asian)"), ("ebi-a-GCST90025983", "direct bilirubin")]
    for key in ["IL6R:ebi-a-GCST90012025", "IL6R:prot-b-23", "IL6R:prot-a-1540", "IL6R:ebi-a-GCST90019398",
                "SAA:ebi-a-GCST90019394", "SAA:ebi-a-GCST90019411"]:
        if key not in ex or not ex[key]:
            continue
        P("")
        P(" > %s" % key)
        for gid, lab in OUTS_BIL:
            if gid not in meta:
                continue
            om = fetch(gid, [x["rsid"] for x in ex[key]])
            rr = run_mr(ex[key], om, "%s -> %s" % (key.split(":")[0], lab))
            if rr:
                results.append(rr)
    P("")

    P("=" * 100)
    P("[5] Secondary analysis: IL6R / SAA -> neonatal outcomes (FinnGen)")
    P("=" * 100)
    for key in ["IL6R:ebi-a-GCST90012025", "SAA:ebi-a-GCST90019394"]:
        if key not in ex or not ex[key]:
            continue
        P("")
        P(" > %s" % key)
        for gid in ["finn-b-P16_NEONTAL_JAUND_OTH_UNSP_CAUSES", "finn-b-P16_BACTERIAL_SEPSIS_NEWBO",
                    "finn-b-P16_DISORD_RELATED_GESTATION_LOW_BIRTHWGHTT_NECIFIED",
                    "finn-b-P16_RESPI_DISTR_NEWBO"]:
            if gid not in meta:
                continue
            om = fetch(gid, [x["rsid"] for x in ex[key]])
            rr = run_mr(ex[key], om, "%s -> %s" % (key.split(":")[0], gid.replace("finn-b-", "")))
            if rr:
                results.append(rr)
    P("")

    P("=" * 100)
    P("[6] Reverse MR: bilirubin -> IL6R / SAA (using the UGT1A1 cis instrument)")
    P("=" * 100)
    for gid in ["ebi-a-GCST90012025", "ebi-a-GCST90019394", "ebi-a-GCST90029070"]:
        if gid not in meta:
            continue
        om = fetch(gid, [x["rsid"] for x in bil_ins])
        rr = run_mr(bil_ins, om, "Bilirubin -> %s" % gid)
        if rr:
            results.append(rr)
    P("")

    # Reverse: CRP -> IL6R (should show a strong negative correlation, further validating the pipeline)
    P("=" * 100)
    P("[7] Additional check: CRP -> IL6R (expected strongly negative; IL6R variants lower both CRP and IL-6 levels)")
    P("=" * 100)
    for gid in ["ebi-a-GCST90012025", "prot-b-23"]:
        if gid not in meta:
            continue
        om = fetch(gid, [x["rsid"] for x in crp_ins])
        rr = run_mr(crp_ins, om, "CRP -> %s" % gid)
        if rr:
            results.append(rr)
    P("")

    os.makedirs(OUT, exist_ok=True)
    p1 = os.path.join(OUT, "MR_il6_saa_bilirubin_result.csv")
    with open(p1, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["label", "method", "nsnp", "b", "se", "p", "snps"])
        w.writeheader(); w.writerows(results)
    with open(os.path.join(OUT, "MR_il6_saa_instruments.json"), "w", encoding="utf-8") as f:
        json.dump({"instruments": ex, "bilirubin": bil_ins, "crp": crp_ins,
                   "coord": {k: list(v) for k, v in co.items()}}, f, ensure_ascii=False, indent=1)
    P("")
    P("=" * 100)
    P("%d MR estimates in total; results saved to %s" % (len(results), p1))
    P("Interpretation: b is the change in the outcome per **1-SD increase in the exposure** (see each dataset for outcome units; bilirubin is mostly in mg/dL or SD)")
    P("Total elapsed %.1f min" % ((time.time() - t0) / 60))
    P("=" * 100)


if __name__ == "__main__":
    main()
