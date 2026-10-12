# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Main analysis of the new direction: bidirectional cis-MR between CRP (inflammation) <-> serum bilirubin

WARNING: the token is read from the environment variable OGWAS_JWT and is not written to this file.
   Usage: export OGWAS_JWT='...'  ;  python this_script

Method notes (all choices made under empirical constraints)
  - The instrument uses the **cis region** (gene +/-50kb), because SERVICES__PIPELINE=Error makes /api/tophits unavailable,
    so genome-wide clumping is not possible; cis-MR is itself the standard approach for pQTLs such as CRP/bilirubin
  - Distance pruning (500 kb) approximates LD clumping -- this approximation must be declared in the Methods
  - Gene coordinates are **resolved dynamically** (MyGene) and never hard-coded: this project has already hit the "writing Ensembl IDs from memory" pitfall twice
  - Sample overlap: the two datasets in the main analysis may both derive from UKB; the **sensitivity analysis switches to Ligthart2018 (CRP) x BioBankJapan (bilirubin)**,
    different studies and different populations (East Asian), with no overlap
"""
import os, sys, json, math, time, urllib.request, urllib.error, urllib.parse
sys.stdout.reconfigure(encoding="utf-8")

JWT = os.environ.get("OGWAS_JWT", "")
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"
P_CIS = 5e-8
PAD = 50_000          # cis window extension
PRUNE = 500_000       # distance-pruning window


def og(path, body=None, timeout=120, tries=4):
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
                print("   !! %s" % str(e)[:100]); return None
            time.sleep(2)


def web(u, timeout=30):
    req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def gene_coord(symbol):
    """Dynamically resolve the gene's Ensembl ID and hg19 coordinates"""
    d = web("https://mygene.info/v3/query?" + urllib.parse.urlencode(
        {"q": "symbol:%s AND taxid:9606" % symbol,
         "fields": "ensembl.gene,symbol,name,genomic_pos_hg19", "size": 5}))
    for h in d.get("hits", []):
        if h.get("symbol") != symbol:
            continue
        eg = h.get("ensembl", {}).get("gene")
        if isinstance(eg, list):
            eg = eg[0]
        g19 = h.get("genomic_pos_hg19")
        if isinstance(g19, list):
            g19 = g19[0]
        if eg and g19:
            return eg, g19.get("chr"), int(g19["start"]), int(g19["end"]), h.get("name")
    return None, None, None, None, None


def get_cis(gwas_id, ch, s, e, label):
    """Fetch cis-region SNPs, filter genome-wide significant, distance-prune"""
    var = "%s:%d-%d" % (ch, s - PAD, e + PAD)
    r = og("/api/associations", {"id": [gwas_id], "variant": [var]})
    rows = r if isinstance(r, list) else []
    sig = [x for x in rows if (x.get("p") or 1) < P_CIS]
    sig.sort(key=lambda x: x["p"])
    kept = []
    for x in sig:
        if all(abs(x["position"] - k["position"]) > PRUNE for k in kept):
            kept.append(x)
    print("   %-12s %s  %s  returned %-4d  p<5e-8: %-4d  after pruning: %d"
          % (label, gwas_id, var, len(rows), len(sig), len(kept)))
    return kept


def fetch_outcome(gwas_id, snps):
    if not snps:
        return {}
    r = og("/api/associations", {"id": [gwas_id], "variant": snps})
    return {x["rsid"]: x for x in (r if isinstance(r, list) else [])}


def align(a_ea, a_nea, b_ea, b_nea, bx, bxse):
    """Align the exposure effect to the outcome allele direction; returns (bx, bxse, status)"""
    if a_ea == b_ea and a_nea == b_nea:
        return bx, bxse, "same"
    if a_ea == b_nea and a_nea == b_ea:
        return -bx, bxse, "flip"
    return None, None, "mismatch"


def run_mr(expo, outmap, label):
    """IVW (multiple SNPs) / Wald ratio (single SNP)"""
    bx, bxse, by, byse, used, miss = [], [], [], [], [], []
    for e in expo:
        rs = e.get("rsid")
        o = outmap.get(rs)
        if not o:
            miss.append(rs)
            continue
        # * FIX: the beta/se passed in must be the **exposure's** (previously the outcome's were passed by mistake, making b constantly equal to 1.0000)
        nx, nsx, st = align(e["ea"], e["nea"], o.get("ea"), o.get("nea"), e["beta"], e["se"])
        if st == "mismatch" or nx is None:
            miss.append(rs + "(allele mismatch)")
            continue
        bx.append(nx); bxse.append(nsx); by.append(o["beta"]); byse.append(o["se"]); used.append(rs)
        print("        SNP %-13s exposure b=%+.4f/se=%.4f (%s/%s) | outcome b=%+.4f/se=%.4f (%s/%s) | %s"
              % (rs, e["beta"], e["se"], e["ea"], e["nea"],
                 o["beta"], o["se"], o.get("ea"), o.get("nea"), st))
    if miss:
        print("        Unmatched: %s" % ", ".join(miss[:8]))
    if not bx:
        print("   %-34s no usable SNP" % label)
        return None
    bx = [float(v) for v in bx]; bxse = [float(v) for v in bxse]
    by = [float(v) for v in by]; byse = [float(v) for v in byse]
    if len(bx) == 1:
        b = by[0] / bx[0]; se = abs(byse[0] / bx[0]); meth = "Wald (1 SNP)"
    else:
        w = [1.0 / (s_ ** 2 / b_ ** 2) for b_, s_ in zip(bx, bxse)]
        b = sum(wi * (byi / bxi) for wi, bxi, byi in zip(w, bx, by)) / sum(w)
        se = math.sqrt(1.0 / sum(w)); meth = "IVW (%d SNP)" % len(bx)
    z = b / se if se else float("nan")
    p = math.erfc(abs(z) / math.sqrt(2)) if se else float("nan")
    print("   %-34s %-16s b=%+.4f  SE=%.4f  p=%.4g" % (label, meth, b, se, p))
    return dict(label=label, method=meth, nsnp=len(bx), b=b, se=se, p=p, snps=",".join(used))


def main():
    print("=" * 96)
    print("CRP <-> serum bilirubin  bidirectional cis-MR   %s" % time.strftime("%Y-%m-%d %H:%M"))
    print("=" * 96)
    print()
    print("[0] Dynamically resolve gene coordinates (not from memory)")
    crops = {}
    for sym in ("CRP", "UGT1A1"):
        ens, ch, s, e, nm = gene_coord(sym)
        crops[sym] = (ens, ch, s, e)
        print("   %-8s -> %-18s %-34s hg19:chr%s:%s-%s" % (sym, ens, str(nm)[:34], ch, s, e))
    print()

    print("[1] Build cis instruments")
    crp_ins = get_cis("ebi-a-GCST90029070", crops["CRP"][1], crops["CRP"][2], crops["CRP"][3], "CRP(Said22)")
    crp_ins2 = get_cis("ieu-b-35", crops["CRP"][1], crops["CRP"][2], crops["CRP"][3], "CRP(Ligthart18)")
    bil_ins = get_cis("ebi-a-GCST90025973", crops["UGT1A1"][1], crops["UGT1A1"][2], crops["UGT1A1"][3], "Bil(Barton21)")
    bil_ins2 = get_cis("bbj-a-53", crops["UGT1A1"][1], crops["UGT1A1"][2], crops["UGT1A1"][3], "Bil(BBJ,East Asian)")
    print()

    print("[2] Forward MR: CRP (exposure) -> bilirubin (outcome)")
    o1 = fetch_outcome("ebi-a-GCST90025973", [x["rsid"] for x in crp_ins])
    r1 = run_mr(crp_ins, o1, "Said22 CRP -> Barton21 Bil")
    o1b = fetch_outcome("bbj-a-53", [x["rsid"] for x in crp_ins2])
    r1b = run_mr(crp_ins2, o1b, "Ligthart18 CRP -> BBJ Bil (non-overlapping)")
    # Cross combinations: improve the instrument-outcome pairing rate (different GWASs cover different variants)
    o1c = fetch_outcome("bbj-a-53", [x["rsid"] for x in crp_ins])
    r1c = run_mr(crp_ins, o1c, "[cross] Said22 CRP -> BBJ Bil")
    o1d = fetch_outcome("ebi-a-GCST90025973", [x["rsid"] for x in crp_ins2])
    r1d = run_mr(crp_ins2, o1d, "[cross] Ligthart18 CRP -> Barton21")
    print()

    print("[3] Reverse MR: bilirubin (exposure) -> CRP (outcome)")
    o2 = fetch_outcome("ebi-a-GCST90029070", [x["rsid"] for x in bil_ins])
    r2 = run_mr(bil_ins, o2, "Barton21 Bil -> Said22 CRP")
    o2b = fetch_outcome("ieu-b-35", [x["rsid"] for x in bil_ins2])
    r2b = run_mr(bil_ins2, o2b, "BBJ Bil -> Ligthart18 CRP (non-overlapping)")
    o2c = fetch_outcome("ieu-b-35", [x["rsid"] for x in bil_ins])
    r2c = run_mr(bil_ins, o2c, "[cross] Barton21 Bil -> Ligthart18")
    o2d = fetch_outcome("ebi-a-GCST90029070", [x["rsid"] for x in bil_ins2])
    r2d = run_mr(bil_ins2, o2d, "[cross] BBJ Bil -> Said22 CRP")
    print()

    res = [x for x in (r1, r1b, r1c, r1d, r2, r2b, r2c, r2d) if x]
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, "MR_crp_bilirubin_result.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv = None
        import csv as _csv
        wr = _csv.DictWriter(f, fieldnames=["label", "method", "nsnp", "b", "se", "p", "snps"])
        wr.writeheader(); wr.writerows(res)
    with open(os.path.join(OUT, "MR_crp_bilirubin_instruments.json"), "w", encoding="utf-8") as f:
        json.dump({"CRP_Said22": crp_ins, "CRP_Ligthart18": crp_ins2,
                   "Bil_Barton21": bil_ins, "Bil_BBJ": bil_ins2}, f, ensure_ascii=False, indent=1)
    print("Saved:", p)
    print()
    print("Interpretation: b>0 means higher CRP -> higher bilirubin; b<0 means higher CRP -> lower bilirubin")
    print("      p<0.05 is nominally significant; cis-MR is the primary evidence for causal direction and must be interpreted together with the reverse MR")


if __name__ == "__main__":
    main()
