# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Supplementary analysis: cross-dataset pairing of IL6R / SAA instruments + neonatal outcomes (with proxy SNV lookup)

Motivation: in the main script each exposure uses only its own lead SNP, so some combinations are "missing in the outcome" (different datasets cover different variants).
This script tests the 4 IL6R instruments + 2 SAA instruments **one by one** against every outcome, and tries proxies=1 when a variant is missing.
"""
import os, sys, json, math, time, csv, urllib.request
sys.stdout.reconfigure(encoding="utf-8")

JWT = os.environ.get("OGWAS_JWT", "")
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"
LOG = os.path.join(OUT, "_mr_secondary.log")
_L = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    _L.append(s)
    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(_L))


def og(path, body=None, timeout=180, tries=4):
    for i in range(tries):
        try:
            d = json.dumps(body).encode() if body is not None else None
            rq = urllib.request.Request("https://api.opengwas.io" + path, data=d, headers={
                "Authorization": "Bearer " + JWT, "Content-Type": "application/json",
                "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(rq, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i == tries - 1:
                P("   !! %s" % str(e)[:100]); return None
            time.sleep(3)


# Read the verified instruments from the main script's results (re-fetched here to stay self-consistent)
IL6R = {  # dataset -> (label, N)
    "ebi-a-GCST90012025": "IL6R Folkersen2020 (N=21758)",
    "ebi-a-GCST90019398": "IL6R Pietzner2020 (N=10708)",
    "prot-a-1540": "IL6R Sun2018 INTERVAL (N=3301)",
    "prot-b-23": "IL6R Folkersen2018 (N=3394)",
}
SAA = {
    "ebi-a-GCST90019394": "SAA1 Pietzner2020 (N=10708)",
    "ebi-a-GCST90019411": "SAA2 Pietzner2020 (N=10708)",
}
# lead SNPs obtained empirically in the previous step (including the information needed for effect alignment); re-fetched here to guard against typos
COORD = {"IL6R": ("1", 154377669, 154441926), "SAA": ("11", 18210770, 18341524)}
P_CIS, PAD, PRUNE = 5e-8, 50_000, 500_000


def get_cis(gid, ch, s, e):
    var = "%s:%d-%d" % (ch, s - PAD, e + PAD)
    r = og("/api/associations", {"id": [gid], "variant": [var]})
    rows = r if isinstance(r, list) else []
    sig = sorted([x for x in rows if (x.get("p") or 1) < P_CIS], key=lambda x: x["p"])
    kept = []
    for x in sig:
        if all(abs(x["position"] - k["position"]) > PRUNE for k in kept):
            kept.append(x)
    return kept


def fetch(gid, snps, proxies=0):
    body = {"id": [gid], "variant": list(snps)}
    if proxies:
        body["proxies"] = 1
    r = og("/api/associations", body)
    return r if isinstance(r, list) else []


def wald(e, o, tag=""):
    """Return (b, se, p, note); e/o are the exposure/outcome records"""
    ea, nea = e.get("ea"), e.get("nea")
    oea, onea = o.get("ea"), o.get("nea")
    if ea is None or oea is None:
        return None, None, None, "allele missing"
    by, byse = float(o["beta"]), float(o["se"])
    if ea == onea and nea == oea:
        by, byse = -by, byse
    elif not (ea == oea and nea == onea):
        return None, None, None, "allele mismatch (%s/%s vs %s/%s)" % (ea, nea, oea, onea)
    bx = float(e["beta"])
    b = by / bx; se = abs(byse / bx)
    p = math.erfc(abs(b / se) / math.sqrt(2)) if se else float("nan")
    return b, se, p, "ok"


def main():
    t0 = time.time()
    P("=" * 100)
    P("Supplementary analysis: cross-dataset instrument pairing + neonatal outcomes (with proxy)   %s" % time.strftime("%Y-%m-%d %H:%M"))
    P("=" * 100)
    P("")
    P("[1] Extract all instruments (cis, after pruning)")
    inst = {}
    for gid in IL6R:
        inst["IL6R:" + gid] = get_cis(gid, *COORD["IL6R"])
        for k in inst["IL6R:" + gid]:
            P("   %-22s %-13s p=%.3e beta=%+.4f se=%.4f ea=%s nea=%s F=%.0f"
              % (gid, k["rsid"], k["p"], k["beta"], k["se"], k.get("ea"), k.get("nea"),
                 (k["beta"] / k["se"]) ** 2))
    for gid in SAA:
        inst["SAA:" + gid] = get_cis(gid, *COORD["SAA"])
        for k in inst["SAA:" + gid]:
            P("   %-22s %-13s p=%.3e beta=%+.4f se=%.4f ea=%s nea=%s F=%.0f"
              % (gid, k["rsid"], k["p"], k["beta"], k["se"], k.get("ea"), k.get("nea"),
                 (k["beta"] / k["se"]) ** 2))
    P("")

    BIL = {"ebi-a-GCST90025973": "total bilirubin Barton21 (N=436055)",
           "ukb-d-30840_irnt": "total bilirubin UKB-Neale",
           "bbj-a-53": "bilirubin BBJ East Asian (N=110207)",
           "ebi-a-GCST90025983": "direct bilirubin Barton21 (N=372420)"}
    NEO = {"finn-b-P16_NEONTAL_JAUND_OTH_UNSP_CAUSES": "neonatal jaundice (ncase=133)",
           "finn-b-P16_DISORD_RELATED_GESTATION_LOW_BIRTHWGHTT_NECIFIED": "preterm birth/low birth weight (258)",
           "finn-b-P16_RESPI_DISTR_NEWBO": "neonatal respiratory distress (171)"}
    results = []

    P("=" * 100)
    P("[2] * Cross-dataset pairing: every IL6R/SAA instrument x every bilirubin outcome (tries proxies when missing)")
    P("=" * 100)
    for key, ins in inst.items():
        if not ins:
            continue
        P("")
        P(" > %s   instruments %s" % (key, ",".join(x["rsid"] for x in ins)))
        for gid, lab in BIL.items():
            recs = fetch(gid, [x["rsid"] for x in ins])
            got = {r["rsid"]: r for r in recs}
            hit = False
            for e in ins:
                o = got.get(e["rsid"])
                if o is None:
                    continue
                b, se, p, note = wald(e, o)
                if b is None:
                    P("        %-13s %-30s %s" % (e["rsid"], lab, note)); continue
                hit = True
                P("        %-13s %-30s b=%+.4f SE=%.4f p=%.4g%s" % (e["rsid"], lab, b, se, p,
                                                                   "  **" if p < 0.05 else ""))
                results.append(dict(exposure=key, snp=e["rsid"], outcome=lab, b=b, se=se, p=p,
                                    method="Wald", proxy=""))
            if not hit:
                # try a proxy
                recs = fetch(gid, [x["rsid"] for x in ins], proxies=1)
                seen = {}
                for r in recs:
                    t = r.get("target_snp") or r.get("target")
                    if t and r.get("proxy") in ("yes", True, 1):
                        seen.setdefault(t, r)
                if seen:
                    P("        (direct match missing, %d proxy hit(s))" % len(seen))
                    for e in ins:
                        pr = seen.get(e["rsid"])
                        if not pr:
                            continue
                        b, se, p, note = wald(e, pr)
                        if b is None:
                            continue
                        P("        %-13s %-30s b=%+.4f SE=%.4f p=%.4g  [proxy %s r2=%s]"
                          % (e["rsid"], lab, b, se, p, pr.get("rsid"), pr.get("r2")))
                        results.append(dict(exposure=key, snp=e["rsid"], outcome=lab, b=b, se=se, p=p,
                                            method="Wald+proxy", proxy=str(pr.get("rsid"))))
                else:
                    P("        %-13s %-30s missing and no proxy" % ("all", lab))
    P("")

    P("=" * 100)
    P("[3] Neonatal outcomes (FinnGen): all IL6R / SAA instruments paired one by one (with proxy)")
    P("=" * 100)
    for key, ins in inst.items():
        if not ins:
            continue
        P("")
        P(" > %s" % key)
        for gid, lab in NEO.items():
            recs = fetch(gid, [x["rsid"] for x in ins])
            got = {r["rsid"]: r for r in recs}
            any_hit = False
            for e in ins:
                o = got.get(e["rsid"])
                if o is None:
                    continue
                b, se, p, note = wald(e, o)
                if b is None:
                    continue
                any_hit = True
                P("        %-13s %-26s b=%+.4f SE=%.4f p=%.4g%s" % (e["rsid"], lab, b, se, p,
                                                                     "  **" if p < 0.05 else ""))
                results.append(dict(exposure=key, snp=e["rsid"], outcome=lab, b=b, se=se, p=p,
                                    method="Wald", proxy=""))
            if not any_hit:
                recs = fetch(gid, [x["rsid"] for x in ins], proxies=1)
                seen = {}
                for r in recs:
                    t = r.get("target_snp") or r.get("target")
                    if t and r.get("proxy") in ("yes", True, 1):
                        seen.setdefault(t, r)
                if seen:
                    for e in ins:
                        pr = seen.get(e["rsid"])
                        if not pr:
                            continue
                        b, se, p, note = wald(e, pr)
                        if b is None:
                            continue
                        P("        %-13s %-26s b=%+.4f SE=%.4f p=%.4g  [proxy %s r2=%s]"
                          % (e["rsid"], lab, b, se, p, pr.get("rsid"), pr.get("r2")))
                        results.append(dict(exposure=key, snp=e["rsid"], outcome=lab, b=b, se=se, p=p,
                                            method="Wald+proxy", proxy=str(pr.get("rsid"))))
                else:
                    P("        %-13s %-26s missing and no proxy (cannot be evaluated)" % ("all instruments", lab))
    P("")

    p = os.path.join(OUT, "MR_il6_saa_secondary_result.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["exposure", "snp", "outcome", "b", "se", "p", "method", "proxy"])
        w.writeheader(); w.writerows(results)
    P("%d estimates in total -> %s" % (len(results), p))
    P("Total elapsed %.1f min" % ((time.time() - t0) / 60))


if __name__ == "__main__":
    main()
