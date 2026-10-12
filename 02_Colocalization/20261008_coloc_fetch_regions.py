# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Fetch regional summary statistics for the IL6R / SAA loci (for coloc)

Regions: IL6R +/-150kb (hg19 chr1:154,227,669-154,591,926)
         SAA1/SAA2 +/-150kb (hg19 chr11:18,110,770-18,441,524)
One CSV per dataset, fields consistent with what coloc requires.
WARNING: the token is read from the environment variable OGWAS_JWT and is never written to disk.
"""
import os, sys, csv, time, json, urllib.request
sys.stdout.reconfigure(encoding="utf-8")

JWT = os.environ.get("OGWAS_JWT", "")
OUT = r"D:\workbuddy_workspace\03-data\MR_finngen\coloc"
LOG = os.path.join(OUT, "_fetch.log")
_L = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    _L.append(s)
    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(_L))


def og(path, body=None, timeout=240, tries=4):
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
                P("   !! %s" % str(e)[:110]); return None
            time.sleep(3)


os.makedirs(OUT, exist_ok=True)
open(LOG, "w").close()

REGIONS = {
    "IL6R": "1:154227669-154591926",
    "SAA": "11:18110770-18441524",
}
# key -> (dataset_id, region, note)
JOBS = {
    "exp_IL6R_folkersen2018": ("prot-b-23", "IL6R", "IL6R pQTL Folkersen2018 N=3394"),
    "exp_IL6R_folkersen2020": ("ebi-a-GCST90012025", "IL6R", "IL6R pQTL Folkersen2020 N=21758"),
    "exp_IL6R_pietzner2020": ("ebi-a-GCST90019398", "IL6R", "IL6R pQTL Pietzner2020 N=10708"),
    "exp_SAA1_pietzner2020": ("ebi-a-GCST90019394", "SAA", "SAA1 pQTL Pietzner2020 N=10708"),
    "out_bil_barton2021": ("ebi-a-GCST90025973", "IL6R", "Total bilirubin Barton2021 N=436055"),
    "out_bil_barton2021_SAAregion": ("ebi-a-GCST90025973", "SAA", "Total bilirubin Barton2021 (SAA region)"),
    "out_bil_ukb_neale": ("ukb-d-30840_irnt", "IL6R", "Total bilirubin UKB-Neale"),
    "out_bil_ukb_neale_SAAregion": ("ukb-d-30840_irnt", "SAA", "Total bilirubin UKB-Neale (SAA region)"),
    "out_bil_bbj": ("bbj-a-53", "IL6R", "Bilirubin BBJ East Asian N=110207"),
    "out_crp_said2022": ("ebi-a-GCST90029070", "IL6R", "CRP Said2022 N=575531 (positive control)"),
}

FIELDS = ["rsid", "chr", "position", "ea", "nea", "eaf", "beta", "se", "p", "n"]
summary = {}
for key, (gid, reg, note) in JOBS.items():
    var = REGIONS[reg]
    r = og("/api/associations", {"id": [gid], "variant": [var]})
    rows = r if isinstance(r, list) else []
    keep = []
    for x in rows:
        if not x.get("rsid") or x.get("beta") is None or x.get("se") is None:
            continue
        if x.get("ea") is None or x.get("nea") is None:
            continue
        keep.append({k: x.get(k) for k in FIELDS})
    fp = os.path.join(OUT, key + ".csv")
    with open(fp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader(); w.writerows(keep)
    ns = [x["n"] for x in keep if x.get("n")]
    P("   %-30s %-46s %s  variants %-5d\n        region %s  N=%s  min p=%.3g"
      % (key, gid, note, len(keep), var, (max(ns) if ns else "?"),
         min((x["p"] for x in keep if x.get("p")), default=float("nan"))))
    summary[key] = dict(dataset=gid, region=reg, note=note, nvar=len(keep),
                        n=(max(ns) if ns else None), file=fp)
    time.sleep(0.5)

json.dump(summary, open(os.path.join(OUT, "_summary.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
P("")
P("Done: %d datasets -> %s" % (len(summary), OUT))
