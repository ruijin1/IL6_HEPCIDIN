# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Mechanism pillar step 2: inspect the available files (matrix / supplementary) and sample structure of candidate GEO datasets"""
import sys, re, json, urllib.request, urllib.error
sys.stdout.reconfigure(encoding="utf-8")

CAND = {
    "GSE202045": "IL-6-stimulated hepatocytes (n=108, human+mouse)",
    "GSE218242": "circadian rhythm / inflammation / drug metabolism primary hepatocytes (n=49)",
    "GSE218243": "same as above, previous batch (n=16)",
    "GSE25504":  "neonatal sepsis whole blood mRNA (n=170)",
    "GSE189320": "human in vitro inflammatory liver model (n=35)",
}


def head(url, to=30):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=to) as r:
            return r.status, int(r.headers.get("Content-Length") or 0)
    except urllib.error.HTTPError as e:
        return e.code, 0
    except Exception:
        return 0, 0


def listing(url, to=30):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=to) as r:
            html = r.read().decode("utf-8", "replace")
    except Exception as e:
        return None, str(e)[:80]
    files = re.findall(r'href="([^"]+)"', html)
    return [f for f in files if not f.startswith("?") and f not in ("/geo/", "/geo/series/")], None


BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/"
out = {}
for acc, note in CAND.items():
    stem, nn = acc[:-3] + "nnn", acc
    print("=" * 100)
    print("[%s] %s" % (acc, note))
    print("=" * 100)
    rec = {"note": note}
    for sub in ("matrix", "suppl"):
        url = "%s%s/%s/" % (BASE, stem, nn if sub == "matrix" else nn)
        if sub == "matrix":
            url = "%s%s/%s/matrix/" % (BASE, stem, nn)
        else:
            url = "%s%s/%s/suppl/" % (BASE, stem, nn)
        fl, err = listing(url)
        if err:
            print("   [%s] directory listing failed: %s" % (sub, err)); continue
        if not fl:
            print("   [%s] no files" % sub); continue
        print("   [%s]" % sub)
        for f in fl:
            if not f.endswith((".gz", ".txt", ".csv", ".tsv", ".xlsx", ".zip", ".tar")):
                continue
            size = 0
            st, size = head(url + f)
            print("      %-58s %8.2f MB  HTTP %s" % (f, size / 2**20, st))
            rec.setdefault(sub, []).append(dict(file=f, size=size))
    out[acc] = rec
    print()

json.dump(out, open(r"D:\workbuddy_workspace\03-data\MR_finngen\geo_files.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("saved geo_files.json")
