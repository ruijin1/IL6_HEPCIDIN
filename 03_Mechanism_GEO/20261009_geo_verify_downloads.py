# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Verify GEO download integrity (compare Content-Length) and automatically re-download incomplete files.

Why this is needed: a real trap hit on 2026-10-08 -- when a foreground command is killed by SIGTERM it leaves behind a file that **appears to exist but is truncated**,
   and the "if the file exists and is >1000 bytes, skip it" logic would treat a truncated file as already downloaded.
   For any intermediate data written to disk from the network, always verify with Content-Length (or a gzip integrity check).
"""
import os, json, sys, gzip, urllib.request, urllib.error
sys.stdout.reconfigure(encoding="utf-8")

DST = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
URLS = {
    "GSE202045_humanPrimaryHep_mRNA_counts.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202045/suppl/GSE202045_Lukowski_2015_HumanPrimaryHep_IL6TC_mRNA_rawCounts.csv.gz",
    "GSE202045_HepG2_mRNA_counts.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202045/suppl/GSE202045_Lukowski_2015_HepG2_IL6TC_mRNA_rawCounts.csv.gz",
    "GSE202045_mousePrimaryHep_mRNA_counts.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202045/suppl/GSE202045_Lukowski_2015_MousePrimaryHep_IL6TC_mRNA_rawCounts.csv.gz",
    "GSE25504_fixed.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE25504/suppl/GSE25504_fixed.txt.gz",
}


def remote_size(url, to=40):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=to) as r:
            return int(r.headers.get("Content-Length") or 0)
    except Exception as e:
        print("   HEAD failed:", str(e)[:80]); return 0


def gzip_ok(fp):
    """Decompress in chunks to verify gzip integrity (checking only the trailing CRC misses mid-file truncation)"""
    try:
        with gzip.open(fp, "rb") as f:
            while True:
                if not f.read(1 << 22):
                    break
        return True
    except Exception as e:
        return False


def download(url, fp):
    tmp = fp + ".part"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as f:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            f.write(b)
    os.replace(tmp, fp)


print("=" * 100)
print("GEO download integrity check and repair")
print("=" * 100)
for name, url in URLS.items():
    fp = os.path.join(DST, name)
    rs = remote_size(url)
    ls = os.path.getsize(fp) if os.path.exists(fp) else 0
    status = "OK" if (rs and ls == rs) else ("missing" if ls == 0 else "**incomplete / size mismatch**")
    print("  %-46s local %9d / remote %9d  %s" % (name, ls, rs, status))
    need = (rs == 0) or (ls != rs)
    if need:
        print("      -> re-downloading ...", end=" ", flush=True)
        try:
            os.path.exists(fp) and os.remove(fp)
            download(url, fp)
            ls2 = os.path.getsize(fp)
            gz = gzip_ok(fp) if name.endswith(".gz") else True
            ok = (rs == 0 or ls2 == rs) and gz
            print("%.2f MB  size%s  gzip%s" % (ls2 / 2**20, "✓" if (rs == 0 or ls2 == rs) else "✗",
                                               "✓" if gz else "✗"))
            if not ok:
                print("      !! still incomplete, please re-run this script")
        except Exception as e:
            print("failed:", str(e)[:120])
    else:
        gz = gzip_ok(fp) if name.endswith(".gz") else True
        print("      gzip integrity: %s" % ("✓" if gz else "✗ re-download needed"))

print()
print("Done. Final files:")
for name in URLS:
    fp = os.path.join(DST, name)
    if os.path.exists(fp):
        print("  %-46s %.2f MB" % (name, os.path.getsize(fp) / 2**20))
