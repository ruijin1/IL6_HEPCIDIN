# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Download the remaining datasets needed for the mechanism pillar (URLs corrected per GEO directory rules: stem = GSE + first 2 digits + nnn)

* Lesson: GEO FTP directory names are `GSE` + some leading digits of the number + `nnn` --
   GSE25504 -> `GSE25nnn` (**not** GSE255nnn); GSE218242 -> `GSE218nnn`.
   Getting this wrong only yields a 404, which is easily misread as "the file does not exist".
"""
import os, sys, gzip, urllib.request, urllib.error
sys.stdout.reconfigure(encoding="utf-8")

DST = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
os.makedirs(DST, exist_ok=True)

B = "https://ftp.ncbi.nlm.nih.gov/geo/series/"
JOBS = [
    # GSE25504 neonatal sepsis whole blood (merged matrix)
    (B + "GSE25nnn/GSE25504/suppl/GSE25504_fixed.txt.gz", "GSE25504_fixed.txt.gz"),
    # GSE218242 primary hepatocytes: circadian rhythm x inflammation x drug metabolism (human+mouse, n=49) -- includes TPM and sample annotation
    (B + "GSE218nnn/GSE218242/suppl/GSE218242_190712Bha_set1_l2tpm.txt.gz", "GSE218242_l2tpm.txt.gz"),
    (B + "GSE218nnn/GSE218242/matrix/GSE218242_series_matrix.txt.gz", "GSE218242_series_matrix.txt.gz"),
    # GSE218243 second batch of the same study (n=16)
    (B + "GSE218nnn/GSE218243/suppl/GSE218243_210308Bha_set3_l2tpm.txt.gz", "GSE218243_l2tpm.txt.gz"),
    (B + "GSE218nnn/GSE218243/matrix/GSE218243_series_matrix.txt.gz", "GSE218243_series_matrix.txt.gz"),
]


def remote_size(url, to=40):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=to) as r:
            return int(r.headers.get("Content-Length") or 0)
    except Exception as e:
        print("  HEAD failed %s : %s" % (url.split("/")[-1], str(e)[:60])); return 0


def gzip_ok(fp):
    try:
        with gzip.open(fp, "rb") as f:
            while f.read(1 << 22):
                pass
        return True
    except Exception:
        return False


for url, name in JOBS:
    fp = os.path.join(DST, name)
    rs = remote_size(url)
    ls = os.path.getsize(fp) if os.path.exists(fp) else 0
    if rs and ls == rs and (not name.endswith(".gz") or gzip_ok(fp)):
        print("  OK(already exists) %-42s %.2f MB" % (name, ls / 2**20)); continue
    if rs == 0:
        print("  skipping (remote unreachable) %-42s" % name); continue
    print("  downloading %-42s (%.2f MB) ..." % (name, rs / 2**20), end=" ", flush=True)
    try:
        tmp = fp + ".part"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b)
        os.replace(tmp, fp)
        ls2 = os.path.getsize(fp)
        ok = (ls2 == rs) and (not name.endswith(".gz") or gzip_ok(fp))
        print("%.2f MB  size%s gzip%s" % (ls2 / 2**20, "✓" if ls2 == rs else "✗",
                                          "✓" if (not name.endswith(".gz") or gzip_ok(fp)) else "✗"))
    except Exception as e:
        print("failed:", str(e)[:110])

print()
print("Final:")
for _, name in JOBS:
    fp = os.path.join(DST, name)
    if os.path.exists(fp):
        print("  %-42s %.2f MB" % (name, os.path.getsize(fp) / 2**20))
