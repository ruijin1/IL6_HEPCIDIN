# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Mechanism pillar step 3: download GEO processed files and check sample naming (to determine grouping)"""
import os, gzip, sys, io, urllib.request
sys.stdout.reconfigure(encoding="utf-8")

DST = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
os.makedirs(DST, exist_ok=True)

FILES = [
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202045/suppl/"
     "GSE202045_Lukowski_2015_HumanPrimaryHep_IL6TC_mRNA_rawCounts.csv.gz",
     "GSE202045_humanPrimaryHep_mRNA_counts.csv.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202045/suppl/"
     "GSE202045_Lukowski_2015_HepG2_IL6TC_mRNA_rawCounts.csv.gz",
     "GSE202045_HepG2_mRNA_counts.csv.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202045/suppl/"
     "GSE202045_Lukowski_2015_MousePrimaryHep_IL6TC_mRNA_rawCounts.csv.gz",
     "GSE202045_mousePrimaryHep_mRNA_counts.csv.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE25504/suppl/GSE25504_fixed.txt.gz",
     "GSE25504_fixed.txt.gz"),
]

for url, name in FILES:
    fp = os.path.join(DST, name)
    if os.path.exists(fp) and os.path.getsize(fp) > 1000:
        print("already exists, skipping:", name, os.path.getsize(fp))
        continue
    print("downloading %-52s ..." % name, end=" ", flush=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=300) as r, open(fp, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b)
        print("OK %.2f MB" % (os.path.getsize(fp) / 2**20))
    except Exception as e:
        print("failed:", str(e)[:120])

print()
print("=" * 104)
print("Check header / sample names / row counts for each file")
print("=" * 104)
for name in [n for _, n in FILES]:
    fp = os.path.join(DST, name)
    if not os.path.exists(fp):
        continue
    print()
    print("########## %s (%.2f MB) ##########" % (name, os.path.getsize(fp) / 2**20))
    with gzip.open(fp, "rt", encoding="utf-8", errors="replace") as f:
        head = f.readline().rstrip("\n")
        cols = head.split("\t") if "\t" in head else head.split(",")
        print("  fields in header row: %d" % len(cols))
        print("  first 24 fields: %s" % cols[:24])
        n = 0
        for _ in f:
            n += 1
        print("  data rows: %s" % format(n, ","))
