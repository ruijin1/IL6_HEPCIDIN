# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Download the **independent replication** datasets for the mechanism pillar (not the two circadian ones)

Background and findings (2026-10-09):
  Although GSE218242 / GSE218243 belong to the "primary hepatocytes, inflammation x drug metabolism" superseries GSE218244 (PMID 38657074),
  they have **no IL-6 stimulation arm**:
    · GSE218242 = circadian time course (24-72h, every 3h, n=49), single factor "time"
    · GSE218243 = ARNTL(BMAL1) siRNA vs NT control (36/42/48/54h, 2 replicates each, n=16)
  -> cannot directly replicate IL-6 -> UGT1A1. Instead we download three human hepatocyte datasets that **do have an inflammatory stimulation arm**:

  GSE230325  dHepaRG (differentiated human liver progenitor cell line) LPS 0/6/24h, 3 replicates -- human hepatocyte acute-phase response
  GSE275342  Hep3B + Huh7 cell lines, IL-6-induced RNA-seq -- direct human IL-6 stimulation
  GSE102006  HepG2 + cytokine mix (TNF/IFNγ/IL-1α/IL-6/LPS) 6/12/24h -- large-sample human hepatocyte inflammation

* Directory name rule (already learned the hard way): GSE number -> `GSE` + first 3 digits + `nnn`
"""
import os, sys, gzip, urllib.request
sys.stdout.reconfigure(encoding="utf-8")

DST = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
os.makedirs(DST, exist_ok=True)
B = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

JOBS = [
    # ---- GSE230325 human dHepaRG + LPS acute phase (first choice: has normalized counts) ----
    (B + "GSE230nnn/GSE230325/suppl/GSE230325_normalized_counts_file_DESeq2_output.csv.gz",
     "GSE230325_normcounts.csv.gz"),
    (B + "GSE230nnn/GSE230325/matrix/GSE230325_series_matrix.txt.gz",
     "GSE230325_series_matrix.txt.gz"),
    # ---- GSE275342 human Hep3B/Huh7 + IL-6 (RAW.tar contains per-sample counts) ----
    (B + "GSE275nnn/GSE275342/matrix/GSE275342_series_matrix.txt.gz",
     "GSE275342_series_matrix.txt.gz"),
    (B + "GSE275nnn/GSE275342/suppl/GSE275342_RAW.tar",
     "GSE275342_RAW.tar"),
    # ---- GSE102006 human HepG2 + cytokine mix (array; take series matrix only) ----
    (B + "GSE102nnn/GSE102006/matrix/GSE102006_series_matrix.txt.gz",
     "GSE102006_series_matrix.txt.gz"),
]


def remote_size(url, to=40):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=to) as r:
            return int(r.headers.get("Content-Length") or 0)
    except Exception as e:
        print("  HEAD failed %s : %s" % (url.split("/")[-1], str(e)[:60]))
        return 0


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
    isgz = name.endswith(".gz")
    if rs and ls == rs and (not isgz or gzip_ok(fp)):
        print("  OK(already exists) %-42s %.2f MB" % (name, ls / 2**20))
        continue
    if rs == 0:
        print("  skipping (remote unreachable) %-42s" % name)
        continue
    print("  downloading %-42s (%.2f MB) ..." % (name, rs / 2**20), end=" ", flush=True)
    try:
        tmp = fp + ".part"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=900) as r, open(tmp, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b)
        os.replace(tmp, fp)
        ls2 = os.path.getsize(fp)
        print("%.2f MB  size%s gzip%s" % (ls2 / 2**20, "✓" if ls2 == rs else "✗",
                                          "✓" if (not isgz or gzip_ok(fp)) else "✗"))
    except Exception as e:
        print("failed:", str(e)[:110])

print("\nFinal:")
for _, name in JOBS:
    fp = os.path.join(DST, name)
    if os.path.exists(fp):
        print("  %-42s %.2f MB" % (name, os.path.getsize(fp) / 2**20))
