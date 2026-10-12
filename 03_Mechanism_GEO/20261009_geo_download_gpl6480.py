# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Download the GPL6480 platform annotation (GSE102006 = HepG2 + cytokine mix, Agilent 4x44K)
Take only the GENE_SYMBOL column for the probe -> symbol mapping.
"""
import os, sys, time, urllib.request
sys.stdout.reconfigure(encoding="utf-8")

DST = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz"
DSTF = os.path.join(DST, "GPL6480.annot.gz")
CHUNK = 1 << 20


def size(url):
    for _ in range(4):
        try:
            r = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(r, timeout=45) as resp:
                return int(resp.headers.get("Content-Length") or 0)
        except Exception as e:
            print("  HEAD retry (%s)" % str(e)[:50]); time.sleep(3)
    return 0


def fetch(url, dst):
    total = size(url)
    if not total:
        print("  remote unreachable"); return False
    have = os.path.getsize(dst) if os.path.exists(dst) else 0
    if have == total:
        print("  already complete %.2f MB" % (total / 2**20)); return True
    if have > total:
        os.remove(dst); have = 0
    print("  starting %.2f MB, have %.2f MB" % (total / 2**20, have / 2**20))
    fails = 0
    with open(dst, "ab" if have else "wb") as f:
        while have < total:
            end = min(have + CHUNK, total) - 1
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0",
                                                           "Range": "bytes=%d-%d" % (have, end)})
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = r.read()
                if not data:
                    raise IOError("empty response")
                f.write(data); f.flush(); have += len(data); fails = 0
                if have % (5 * CHUNK) < CHUNK or have == total:
                    print("    %5.1f%%  %.2f / %.2f MB" % (100 * have / total, have / 2**20, total / 2**20))
            except Exception as e:
                fails += 1
                print("    chunk failed (#%d) %s" % (fails, str(e)[:60]))
                if fails > 25:
                    return False
                time.sleep(3)
                try:
                    have = os.path.getsize(dst)
                except OSError:
                    have = 0
    print("  complete %s" % ("✓" if os.path.getsize(dst) == total else "incomplete"))
    return os.path.getsize(dst) == total


if __name__ == "__main__":
    fetch(URL, DSTF)
