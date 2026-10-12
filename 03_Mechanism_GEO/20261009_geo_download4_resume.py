# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Resumable downloader (GEO FTP silently hangs, so progress + chunked retry are mandatory)
* Lesson: a single urllib connection reading a large GEO file stalls halfway without raising an error (this time GSE275342_RAW.tar stalled at 1.00 MB for over 10 minutes),
   so we switched to: range-chunked requests + per-chunk retry + real-time logging, and judge **progress per chunk** rather than as a whole.
"""
import os, sys, time, urllib.request
sys.stdout.reconfigure(encoding="utf-8")

DST = r"D:\workbuddy_workspace\03-data\MR_finngen\geo"
B = "https://ftp.ncbi.nlm.nih.gov/geo/series/"
JOBS = [
    (B + "GSE275nnn/GSE275342/suppl/GSE275342_RAW.tar", "GSE275342_RAW.tar"),
    (B + "GSE102nnn/GSE102006/matrix/GSE102006_series_matrix.txt.gz", "GSE102006_series_matrix.txt.gz"),
]
CHUNK = 1 << 20  # 1 MB


def size(url):
    for _ in range(4):
        try:
            r = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(r, timeout=45) as resp:
                return int(resp.headers.get("Content-Length") or 0)
        except Exception as e:
            print("   HEAD retry (%s)" % str(e)[:50]); time.sleep(3)
    return 0


def fetch(url, dst):
    total = size(url)
    if not total:
        print("  remote unreachable, giving up %s" % os.path.basename(dst)); return False
    have = os.path.getsize(dst) if os.path.exists(dst) else 0
    if have == total:
        print("  already complete %s %.2f MB" % (os.path.basename(dst), total / 2**20)); return True
    if have > total:
        os.remove(dst); have = 0
    mode = "ab" if have else "wb"
    print("  starting %s  total %.2f MB, have %.2f MB" % (os.path.basename(dst), total / 2**20, have / 2**20))
    fails = 0
    with open(dst, mode) as f:
        while have < total:
            # * per-chunk timeout guard: GEO hangs silently, so reconnect if nothing is read
            end = min(have + CHUNK, total) - 1
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0",
                                                           "Range": "bytes=%d-%d" % (have, end)})
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = r.read()
                if not data:
                    raise IOError("empty response")
                f.write(data); f.flush()
                have += len(data)
                fails = 0
                if have % (8 * CHUNK) < CHUNK or have == total:
                    print("    %5.1f%%  %.2f / %.2f MB" % (100 * have / total, have / 2**20, total / 2**20))
            except Exception as e:
                fails += 1
                print("    chunk failed (#%d) %s; reconnecting in 3 s" % (fails, str(e)[:60]))
                if fails > 25:
                    print("    too many consecutive failures, aborting (the next run will resume automatically)"); return False
                time.sleep(3)
                try:
                    have = os.path.getsize(dst)
                except OSError:
                    have = 0
    ok = os.path.getsize(dst) == total
    print("  %s %s %.2f MB" % ("complete" if ok else "incomplete", os.path.basename(dst), os.path.getsize(dst) / 2**20))
    return ok


def main():
    print("Resumable download  %s" % DST)
    for url, name in JOBS:
        fp = os.path.join(DST, name)
        for attempt in (1, 2, 3):
            if os.path.exists(fp) and os.path.exists(fp + ".part"):
                break
            if fetch(url, fp):
                break
            print("  attempt %d incomplete, retrying..." % attempt)
    print("\nFinal:")
    for _, name in JOBS:
        fp = os.path.join(DST, name)
        print("  %-42s %s" % (name, "%.2f MB" % (os.path.getsize(fp) / 2**20) if os.path.exists(fp) else "missing"))


if __name__ == "__main__":
    main()
