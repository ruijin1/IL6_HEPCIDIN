# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Batch download of FinnGen R12 neonatal summary statistics
  - Resumable downloads (HTTP Range); can be re-run and skip completed files
  - Each file is checked for integrity against the official Content-Length
  - Measured bandwidth ~113 KB/s (~2 hours per file), so files are ordered by relevance and downloaded serially
"""
import os, sys, time, json, urllib.request, urllib.error

sys.stdout.reconfigure(encoding="utf-8")

OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"
BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/"
LOG = os.path.join(OUT, "_download_progress.log")
os.makedirs(OUT, exist_ok=True)

# Ordered by relevance to the new "inflammation -> bilirubin" direction
TARGETS = [
    ("P16_NEONTAL_JAUND_OTH_UNSP_CAUSES",              "neonatal jaundice (ncase 499)"),
    ("P16_BACTERIAL_SEPSIS_NEWBO",                     "neonatal bacterial sepsis (ncase 203)"),
    ("P16_DISORD_RELATED_GESTATION_LOW_BIRTHWGHTT_NECIFIED", "preterm birth/low birth weight (ncase 782)"),
    ("P16_RESPI_DISTR_NEWBO",                          "neonatal respiratory distress (ncase 498)"),
    ("P16_CHRONIC_RESPI_DISEA_ORIGINA_PERINA_PERIOD",  "chronic respiratory disease of perinatal origin = BPD (ncase 69)"),
]


def log(msg):
    line = "%s  %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return int(r.headers.get("Content-Length", 0))


def fetch(ph, label):
    url = BASE + "finngen_R12_%s.gz" % ph
    dst = os.path.join(OUT, "finngen_R12_%s.gz" % ph)
    try:
        total = remote_size(url)
    except Exception as e:
        log("  !! %s failed to get size: %s" % (ph, str(e)[:80]))
        return False
    if total <= 0:
        log("  !! %s remote size is 0" % ph)
        return False
    have = os.path.getsize(dst) if os.path.exists(dst) else 0
    if have == total:
        log("  [OK] %s already complete (%.1f MB), skipping" % (label, total / 2 ** 20))
        return True
    if have > total:
        log("  ! %s local file is larger than remote, deleting and re-downloading" % ph)
        os.remove(dst); have = 0
    log("  down %s  start %.1f/%.1f MB" % (label, have / 2 ** 20, total / 2 ** 20))
    t0 = time.time(); last = t0; last_n = have
    while have < total:
        try:
            hdr = {"User-Agent": "Mozilla/5.0"}
            if have:
                hdr["Range"] = "bytes=%d-" % have
            req = urllib.request.Request(url, headers=hdr)
            with urllib.request.urlopen(req, timeout=60) as r:
                mode = "ab" if have else "wb"
                with open(dst, mode) as f:
                    while True:
                        b = r.read(1 << 18)
                        if not b:
                            break
                        f.write(b); have += len(b)
                        now = time.time()
                        if now - last >= 20:
                            sp = (have - last_n) / (now - last)
                            eta = (total - have) / sp / 60 if sp > 0 else 0
                            log("     %5.1f%%  %.0f KB/s  ~%.0f minutes remaining"
                                % (100 * have / total, sp / 1024, eta))
                            last, last_n = now, have
        except Exception as e:
            log("     interrupted (%s), resuming in 5 s" % str(e)[:60])
            time.sleep(5)
    ok = os.path.getsize(dst) == total
    log("  %s %s  elapsed %.1f min  size %.1f MB"
        % ("[OK]" if ok else "[FAIL]", label, (time.time() - t0) / 60, os.path.getsize(dst) / 2 ** 20))
    return ok


def main():
    log("=" * 70)
    log("Starting FinnGen R12 neonatal summary statistics download (%d files, serial)" % len(TARGETS))
    log("=" * 70)
    done = []
    for ph, label in TARGETS:
        if fetch(ph, label):
            done.append(ph)
    log("")
    log("All done: %d / %d succeeded" % (len(done), len(TARGETS)))
    for ph, label in TARGETS:
        p = os.path.join(OUT, "finngen_R12_%s.gz" % ph)
        if os.path.exists(p):
            log("   %-58s %8.1f MB" % (label, os.path.getsize(p) / 2 ** 20))


if __name__ == "__main__":
    main()
