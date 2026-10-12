# -*- coding: utf-8 -*-
r"""FinnGen R12 tabix region-fetch utility (pure-Python implementation for environments without the tabix executable)

Key points (all from the tested recipe of the skill `eqtl-catalogue-region-fetch`):
  1. `.tbi` itself is BGZF-compressed -> it must be decompressed block by block first (do not use the gzip module: it raises EOFError on truncation and loses the buffer)
  2. Index structure: magic(4) + 8 x int32 + names -> per ref: n_bin + bins(chunks) + n_intv + linear index
  3. reg2bins: min_shift=14, n_lvls=5, s = 14+10+5 = 29, per level t += 1<<((l<<1)+l)
  4. virtual offset = (block_start << 16) | within_block_offset
     -> file offset = voffset >> 16; when fetching an interval the end must be padded by +65536 to cover the block boundary
  5. After fetching with HTTP Range, **decompress BGZF block by block**; line endings are CRLF, so the last column retains a stray \r
"""
import zlib, struct, urllib.request


def http_range(url, start, end, timeout=300, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={
                "Range": "bytes=%d-%d" % (start, end), "User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            if i == tries - 1:
                raise
    return b""


def bgzf_decompress(raw, skip_broken_head=True):
    """Decompress BGZF block by block (a plain gzip stream is safely truncated here)"""
    i, n, out = 0, len(raw), []
    while i + 18 <= n:
        if raw[i:i + 2] != b"\x1f\x8b":
            if skip_broken_head and not out:
                i += 1
                continue
            break
        xlen = int.from_bytes(raw[i + 10:i + 12], "little")
        bsize, j, e = None, i + 12, i + 12 + xlen
        while j + 4 <= e:
            if raw[j] == 66 and raw[j + 1] == 67:          # 'BC' subfield
                bsize = int.from_bytes(raw[j + 4:j + 6], "little")
                break
            j += 4 + int.from_bytes(raw[j + 2:j + 4], "little")
        if bsize is None:
            break
        total = bsize + 1
        if i + total > n:
            break
        try:
            out.append(zlib.decompress(raw[i + 12 + xlen:i + total - 8], -15))
        except Exception:
            break
        i += total
    return b"".join(out)


def read_tbi(tbi_url, cache_path=None):
    import os
    if cache_path and os.path.exists(cache_path):
        data = open(cache_path, "rb").read()
    else:
        raw = http_range(tbi_url, 0, 8 << 20)
        data = bgzf_decompress(raw)
        if cache_path:
            open(cache_path, "wb").write(data)
    if data[:4] != b"TBI\x01":
        raise ValueError("not a valid .tbi (magic=%r)" % data[:4])
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", data[4:36])
    off = 36
    names = data[off:off + l_nm].split(b"\x00")
    names = [x.decode() for x in names if x]
    off += l_nm
    refs = {}
    for k in range(n_ref):
        n_bin, = struct.unpack("<i", data[off:off + 4]); off += 4
        bins = {}
        for _ in range(n_bin):
            bid, nc = struct.unpack("<Ii", data[off:off + 8]); off += 8
            ch = []
            for _ in range(nc):
                cb, ce = struct.unpack("<QQ", data[off:off + 16]); off += 16
                ch.append((cb, ce))
            bins[bid] = ch
        n_intv, = struct.unpack("<i", data[off:off + 4]); off += 4
        lin = struct.unpack("<%dQ" % n_intv, data[off:off + 8 * n_intv]); off += 8 * n_intv
        refs[names[k]] = dict(bins=bins, lin=lin)
    return dict(n_ref=n_ref, col_seq=col_seq, col_beg=col_beg, col_end=col_end,
                skip=skip, refs=refs)


def reg2bins(beg, end, min_shift=14, n_lvls=5):
    """htslib reg2bins: 0-based half-open"""
    lst = []
    end -= 1
    s = min_shift + (n_lvls << 1) + n_lvls          # 14 + 10 + 5 = 29
    t = 0
    for l in range(n_lvls + 1):
        lst.extend(range(t + (beg >> s), t + (end >> s) + 1))
        s -= 3
        t += 1 << ((l << 1) + l)
    return lst


def _merge_chunks(chunks):
    chunks = sorted(chunks)
    out = []
    for s, e in chunks:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def fetch_region(data_url, tbi, chrom, start, end, cseq=None, cpos=None, cache_path=None, tail_pad=65536):
    """Return the raw text lines within the region (list[str])

    * The column positions must agree with col_seq / col_beg in the .tbi and **must not be hard-coded**:
        FinnGen summary_stats (including lab_values): col_seq=1, col_beg=2  -> f[1]=chrom, f[2]=pos
        FinnGen FINEMAP.snp                    : col_seq=6, col_beg=7  -> f[6]=chr1..., f[7]=pos
        FinnGen LD matrix                       : col_seq=1, col_beg=2  -> f[1]=pos (note f[2] is the string variant)
      Using the wrong columns makes the filter **silently drop every row** (hit this in the current round: finemap returned "no data in region" everywhere).
    """
    # ** Column positions are always derived automatically from the tbi and **must not be hard-coded** (the biggest pitfall this round):
    #    the col_seq / col_beg that htslib writes into .tbi are **1-based**; subtract 1 when using them:
    #      summary_stats  col_seq=1, col_beg=2  -> cseq=0, cpos=1
    #      FINEMAP.snp    col_seq=6, col_beg=7  -> cseq=5(chromosome), cpos=6(position)
    #    treating them as 0-based makes the filter **silently drop every row** (manifesting as the false negative "no data in region").
    if cseq is None:
        cseq = int(tbi.get("col_seq", 1)) - 1
    if cpos is None:
        cpos = int(tbi.get("col_beg", 2)) - 1
    ref = tbi["refs"].get(chrom)
    if ref is None:
        return []
    chunks = []
    for b in reg2bins(start - 1, end):
        if b in ref["bins"]:
            chunks.extend(ref["bins"][b])
    if not chunks:
        return []
    chunks = _merge_chunks(chunks)
    lo = min(c[0] for c in chunks)
    hi = max(c[1] for c in chunks)
    fs = (lo >> 16)
    # * End padding: the default +64 KB suffices for FinnGen summary_stats;
    #   but `*.FINEMAP.snp.bgz` in practice **under-covers** (the content fetched per tbi chunk stops before the target region,
    #   causing the false negative "0 rows in region") -> such files need a larger padding (use the tail_pad argument).
    fe = (hi >> 16) + tail_pad
    raw = b""
    CH = 4 << 20
    p = fs
    while p < fe:
        q = min(p + CH - 1, fe)
        raw += http_range(data_url, p, q)
        p = q + 1
    txt = bgzf_decompress(raw).decode("utf-8", "replace")
    out = []
    for ln in txt.splitlines():
        f = ln.rstrip("\r").split("\t")
        if len(f) <= max(cseq, cpos):
            continue
        if f[cseq] != chrom:
            continue
        try:
            pos = int(f[cpos])
        except ValueError:
            continue
        if start <= pos <= end:
            out.append(ln.rstrip("\r"))
    return out


def fetch_all(data_url, timeout=600):
    """Whole-package download + BGZF decompression (for **small files**; the safe route when the tbi chunk coverage is incomplete).

    * Tested (2026-10-09): fetching FinnGen `*.FINEMAP.snp.bgz` (1-3 MB) via the tbi **under-covers**
      -- the chunk retrieved per bin ends before the target region, causing the false negative "no data in region";
      while the same tbi returns 0 rows for a large-range query, indicating the index and the file are not fully consistent.
      For such small files, a direct whole-package download is the most reliable (a few MB, about 1 minute).
    """
    req = urllib.request.Request(data_url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        raw = r.read()
    return bgzf_decompress(raw).decode("utf-8", "replace")
