# -*- coding: utf-8 -*-
r"""Search transcriptome datasets related to "fetal / neonatal hepatocytes" -- to test 13.9#2 (developmental-stage explanation)

Questions to answer (not just fetching data at random):
  Q1 How do fetal/neonatal-derived hepatocytes (or iPS-liver organoids) respond in UGT1A1 under **+/-IL-6 / +/- inflammatory stimulation**?
     -> directly tests the leading explanation that "adult UGT1A1 is non-rate-limiting, neonatal is near rate-limiting"
  Q2 If there is no stimulation arm: in the **baseline** transcriptome of fetal vs adult hepatocytes,
     how different are the expression of IL-6 signaling pathway components (IL6R / IL6ST / JAK1 / STAT3) and the bilirubin pathway (UGT1A1 / HNF4A / ONECUT1)?
     -> judges whether fetal hepatocytes **possess** the molecular prerequisites to mount a UGT1A1 downregulation response to IL-6

* Same discipline as previous rounds: **read `!Series_overall_design` carefully first before deciding whether to download** (trap 7).
"""
import urllib.request, urllib.parse, json, re, time, sys
sys.stdout.reconfigure(encoding="utf-8")


def get(u, tries=4):
    for i in range(tries):
        try:
            r = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
            return urllib.request.urlopen(r, timeout=90).read().decode("utf-8", "replace")
        except Exception as e:
            if i == tries - 1:
                return "ERR:" + str(e)
            time.sleep(2)


def search(term, retmax=40):
    u = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%s&retmax=%d&retmode=json"
         % (urllib.parse.quote(term), retmax))
    t = get(u)
    if t.startswith("ERR:"):
        print("  search failed", t[:80]); return []
    ids = json.loads(t)["esearchresult"]["idlist"]
    if not ids:
        return []
    sj = json.loads(get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=%s&retmode=json"
                        % ",".join(ids)))
    out = []
    for i in ids:
        r = sj["result"].get(i, {})
        out.append(dict(acc=r.get("accession"), n=r.get("n_samples"), taxon=r.get("taxon", ""),
                        title=(r.get("title") or ""), gds=r.get("gdsType", "")))
    return out


QUERIES = {
    "A fetal hepatocyte": '(fetal hepatocyte[All Fields] OR "fetal liver"[Title] AND hepatocyte[All Fields]) AND GSE[ETYP] AND "Homo sapiens"[Organism]',
    "B neonatal hepatocyte": '(neonatal[All Fields] AND hepatocyte[All Fields]) AND GSE[ETYP] AND "Homo sapiens"[Organism]',
    "C iPS/ESC liver organoid + inflammation/IL-6": '("hepatic organoid"[All Fields] OR "liver organoid"[All Fields] OR iHep[All Fields] OR "hepatocyte-like"[All Fields]) AND (IL6[All Fields] OR "IL-6"[All Fields] OR inflammation[All Fields] OR LPS[All Fields] OR "acute phase"[All Fields]) AND GSE[ETYP] AND "Homo sapiens"[Organism]',
    "D fetal vs adult liver transcriptome comparison": '(fetal[Title] AND adult[Title] AND liver[Title]) AND GSE[ETYP] AND "Homo sapiens"[Organism]',
    "E fetal hepatocytes + development/maturity": '("fetal liver"[All Fields] AND (maturation[All Fields] OR ontogeny[All Fields] OR developmental[All Fields])) AND GSE[ETYP] AND "Homo sapiens"[Organism]',
    "F hepatocyte development UGT1A1/HNF4A": '(UGT1A1[All Fields] OR "UDP-glucuronosyltransferase 1"[All Fields] OR HNF4A[All Fields]) AND (fetal[All Fields] OR neonatal[All Fields] OR developmental[All Fields]) AND GSE[ETYP]',
}

ALL = {}
for name, q in QUERIES.items():
    hits = search(q)
    print("\n" + "=" * 110)
    print("%s  -> hits %d" % (name, len(hits)))
    print("=" * 110)
    for h in hits:
        key = h["acc"]
        if key and key not in ALL:
            ALL[key] = (name, h)
        print("  %-12s n=%-5s %-26s %s" % (h["acc"], h["n"], (h["taxon"] or "")[:26], h["title"][:74]))

print("\n\n%d unique datasets in total" % len(ALL))

# ---- for the candidates most likely to have an intervention arm, fetch overall_design + platform to verify (trap 7 + trap 14) ----
CAND = [a for a in ALL]
print("\n\n" + "#" * 110)
print("# Verify !Series_overall_design and platform for each candidate (fetch just these few hundred bytes, no large file downloads)")
print("#" * 110)
for acc in CAND[:40]:
    t = get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&form=text&view=brief" % acc)
    if t.startswith("ERR:"):
        print("\n%s  fetch failed %s" % (acc, t[:60])); continue
    d = {}
    for line in t.splitlines():
        m = re.match(r"^!(Series_title|Series_overall_design|Series_platform_id|Series_type)\s*=\s*(.*)", line)
        if m:
            d[m.group(1)] = m.group(2)
    body = " ".join(d.values()).lower()
    # scoring: whether there is an intervention arm + whether it is fetal/neonatal
    flag = []
    if re.search(r"il-?6|interleukin", body): flag.append("IL6")
    if re.search(r"lps|lipopolysaccharide", body): flag.append("LPS")
    if re.search(r"inflammat|cytokine|acute.?phase", body): flag.append("INFLAM")
    if re.search(r"fetal|foetal|neonat|infant|newborn|pediatric", body): flag.append("PEDI")
    if re.search(r"organoid|ipsc|ips\b|esc\b|stem", body): flag.append("ORGANOID")
    mark = "*" if ("IL6" in flag or "LPS" in flag or "INFLAM" in flag) and "PEDI" in flag else ("·" if flag else " ")
    print("\n%s [%s] %-12s %s" % (mark, ",".join(flag) or "-", acc, d.get("Series_title", "")[:80]))
    print("     platform %s | %s" % (d.get("Series_platform_id", "?"), d.get("Series_type", "")[:60]))
    print("     design %s" % d.get("Series_overall_design", "")[:300])
    time.sleep(0.35)
