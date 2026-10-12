# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
"""Mechanism pillar step 1: search GEO for usable datasets on "inflammation / IL-6 stimulation -> bilirubin metabolism pathway genes"

Search strategy (multiple angles, to avoid relying on a single keyword):
  A. IL-6 / cytokine-stimulated hepatocytes (direct perturbation experiments, strongest mechanistic evidence)
  B. sepsis / neonatal infection whole-blood or liver tissue transcriptomes (phenotype homologous to the clinical cohort)
  C. hyperbilirubinemia / UGT1A1 (positive phenotype)
  D. oxidative stress / antioxidant depletion related to bilirubin
Output: hit count for each query + accession / title / species / sample count / availability of a matrix file for the top datasets
"""
import sys, json, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8")

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(url, to=45):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=to) as r:
        return r.read().decode("utf-8", "replace")


def esearch(db, term, retmax=25):
    u = EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmax": retmax, "retmode": "json"})
    return json.loads(get(u))["esearchresult"]


def esummary(db, ids):
    if not ids:
        return {}
    u = EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
        {"db": db, "id": ",".join(ids), "retmode": "json"})
    return json.loads(get(u)).get("result", {})


TERMS = [
    ("A1", 'interleukin-6 AND hepatocyte AND "Homo sapiens"[Organism] AND gse[Entry Type]'),
    ("A2", 'IL-6 AND (HepG2 OR hepatocyte) AND expression profiling by array[DataSet Type]'),
    ("A3", 'STAT3 AND hepatocyte AND expression profiling by high throughput sequencing[DataSet Type]'),
    ("B1", 'neonatal sepsis AND expression profiling by array[DataSet Type]'),
    ("B2", 'sepsis AND whole blood AND "Homo sapiens"[Organism] AND gse[Entry Type]'),
    ("B3", 'neonatal infection AND transcriptome AND "Homo sapiens"[Organism]'),
    ("C1", 'hyperbilirubinemia AND "Homo sapiens"[Organism] AND gse[Entry Type]'),
    ("C2", 'UGT1A1 AND expression profiling'),
    ("D1", 'oxidative stress AND liver AND "Homo sapiens"[Organism] AND gse[Entry Type]'),
]

seen = {}
for tag, term in TERMS:
    try:
        r = esearch("gds", term, 25)
        n = r.get("count", "?")
        ids = r.get("idlist", [])
        print("=" * 104)
        print("[%s] hits %s   %s" % (tag, n, term[:90]))
        print("=" * 104)
        if not ids:
            continue
        sm = esummary("gds", ids)
        for i in ids:
            d = sm.get(i)
            if not isinstance(d, dict):
                continue
            acc = d.get("accession", "")
            if not acc:
                continue
            tax = d.get("taxon", "")
            ns = d.get("n_samples", "?")
            gtype = d.get("gdstype", "")
            title = (d.get("title") or "")[:96]
            print("   %-12s %-26s n=%-4s %s" % (acc, str(tax)[:26], ns, title))
            if gtype:
                print("        %s" % gtype)
            seen.setdefault(acc, dict(tag=tag, title=d.get("title"), taxon=tax,
                                      n=ns, type=gtype, summary=(d.get("summary") or "")[:400]))
        print()
        time.sleep(0.4)
    except Exception as e:
        print("[%s] search failed %s\n" % (tag, str(e)[:120]))

json.dump(seen, open(r"D:\workbuddy_workspace\03-data\MR_finngen\geo_candidates.json", "w",
                     encoding="utf-8"), ensure_ascii=False, indent=1)
print("collected %d candidate datasets -> geo_candidates.json" % len(seen))
