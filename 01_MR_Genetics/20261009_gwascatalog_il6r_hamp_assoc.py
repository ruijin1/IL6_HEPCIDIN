# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""**Panorama of reported associations** at the IL6R / HAMP loci (new GWAS Catalog API, no authentication required)
   -- Purpose: before using OpenGWAS for formal MR, do a free first pass to screen "whether there are any usable erythroid/iron-metabolism phenotypes"

* New API (the old REST returns 410 Gone):
     https://www.ebi.ac.uk/gwas/api/v2/variants/<rsid>/associations
     https://www.ebi.ac.uk/gwas/api/v2/genes/<sym>/associations
* This endpoint returns **associations already reported in the literature** (not full summary stats) -> use it only for **screening**, never as MR input.
* Phenotype keywords of interest: erythroid (reticulocyte / erythrocyte / haemoglobin / hematocrit / MCV),
   iron metabolism (ferritin / iron / transferrin / hepcidin), inflammation (CRP / IL-6 / sIL6R), bilirubin.
"""
import urllib.request, urllib.parse, json, re, time, sys, os
sys.stdout.reconfigure(encoding="utf-8")

OUT = r"D:\workbuddy_workspace\03-data\MR_finngen"

VARIANTS = ["rs4129267", "rs4601580", "rs7514452", "rs6690230",   # IL6R cis-pQTL
            "rs10421768", "rs855791", "rs1799945", "rs1800562",   # HAMP promoter, TMPRSS6, HFE C282Y/H63D
            "rs4820268", "rs4895441"]                             # TMPRSS6 / HBB loci (erythroid reference)

KEY_PAT = re.compile(r"reticulocyt|erythrocyt|haemoglob|hemoglob|hematocrit|haematocrit|"
                     r"red blood cell|ferritin|\biron\b|transferrin|hepcidin|"
                     r"C-reactive|interleukin|c-reactive|bilirubin|erythropoiet|"
                     r"mean corpuscular|anemia|anaemia|platelet|leukocyte|monocyte", re.I)


def get(u, tries=4, to=70):
    for i in range(tries):
        try:
            r = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
            with urllib.request.urlopen(r, timeout=to) as resp:
                return json.loads(resp.read().decode("utf-8", "replace"))
        except Exception as e:
            if i == tries - 1:
                return {"__err__": str(e)[:100]}
            time.sleep(2.5)


def assoc(rs, size=500):
    """* The v2 API field names differ from the old REST: traitName(list) / efoTraits[{label}] /
       mappedGenes(list of str) / beta(str, e.g. "0.1265 unit increase") / accessionId / pubmedId.
       Using the old field names (efoTraits[].trait, .get on mappedGenes) silently yields a wholly empty table -- hit this in the current round.
    """
    u = "https://www.ebi.ac.uk/gwas/api/v2/variants/%s/associations?size=%d" % (rs, size)
    d = get(u)
    if "__err__" in d:
        return None, d["__err__"]
    emb = d.get("_embedded", {}).get("associations", [])
    total = (d.get("page") or {}).get("totalElements")
    out = []
    for a in emb:
        ra = a.get("riskAllele") or [{}]
        tn = a.get("traitName") or []
        if not tn:
            tn = [t.get("label") for t in (a.get("efoTraits") or []) if t.get("label")]
        p = None
        pe = a.get("pValueExponent")
        if pe is not None:
            mant = a.get("pValue") or 1
            p = float(mant) * (10 ** pe)
        out.append(dict(rsid=rs,
                        risk_allele=(ra[0].get("label") or ""),
                        traits="; ".join(str(x) for x in tn),
                        mapped=",".join(a.get("mappedGenes") or []),
                        beta=str(a.get("beta") or ""),
                        raf=a.get("riskFrequency"),
                        p=p, study=a.get("accessionId") or "", pubmed=a.get("pubmedId") or ""))
    return out, total


def main():
    allrows = []
    for rs in VARIANTS:
        rows, total = assoc(rs)
        print("\n" + "=" * 108)
        print("%s   associations %s" % (rs, ("returned %d entries (of %s)" % (len(rows), total)) if rows is not None else "fetch failed: " + str(total)))
        print("=" * 108)
        if not rows:
            continue
        allrows += rows
        key = [r for r in rows if KEY_PAT.search(r["traits"])]
        key.sort(key=lambda r: (r["p"] if r["p"] else 1))
        print("  keyword hits (erythroid/iron/inflammation/bilirubin) %d / %d entries:" % (len(key), len(rows)))
        for r in key[:22]:
            print("    %-26s %-10s p=%-10s %-8s %s" % (
                r["traits"][:26], r["risk_allele"][:10],
                ("%.2g" % r["p"]) if r["p"] else "NR",
                str(r["beta"])[:10], r["mapped"][:22]))
        if not key:
            print("    (none) all phenotypes: %s" % "; ".join(sorted({r["traits"][:30] for r in rows})[:14]))
        time.sleep(0.5)

    # Summary: which "phenotype categories" have appeared at the IL6R locus
    print("\n\n" + "#" * 108)
    print("# Summary: phenotypes seen at the IL6R locus (rs4129267/rs4601580/rs7514452/rs6690230)")
    print("#" * 108)
    il6r = [r for r in allrows if r["rsid"] in ("rs4129267", "rs4601580", "rs7514452", "rs6690230")]
    tra = {}
    for r in il6r:
        for t in r["traits"].split("; "):
            if t.strip():
                tra[t.strip()] = min(tra.get(t.strip(), 1), r["p"] or 1)
    for t, p in sorted(tra.items(), key=lambda kv: kv[1]):
        star = "*" if KEY_PAT.search(t) else " "
        print("  %s %-70s min p=%.2g" % (star, t[:70], p))

    if allrows:
        import pandas as pd
        p = os.path.join(OUT, "gwas_catalog_IL6R_HAMP_associations.csv")
        pd.DataFrame(allrows).to_csv(p, index=False, encoding="utf-8-sig")
        print("\nFull table saved: %s (%d rows)" % (p, len(allrows)))


if __name__ == "__main__":
    main()
