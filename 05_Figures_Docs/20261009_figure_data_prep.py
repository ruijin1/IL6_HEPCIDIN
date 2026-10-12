# -*- coding: utf-8 -*-
# NOTE: Local working directories were renamed to ASCII equivalents for public release
#       (e.g. D:\workbuddy_workspace\03-data\MR_finngen). Analysis logic is unchanged.
r"""Prepare data for the manuscript figures (using only verified data sources, inventing no numbers)

Outputs under `03-data/Figures/figdata/`:
  fig2A_bilirubin_mimic3.csv    MIMIC-III pre/post paired delta bilirubin (7 models; from D_bilirubin_phototherapy_robust.csv)
  fig2B_bilirubin_mimic4.csv    MIMIC-IV "postnatal-age aligned" comparison: median first bilirubin by postnatal-age bin (positive vs no infection)
  fig3A_iron_forest.csv         MIMIC-III erythrocyte indices x 4 transfusion-robustness models (from C_transfusion_robust_results.csv)
  fig3B_iron_dose.csv           MIMIC-III dose-response (group means for 0/1/>=2 positive culture episodes) -- recomputed from the database and saved
  fig4A_genetic_forest.csv      correlation-aware IVW for the relevant instruments (from MR_il6r_iron_IVW_correlated.csv)
  fig4B_hamp_locus.csv          minimum p in the HAMP locus region (from finngen_lab_locus_scan.csv)
  fig5_hepatocyte.csv           direction of UGT1A1 / HAMP across four hepatocyte systems (from mechanism_summary_panel.csv)

Each panel is written twice:
  <name>.csv          full table, UTF-8 with BOM (for inspection / the record)
  <name>_ascii.csv    the copy the R figure script reads: only the columns that script
                      uses, label columns renamed to *_en, pure ASCII without a BOM
                      (R on this machine cannot translate non-ASCII bytes)
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd

MR = r"D:\workbuddy_workspace\03-data\MR_finngen"
IRON = r"D:\workbuddy_workspace\03-data\MIMIC_NICU_iron"
OUT = r"D:\workbuddy_workspace\03-data\Figures\figdata"
os.makedirs(OUT, exist_ok=True)

# The saved MIMIC-IV exposure table predates the English release and stores the group
# labels in the original form. Map them to the canonical English labels so the filters
# below work against either the old file or a freshly regenerated one.
INFECT_GROUP_ALIAS = {
    "1 血培养阳性": "1 Blood culture positive",
    "4 无感染暴露": "4 No infection exposure",
}

# ---------------------------------------------------------------------------
# Labels used by the R figure script.  The result tables carry the original
# labels; these tables map them to the short English forms the published
# figures use.  The tables (and the R script's axis levels) are the authority
# for the figure wording, which differs from the manuscript table wording.
# ---------------------------------------------------------------------------
FIG2A_LABEL_EN = {
    "主模型（阳性 vs 阴性）": "Unadjusted",
    "+ 调整前后时间间隔与基线值": "+ interval and baseline value",
    "① 调整窗口内光疗": "+ phototherapy in window",
    "② 仅窗口内无光疗": "No phototherapy in window",
    "⑤ 确定病原(非CoNS)阳性 vs 阴性": "Definite pathogen vs negative",
    "⑥ 仅CoNS阳性 vs 阴性": "CoNS only vs negative",
}

FIG3A_MODEL_EN = {
    "主（全部）": "Main analysis",
    "① 调整输血": "+ adjusted for transfusion",
    "② 仅未输血": "Never transfused only",
    "③ 仅输血前化验": "Pre-transfusion only",
}
FIG3A_EXPOSURE_EN = {
    "pos_blood": "Any positive culture",
    "pos_noncons": "Definite pathogen",
}

# Trait labels used by Figures 4a / 4b
TRAIT_EN = {
    "血红蛋白 Hemoglobin": "Haemoglobin",
    "红细胞计数": "Red-cell count",
    "血细胞比容 Hct": "Haematocrit",
    "网织红细胞%": "Reticulocytes (%)",
    "网织红细胞计数": "Reticulocyte count",
    "促红细胞生成素 EPO": "EPO",
    "铁蛋白 Ferritin": "Ferritin",
    "血清铁 Iron": "Serum iron",
    "转铁蛋白 Transferrin": "Transferrin",
    "转铁蛋白饱和度 TSAT": "Transferrin saturation",
    "可溶性转铁蛋白受体 sTfR": "sTfR",
    "总胆红素（参照）": "Total bilirubin",
}

# Hepatocyte-system labels used by Figure 5 (three lines per tile label)
FIG5_DATASET_EN = {
    "GSE202045 HepG2":          "GSE202045\nHepG2\nIL-6",
    "GSE202045 Primary":        "GSE202045\nPrimary hepatocytes\nIL-6",
    "GSE230325 dHepaRG LPS":    "GSE230325\ndHepaRG\nLPS",
    "GSE218243 BMAL1 KD":       "GSE218243\nPrimary hepatocytes\nBMAL1 knockdown",
    "GSE275342 Hep3B":          "GSE275342\nHep3B\nIL-6",
    "GSE275342 Huh7":           "GSE275342\nHuh7\nIL-6",
}


def write_ascii(df, name, cols):
    """Write the plain-ASCII copy that the R figure script reads.

    R on this machine cannot translate non-ASCII bytes into the native encoding, so the
    exported copy keeps only the columns the R script actually uses and is written without
    a BOM in pure ASCII.  A non-ASCII byte would raise here, which is the intended check.
    """
    df[cols].to_csv(os.path.join(OUT, name + "_ascii.csv"), index=False, encoding="ascii")


def fig2A():
    d = pd.read_csv(os.path.join(MR, "coloc", "..", "..", "MIMIC_NICU_iron",
                                 "D_bilirubin_phototherapy_robust.csv"), encoding="utf-8-sig")
    # * 2026-10-10: the source table now also records models that were "insufficient sample, not estimable" (estimable=0, e.g. (3) n=125 positive 5),
    #   whose b/lo/hi/p are NaN. They must be removed before plotting, otherwise fig 2a would gain an empty row that cannot be drawn.
    if "estimable" in d.columns:
        d = d[d["estimable"] == 1]
    d = d.dropna(subset=["b"])
    d = d.rename(columns={"model": "label"})
    d.to_csv(os.path.join(OUT, "fig2A_bilirubin_mimic3.csv"), index=False, encoding="utf-8-sig")
    write_ascii(d.assign(label_en=d.label.replace(FIG2A_LABEL_EN)), "fig2A_bilirubin_mimic3",
                ["label_en", "n", "ncase", "b", "lo", "hi", "p"])
    print("fig2A: %d estimable models (insufficient-sample ones removed)" % len(d)); print(d.to_string(index=False))


def fig2B():
    d = pd.read_csv(os.path.join(MR, "observ_infection_bilirubin_adm.csv"), encoding="utf-8-sig")
    d["infect_group"] = d.infect_group.replace(INFECT_GROUP_ALIAS)
    d = d[d.infect_group.isin(["1 Blood culture positive", "4 No infection exposure"])].copy()
    d["grp"] = np.where(d.infect_group == "1 Blood culture positive", "Culture-positive", "No infection")
    bins = np.arange(0, 169, 12)
    d = d[(d.bili_first_age_h >= 0) & (d.bili_first_age_h <= 168)]
    d["win"] = pd.cut(d.bili_first_age_h, bins, right=False)
    rows = []
    # * Plot only the overlapping window where both groups have data (<=48h): first bilirubin in culture-positive infants is almost always
    #   sampled within 48h, whereas the "no infection" group only has abundant data after 72h (screening at postnatal days 3-5) -- comparing across windows confounds the testing indication.
    for (w, g), sub in d.groupby(["win", "grp"], observed=True):
        if w.left >= 48:
            continue
        if len(sub) < 15:
            continue
        rows.append(dict(window_mid=(w.left + w.right) / 2, group=g, n=len(sub),
                         median=sub.bili_first.median(),
                         q1=sub.bili_first.quantile(.25), q3=sub.bili_first.quantile(.75)))
    R = pd.DataFrame(rows).sort_values(["group", "window_mid"])
    R.to_csv(os.path.join(OUT, "fig2B_bilirubin_mimic4.csv"), index=False, encoding="utf-8-sig")
    write_ascii(R, "fig2B_bilirubin_mimic4", ["window_mid", "group", "n", "median", "q1", "q3"])
    print("\nfig2B: %d rows" % len(R)); print(R.to_string(index=False))


def fig3A():
    d = pd.read_csv(os.path.join(IRON, "C_transfusion_robust_results.csv"), encoding="utf-8-sig")
    d["outcome_base"] = d.outcome.astype(str).str.replace(r"\(preTx\)$", "", regex=True).str.upper()
    d["model_en"] = d.model.replace(FIG3A_MODEL_EN)
    d["exposure_en"] = d.exposure.replace(FIG3A_EXPOSURE_EN)
    d.to_csv(os.path.join(OUT, "fig3A_iron_forest.csv"), index=False, encoding="utf-8-sig")
    write_ascii(d, "fig3A_iron_forest",
                ["outcome_base", "model_en", "exposure_en", "n", "ncase", "b", "lo", "hi", "p"])
    print("\nfig3A: %d rows" % len(d))


def fig3B():
    """Recompute the dose-response from the database (group means for 0/1/>=2 positive culture episodes) and save it to disk."""
    import psycopg2
    c = psycopg2.connect(host="localhost", port=5432, user="postgres", password="1234", dbname="mimiciii")
    base = pd.read_csv(os.path.join(IRON, "nicu_base.csv"), encoding="utf-8-sig",
                       parse_dates=["dob", "intime"])
    lab = pd.read_csv(os.path.join(IRON, "nicu_cbc.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    cul = pd.read_csv(os.path.join(IRON, "nicu_culture.csv"), encoding="utf-8-sig")
    w = pd.read_csv(os.path.join(IRON, "nicu_weight.csv"), encoding="utf-8-sig", parse_dates=["charttime"])
    base = base[base.age_days_at_icu <= 28]
    b = w[w.label.astype(str).str.contains("Birth Weight", case=False, na=False)]
    bw = b[(b.valuenum > 0.2) & (b.valuenum < 7)].sort_values("charttime").groupby("hadm_id").valuenum.first()
    npos = cul[cul.hadm_id.isin(base.hadm_id) & (cul.positive == 1)].groupby("hadm_id").size()
    lab = lab[lab.hadm_id.isin(base.hadm_id)].merge(base[["hadm_id", "dob"]], on="hadm_id")
    lab["dol"] = (lab.charttime - lab.dob).dt.total_seconds() / 86400
    lab = lab[(lab.dol >= 3) & (lab.dol <= 14) & lab["var"].isin(["mcv", "mch", "rdw", "hct"])]
    agg = lab.groupby(["hadm_id", "var"]).valuenum.mean().reset_index()
    wide = agg.pivot(index="hadm_id", columns="var", values="valuenum").join(bw.rename("bw_kg"))
    npos_s = pd.Series(wide.index.map(npos), index=wide.index).fillna(0).astype(int)
    wide["n_pos"] = npos_s.clip(upper=2)
    wide = wide.dropna(subset=["bw_kg"])
    # * The dose-response must use the same sample as the main cross-sectional analysis (A_crosssection_results.csv, n = 710),
    #   otherwise the denominator differs between indices (under the old definition hct reached 957 and the positive count was 120 != 107 in the main analysis).
    core = wide.dropna(subset=["mcv", "mch", "rdw"])
    print("\n[fig3B] unified to the main-analysis sample: n = %d (any positive = %d)"
          % (len(core), int((core.n_pos > 0).sum())))
    rows = []
    for v in ["mcv", "mch", "rdw", "hct"]:
        for k, sub in core.dropna(subset=[v]).groupby("n_pos"):
            if len(sub) < 5:
                continue
            se = sub[v].std(ddof=1) / np.sqrt(len(sub))
            rows.append(dict(outcome=v, n_pos=int(k), n=len(sub), mean=sub[v].mean(),
                             lo=sub[v].mean() - 1.96 * se, hi=sub[v].mean() + 1.96 * se))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "fig3B_iron_dose.csv"), index=False, encoding="utf-8-sig")
    write_ascii(R, "fig3B_iron_dose", ["outcome", "n_pos", "n", "mean", "lo", "hi"])
    print("\nfig3B dose-response:"); print(R.to_string(index=False))


def fig4A():
    d = pd.read_csv(os.path.join(MR, "MR_il6r_iron_IVW_correlated.csv"), encoding="utf-8-sig")
    d["trait_en"] = d.trait.replace(TRAIT_EN)
    d.to_csv(os.path.join(OUT, "fig4A_genetic_forest.csv"), index=False, encoding="utf-8-sig")
    write_ascii(d, "fig4A_genetic_forest",
                ["trait_en", "beta_ivw", "se_ivw", "lo", "hi", "p_ivw", "Q", "Q_p"])
    print("\nfig4A: %d traits" % len(d))


def fig4B():
    d = pd.read_csv(os.path.join(MR, "finngen_lab_locus_scan.csv"), encoding="utf-8-sig")
    h = d[d.region == "HAMP"].groupby("trait").p.min().reset_index()
    h["neglog10p"] = -np.log10(h.p)
    h["trait_en"] = h.trait.replace(TRAIT_EN)
    h = h.sort_values("neglog10p", ascending=False)
    h.to_csv(os.path.join(OUT, "fig4B_hamp_locus.csv"), index=False, encoding="utf-8-sig")
    write_ascii(h, "fig4B_hamp_locus", ["trait_en", "p", "neglog10p"])
    print("\nfig4B minimum p in the HAMP region:"); print(h.to_string(index=False))


def fig5():
    """Direction of the panel genes across the independent human hepatocyte systems.

    One value per gene and system.  For GSE202045 the 24 h contrast is the reported
    timepoint; for GSE230325 the 6 h and 24 h contrasts are averaged (see the code index,
    section on the mechanism layer).  Mouse primary hepatocytes are not part of the figure.
    The gene subset is applied by the R script, so every gene is exported here.
    """
    d = pd.read_csv(os.path.join(MR, "mechanism_summary_panel.csv"), encoding="utf-8-sig")
    d = d[d.found == True].copy()
    d["log2FC"] = pd.to_numeric(d.log2FC, errors="coerce")
    d = d[d.log2FC.notna()]

    cu = d.contrast_unified.astype(str)
    ds = d.dataset.astype(str)
    lab = pd.Series(np.nan, index=d.index, dtype=object)
    lab[cu.str.startswith("GSE202045 HepG2")] = FIG5_DATASET_EN["GSE202045 HepG2"]
    lab[cu.str.startswith("GSE202045 人原代肝细胞")] = FIG5_DATASET_EN["GSE202045 Primary"]
    lab[cu.str.startswith("LPS")] = FIG5_DATASET_EN["GSE230325 dHepaRG LPS"]
    lab[cu.str.startswith("IL-6 vs 载体 (Hep3B)")] = FIG5_DATASET_EN["GSE275342 Hep3B"]
    lab[cu.str.startswith("IL-6 vs 载体 (Huh7)")] = FIG5_DATASET_EN["GSE275342 Huh7"]
    lab[ds.str.startswith("GSE218243")] = FIG5_DATASET_EN["GSE218243 BMAL1 KD"]
    d["dataset_en"] = lab
    d = d[d.dataset_en.notna()]

    # GSE202045 keeps the 24 h contrast only; every other system keeps all its timepoints
    d = d[(~cu.loc[d.index].str.startswith("GSE202045")) | cu.loc[d.index].str.contains("24h")]

    keep = d.groupby(["gene", "dataset_en"], as_index=False).log2FC.mean()
    keep = keep[["gene", "dataset_en", "log2FC"]].sort_values(["gene", "dataset_en"])
    keep.to_csv(os.path.join(OUT, "fig5_hepatocyte.csv"), index=False, encoding="utf-8-sig")
    write_ascii(keep, "fig5_hepatocyte", ["gene", "dataset_en", "log2FC"])
    print("\nfig5: %d rows" % len(keep))
    print(keep.pivot_table(index="gene", columns="dataset_en", values="log2FC").round(2).to_string())


if __name__ == "__main__":
    fig2A(); fig2B(); fig3A(); fig3B(); fig4A(); fig4B(); fig5()
    print("\nAll outputs ->", OUT)
