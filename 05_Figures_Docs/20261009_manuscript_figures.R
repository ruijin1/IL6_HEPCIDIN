# =============================================================================
# Manuscript figures 1-6  (IL-6 / hepcidin / bilirubin triangulation)
#   usage: Rscript 20261009_manuscript_figures.R <figdata_dir> <out_dir> [ascii_only]
#   All labels are read from the *_en columns -> this source is 100% ASCII
#   (R on this machine cannot translate \uXXXX escapes; Chinese paths break ggsave).
#   * fig2A: model "(4) non-CoNS" was a duplicate of "+ interval and baseline" in the
#     original script (the exposure was not actually swapped) -> removed in prep.
#   * fig6 systemic arm: "sTfR down" printed the sign the wrong way round for the active
#     (HAMP-up) arm. Attenuated signalling lowers sTfR (fig4a beta = -0.0263; fig6 supporting
#     line "lower sTfR/RDW"), so the active arm must read "sTfR up", as iron restriction
#     raises erythroid transferrin-receptor expression. Also fixes the self-contradiction
#     with the supporting-evidence line on the same panel.
# =============================================================================
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(scales)})

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]; OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
## Always read *_ascii.csv: R on this machine cannot translate non-ASCII bytes into native encoding,
## reading a CSV that contains Chinese directly throws "invalid input found on input connection" and returns 0 rows.
G <- function(f) read.csv(file.path(DIR, sub("\\.csv$", "_ascii.csv", f)),
                           stringsAsFactors = FALSE, fileEncoding = "ASCII", check.names = FALSE)

## SR: "The thinnest lines in the final figure should be no smaller than one point wide."
## On this device a ggplot linewidth of L mm is emitted as L * 2.134 pt, so LW = 0.48 mm
## (1.02 pt) is applied as a floor to every visible stroke: grid lines, panel border, axis
## ticks, reference lines, error bars, tile borders and point-marker outlines.
LW <- 0.48

th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        line = element_line(linewidth = LW),
        panel.grid.major = element_line(linewidth = LW),
        panel.border = element_rect(colour = "grey20", fill = NA, linewidth = LW),
        axis.ticks = element_line(linewidth = LW),
        axis.line = element_line(linewidth = LW),
        legend.ticks = element_line(colour = "white", linewidth = LW),
        strip.background = element_rect(fill = "grey94", colour = NA, linewidth = LW),
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = unit(0.35, "cm"),
        plot.tag = element_text(face = "bold", size = 11),
        plot.title = element_text(size = 9, face = "bold"))

## Figure 1 and Figure 6 are schematic diagrams drawn entirely in data coordinates.
## `th` is built on theme_bw(), which re-enables the panel border, major grid, axis ticks and
## axis text that theme_void() switches off -- appending `th` after theme_void() therefore
## framed the cartoon in a dark box and printed meaningless 0-100 / 0-60 tick labels on it.
## th_void turns that furniture off again while keeping the white panel background.
th_void <- th + theme(axis.title = element_blank(), axis.text = element_blank(),
                      axis.ticks = element_blank(), axis.line = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())

save_fig <- function(p, name, w = 17, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "cm", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "cm", dpi = 300)
  tiff(file.path(OUT, paste0(name, ".tif")), width = w, height = h, units = "cm", res = 600,
       compression = "lzw", type = "cairo"); print(p); dev.off()
  cat(sprintf("  saved %-30s %.0f x %.0f cm\n", name, w, h))
}

## ------------------------------------------------------------------ Figure 1
fig1 <- function() {
  box <- function(x, y, w, h, lab, sub, fill) list(
    annotate("rect", xmin = x, xmax = x + w, ymin = y, ymax = y + h, fill = fill, colour = "grey35", linewidth = LW),
    annotate("text", x = x + w/2, y = y + h*0.70, label = lab, size = 2.9, fontface = "bold"),
    annotate("text", x = x + w/2, y = y + h*0.30, label = sub, size = 2.30, colour = "grey20"))
  ar <- function(x, y0, y1) annotate("segment", x = x, y = y0, xend = x, yend = y1,
                                     arrow = arrow(length = unit(0.12, "cm")), linewidth = LW, colour = "grey45")
  p <- ggplot() + xlim(0, 100) + ylim(0, 64) +
    annotate("text", x = 50, y = 62, size = 3.4, fontface = "bold", label = "Triangulation of three evidence layers") +
    box(1.5, 46, 30, 12, "Clinical", "MIMIC-IV: 17,244 NICU/SCN admissions\nMIMIC-III: 7940 NICU admissions\n8255 blood-culture episodes", "#D6E4F0") +
    box(35, 46, 30, 12, "Genetic", "FinnGen R12: 5 neonatal endpoints\n15 quantitative lab traits (n 5055-447,271)\ncis-pQTL: IL6R, SAA1/SAA2, CRP", "#E6DED1") +
    box(68.5, 46, 30, 12, "Transcriptomic", "6 human hepatocyte systems (GEO)\nIL-6 / LPS / BMAL1 knockdown\n93 genes across five GEO datasets", "#DCEBDC") +
    ar(16.5, 45.7, 40) + ar(50, 45.7, 40) + ar(83.5, 45.7, 40) +
    box(1.5, 28, 30, 11.5, "1. Clinical question", "Does culture-proven infection\nlower bilirubin, and lower MCV /\nMCH (iron-restricted erythropoiesis)?", "#F4F8FC") +
    box(35, 28, 30, 11.5, "2. Genetic instruments", "Does attenuated IL-6 signalling\nraise MCV / MCH / TSAT / bilirubin\nand lower sTfR and RDW?", "#FAF7F2") +
    box(68.5, 28, 30, 11.5, "3. Hepatic arm", "Does IL-6 suppress hepatic UGT1A1\nin vitro (which would RAISE\nbilirubin)?", "#F6FAF6") +
    annotate("text", x = 50, y = 22.5, size = 3.05, fontface = "bold", label = "Two opposing arms of one pathway") +
    annotate("text", x = 50, y = 18.4, size = 2.45, label = "IL-6 up  ->  hepatic UGT1A1 down   (would raise bilirubin)") +
    annotate("text", x = 50, y = 14.9, size = 2.45, label = "IL-6 up  ->  hepcidin (HAMP) up  ->  iron restriction  ->  less erythropoietic haem turnover   (would lower bilirubin)") +
    annotate("text", x = 50, y = 10.1, size = 3.0, fontface = "bold", colour = "#8C1D18",
             label = "Net effect in vivo: the erythropoietic arm predominates  ->  bilirubin falls") +
    annotate("text", x = 50, y = 5.8, size = 2.3, colour = "grey25",
             label = "Clinical, genetic and HAMP-locus evidence are concordant; the hepatocyte data reveal the opposing hepatic arm.") +
    annotate("text", x = 50, y = 2.8, size = 2.3, colour = "grey40",
             label = "Colocalisation at IL6R: the proximal signal is shared with MCH, RDW, MCV, MCHC and platelet count (PP.H4 up to 0.99)") +
    theme_void() + th_void
  save_fig(p, "Figure1_study_design", 17, 16)
}

## ------------------------------------------------------------------ Figure 2
fig2 <- function() {
  a <- G("fig2A_bilirubin_mimic3.csv")
  lev <- c("CoNS only vs negative", "Unadjusted", "No phototherapy in window",
           "+ phototherapy in window", "+ interval and baseline value", "Definite pathogen vs negative")
  a$lab <- factor(a$label_en, levels = lev[lev %in% a$label_en])
  pA <- ggplot(a, aes(x = b, y = lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = LW) +
    geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = LW, size = 0.28, colour = "#2166AC") +
    geom_text(aes(x = -Inf, label = sprintf("n=%d (%d pos.)", n, ncase)), hjust = -0.04, size = 2.05, colour = "grey30") +
    labs(x = "Delta bilirubin, mg/dL (positive minus negative culture, 72 h)", y = NULL, tag = "a",
         title = "MIMIC-III NICU: within-subject change in bilirubin around each blood culture") +
    coord_cartesian(xlim = c(-4.50, 0.45)) + th
  b <- G("fig2B_bilirubin_mimic4.csv")
  b$group <- factor(b$group, levels = c("No infection", "Culture-positive"))
  pB <- ggplot(b, aes(x = window_mid, y = median, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.16, colour = NA) +
    geom_line(linewidth = LW) + geom_point(size = 1.1, stroke = 0.85) +
    geom_text(aes(label = paste0("n=", n)), vjust = -1.0, size = 2.0, show.legend = FALSE) +
    scale_colour_manual(values = c("No infection" = "#4D4D4D", "Culture-positive" = "#B2182B")) +
    scale_fill_manual(values = c("No infection" = "#4D4D4D", "Culture-positive" = "#B2182B")) +
    scale_x_continuous(breaks = seq(0, 48, 12)) +
    labs(x = "Age at first bilirubin measurement (hours)", y = "Total bilirubin (mg/dL)",
         colour = NULL, fill = NULL, tag = "b",
         title = "MIMIC-IV: median (IQR) first bilirubin by age window (overlapping window only)") +
    th + theme(legend.position = c(0.70, 0.85), legend.background = element_rect(fill = "white", colour = NA))
  save_fig(pA / pB + plot_layout(heights = c(1.05, 1)), "Figure2_bilirubin", 17, 13.5)
}

## ------------------------------------------------------------------ Figure 3
fig3 <- function() {
  x <- G("fig3A_iron_forest.csv")
  x <- x[!is.na(x$model_en) & !is.na(x$exposure_en), ]
  x$outcome_base <- factor(x$outcome_base, levels = c("MCV", "MCH", "RDW", "HCT"))
  x$lab <- paste0(x$model_en, "   [", x$exposure_en, "]")
  ord <- c("Pre-transfusion only   [Definite pathogen]", "Never transfused only   [Definite pathogen]",
           "+ adjusted for transfusion   [Definite pathogen]", "Main analysis   [Definite pathogen]",
           "Pre-transfusion only   [Any positive culture]", "Never transfused only   [Any positive culture]",
           "+ adjusted for transfusion   [Any positive culture]", "Main analysis   [Any positive culture]")
  pA <- ggplot(x, aes(x = b, y = factor(lab, levels = ord))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = LW) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = exposure_en), linewidth = LW, size = 0.3) +
    scale_colour_manual(values = c("Definite pathogen" = "#B2182B", "Any positive culture" = "#4393C3"), name = NULL) +
    facet_wrap(~outcome_base, scales = "free_x", nrow = 2) +
    labs(x = "Adjusted difference (original units)", y = NULL, tag = "a",
         title = "MIMIC-III: erythrocyte indices by culture status,\nwith the three prespecified analyses for transfusion") +
    th + theme(legend.position = "top")
  d <- G("fig3B_iron_dose.csv")
  d$outcome <- factor(toupper(d$outcome), levels = c("MCV", "MCH", "RDW", "HCT"))
  d$lab <- ifelse(d$n_pos == 0, "0", ifelse(d$n_pos == 1, "1", "2+"))
  pB <- ggplot(d, aes(x = factor(lab, levels = c("0", "1", "2+")), y = mean, group = outcome)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = LW, colour = "#B2182B") +
    geom_line(linewidth = LW, colour = "#B2182B") + geom_point(size = 1.2, colour = "#B2182B", stroke = 0.85) +
    geom_text(aes(label = paste0("n=", n)), vjust = -1.6, size = 2.0, colour = "grey30") +
    facet_wrap(~outcome, scales = "free_y", nrow = 1) +
    labs(x = "Number of positive blood cultures", y = "Mean value (95% CI)", tag = "b",
         title = "Monotonic dose-response across the number of positive cultures") + th
  save_fig(pA / pB + plot_layout(heights = c(1.35, 1)), "Figure3_iron_erythropoiesis", 17, 19)
}

## ------------------------------------------------------------------ Figure 4
fig4 <- function() {
  a <- G("fig4A_genetic_forest.csv")
  a$het <- ifelse(a$Q_p < 0.05, "Heterogeneous (not interpreted)", "Homogeneous")
  lv <- c("MCH", "MCV", "RDW", "sTfR", "Transferrin saturation", "Total bilirubin", "Ferritin",
          "Haemoglobin", "Serum iron", "Haematocrit", "Red-cell count", "Reticulocyte count",
          "Reticulocytes (%)", "Transferrin", "EPO")
  a <- a[a$trait_en %in% lv, ]; a$trait_en <- factor(a$trait_en, levels = rev(lv))
  pA <- ggplot(a, aes(x = beta_ivw, y = trait_en)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = LW) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = het), linewidth = LW, size = 0.3) +
    scale_colour_manual(values = c("Homogeneous" = "#B2182B", "Heterogeneous (not interpreted)" = "#9E9E9E"), name = NULL) +
    labs(x = "SD change per SD higher genetically predicted soluble IL-6 receptor", y = NULL, tag = "a",
         title = "IL6R cis-pQTL: correlation-aware IVW across iron and erythropoiesis traits") +
    th + theme(legend.position = "top")
  h <- G("fig4B_hamp_locus.csv")
  h$trait_en <- factor(h$trait_en, levels = rev(h$trait_en[order(h$neglog10p)]))
  pB <- ggplot(h, aes(x = neglog10p, y = trait_en)) +
    geom_col(fill = "#67A9CF", width = 0.68) +
    geom_vline(xintercept = -log10(5e-8), linetype = 2, colour = "#B2182B", linewidth = LW) +
    geom_text(aes(label = sprintf("%.1e", p)), hjust = -0.18, size = 2.0, colour = "grey25") +
    labs(x = "minus log10 minimum p in the HAMP region (chr19:35.1-35.4 Mb, GRCh38)", y = NULL, tag = "b",
         title = "The HAMP locus itself, independent of the IL6R instruments") +
    coord_cartesian(xlim = c(0, 20)) + th
  save_fig(pA / pB + plot_layout(heights = c(1.3, 1)), "Figure4_genetic_evidence", 17, 17)
}

## ------------------------------------------------------------------ Figure 5
fig5 <- function() {
  d <- G("fig5_hepatocyte.csv")
  keep <- c("UGT1A1", "UGT1A6", "UGT1A9", "HAMP", "CYP3A4", "CYP1A2", "ALB", "NR1D1", "CRP", "SOCS3")
  d <- d[d$gene %in% keep, ]
  d$gene <- factor(d$gene, levels = rev(keep))
  ord <- c("GSE202045\nPrimary hepatocytes\nIL-6", "GSE202045\nHepG2\nIL-6",
           "GSE230325\ndHepaRG\nLPS", "GSE218243\nPrimary hepatocytes\nBMAL1 knockdown",
           "GSE275342\nHep3B\nIL-6", "GSE275342\nHuh7\nIL-6")
  d$dataset_en <- factor(d$dataset_en, levels = ord[ord %in% d$dataset_en])
  p <- ggplot(d, aes(x = dataset_en, y = gene)) +
    geom_tile(aes(fill = log2FC), colour = "white", linewidth = LW) +
    geom_text(aes(label = sprintf("%+.2f", log2FC), colour = abs(log2FC) > 3.2), size = 2.0) +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey15"), guide = "none") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-3.2, 3.2), oob = squish, name = "log2FC") +
    labs(x = NULL, y = NULL,
         title = "Direction of UGT1A1, HAMP and reference genes across independent human hepatocyte systems") +
    th + theme(axis.text.x = element_text(size = 6.6))
  save_fig(p, "Figure5_hepatocyte_systems", 17, 10)
}

## ------------------------------------------------------------------ Figure 6
fig6 <- function() {
  p <- ggplot() + xlim(0, 100) + ylim(0, 60) +
    annotate("text", x = 50, y = 58, size = 3.4, fontface = "bold", label = "Proposed two-arm model of IL-6 effects on bilirubin") +
    annotate("rect", xmin = 35, xmax = 65, ymin = 49, ymax = 53.5, fill = "#F0F0F0", colour = "grey40", linewidth = LW) +
    annotate("text", x = 50, y = 51.2, size = 2.9, fontface = "bold", label = "IL-6 signalling") +
    annotate("segment", x = 35, y = 51.2, xend = 20, yend = 42, arrow = arrow(length = unit(0.13, "cm")), linewidth = LW, colour = "#2166AC") +
    annotate("rect", xmin = 2, xmax = 31, ymin = 33.5, ymax = 42, fill = "#E8F0F8", colour = "#2166AC", linewidth = LW) +
    annotate("text", x = 16.5, y = 39.9, size = 2.6, fontface = "bold", colour = "#10467F", label = "Hepatic arm (in vitro)") +
    annotate("text", x = 16.5, y = 36.3, size = 2.3, colour = "grey20", label = "UGT1A1 down (-1.54 to -0.09 in\n5 of 6 systems; absent in Hep3B)") +
    annotate("segment", x = 16.5, y = 33.2, xend = 16.5, yend = 27.5, arrow = arrow(length = unit(0.13, "cm")), linewidth = LW, colour = "#2166AC") +
    annotate("rect", xmin = 2, xmax = 31, ymin = 19.5, ymax = 27.5, fill = "#F7FAFD", colour = "#2166AC", linewidth = LW) +
    annotate("text", x = 16.5, y = 25.4, size = 2.6, fontface = "bold", colour = "#10467F", label = "Would predict bilirubin UP") +
    annotate("text", x = 16.5, y = 22.0, size = 2.2, colour = "grey30", label = "Not observed in vivo") +
    annotate("segment", x = 65, y = 51.2, xend = 80, yend = 42, arrow = arrow(length = unit(0.13, "cm")), linewidth = LW, colour = "#B2182B") +
    annotate("rect", xmin = 69, xmax = 98, ymin = 33.5, ymax = 42, fill = "#FBE9EA", colour = "#B2182B", linewidth = LW) +
    annotate("text", x = 83.5, y = 39.9, size = 2.6, fontface = "bold", colour = "#8C1D18", label = "Systemic arm (in vivo)") +
    annotate("text", x = 83.5, y = 36.3, size = 2.3, colour = "grey20", label = "HAMP up -> iron restriction ->\nMCV/MCH down, RDW up, sTfR up") +
    annotate("segment", x = 83.5, y = 33.2, xend = 83.5, yend = 27.5, arrow = arrow(length = unit(0.13, "cm")), linewidth = LW, colour = "#B2182B") +
    annotate("rect", xmin = 69, xmax = 98, ymin = 19.5, ymax = 27.5, fill = "#FDF5F5", colour = "#B2182B", linewidth = LW) +
    annotate("text", x = 83.5, y = 25.4, size = 2.6, fontface = "bold", colour = "#8C1D18", label = "Observed: bilirubin DOWN") +
    annotate("text", x = 83.5, y = 22.0, size = 2.2, colour = "grey30", label = "MIMIC-IV -0.355 log units;\nMIMIC-III -1.07 mg/dL") +
    annotate("text", x = 50, y = 14.2, size = 2.5, fontface = "bold", label = "Supporting evidence") +
    annotate("text", x = 50, y = 10.6, size = 2.25, colour = "grey20", label = "Clinical: culture-proven infection -> lower bilirubin (MIMIC-IV and MIMIC-III) and lower MCV/MCH, robust to transfusion adjustment") +
    annotate("text", x = 50, y = 7.6, size = 2.25, colour = "grey20", label = "Genetic: IL6R cis-pQTL -> higher MCV/MCH/TSAT/bilirubin, lower sTfR/RDW; HAMP locus independently associated (MCH p = 1.7e-17)") +
    annotate("text", x = 50, y = 4.6, size = 2.25, colour = "grey20", label = "Transcriptomic: IL-6 suppresses UGT1A1 in primary human hepatocytes and dHepaRG (the hepatic arm, opposite direction)") +
    theme_void() + th_void
  save_fig(p, "Figure6_two_arm_model", 17, 16)
}

cat("Figures -> ", OUT, "\n", sep = "")
fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
cat("done\n")
