# Figure 3 only, with the MIMIC-IV panel added.
# usage: Rscript fig3_only.R <in_dir> <out_dir>
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(scales)})

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]; OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
G <- function(f) read.csv(file.path(DIR, sub("\\.csv$", "_ascii.csv", f)),
                          stringsAsFactors = FALSE, fileEncoding = "ASCII", check.names = FALSE)

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

save_fig <- function(p, name, w = 17, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "cm", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "cm", dpi = 300)
  tiff(file.path(OUT, paste0(name, ".tif")), width = w, height = h, units = "cm", res = 600,
       compression = "lzw", type = "cairo"); print(p); dev.off()
  cat(sprintf("  saved %-30s %.0f x %.0f cm\n", name, w, h))
}

fig3 <- function() {
  x <- G("fig3A_iron_forest.csv")
  x <- x[!is.na(x$model_en) & !is.na(x$exposure_en), ]
  x$outcome_base <- factor(x$outcome_base, levels = c("MCV", "MCH", "RDW", "HCT"))
  x$lab <- paste0(x$model_en, "   [", x$exposure_en, "]")
  ord <- c("Pre-transfusion only   [Definite pathogen]", "Never transfused only   [Definite pathogen]",
           "+ adjusted for transfusion   [Definite pathogen]", "Main analysis   [Definite pathogen]",
           "Pre-transfusion only   [Any positive culture]", "Never transfused only   [Any positive culture]",
           "+ adjusted for transfusion   [Any positive culture]", "Main analysis   [Any positive culture]",
           "Main analysis   [CoNS only]")
  pA <- ggplot(x, aes(x = b, y = factor(lab, levels = ord))) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = LW) +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = exposure_en), linewidth = LW, size = 0.3) +
    scale_colour_manual(values = c("Definite pathogen" = "#B2182B", "Any positive culture" = "#4393C3",
                                   "CoNS only" = "#4D4D4D"), name = NULL) +
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

  c3 <- G("fig3C_iron_mimiciv.csv")
  c3$outcome_base <- factor(c3$outcome_base, levels = c("MCV", "MCH", "RDW"))
  c3$exposure_en <- factor(c3$exposure_en, levels = c("CoNS only", "Any positive culture", "Definite pathogen"))
  pC <- ggplot(c3, aes(x = b, y = exposure_en, colour = exposure_en)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = LW) +
    geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = LW, size = 0.3) +
    scale_colour_manual(values = c("Definite pathogen" = "#B2182B", "Any positive culture" = "#4393C3",
                                   "CoNS only" = "#4D4D4D"), guide = "none") +
    facet_wrap(~outcome_base, scales = "free_x", nrow = 1) +
    labs(x = "Adjusted difference (original units)", y = NULL, tag = "c",
         title = "MIMIC-IV: the same pattern, at about half the effect size") +
    th + theme(axis.text.y = element_text(size = 8))

  save_fig(pA / pB / pC + plot_layout(heights = c(1.35, 1.0, 0.8)),
           "Figure3_iron_erythropoiesis", 17, 22)
}

fig3()
