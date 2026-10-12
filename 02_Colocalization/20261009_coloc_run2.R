# =============================================================================
# Colocalization round 2:  IL6R pQTL x MCH / MCV / bilirubin(FinnGen lab) / neonatal jaundice
#   usage:  Rscript 20261009_coloc_run2.R <input_dir> <output_csv>
#
# Why this round exists
#   Round 1 (20261008_coloc_il6r_bilirubin.R) showed IL6R x total bilirubin PP.H4 = 0.98-0.99,
#   but the later self-audit (report section 15.5/15.9) found that
#     (a) the 4 IL6R instruments are NOT independent (r2 0.00-0.437), and
#     (b) the FinnGen FINEMAP credible set for MCH is centred on chr1:154,409,573 and the
#         MCV credible set lies 86-135 kb away from the pQTL lead variant.
#   => colocalisation of IL6R with MCH/MCV had NOT been demonstrated. This script closes that gap
#      with a formal coloc.abf, using
#        exposure = IL6R cis-pQTL (OpenGWAS, GRCh37, round-1 cached files)
#        outcome  = MCH / MCV / bilirubin from FinnGen R12 lab_values (GRCh38, tabix),
#                   harmonised by a constant liftover offset (+27,524 bp) and allele matching
#        plus a NEW pair: IL6R pQTL x neonatal jaundice (FinnGen P16, case-control, rsID-matched).
#
# Data-handling notes carried over from round 1 (all learned the hard way):
#   1. CSV may carry a UTF-8 BOM -> fileEncoding="UTF-8-BOM" + strip BOM, else merge(by="rsid") breaks.
#   2. Never subset with a possibly-NA logical index (m[cond,] inserts all-NA rows and coloc then
#      errors "beta contains missing values"); use which(); donate MAF from another GWAS at the
#      same rsID when eaf is missing.
#   3. Deduplicate rsIDs (keep the most significant).
#   4. Report the number of shared variants: <200 means underpowered and a 50/50 H3/H4 split is
#      uninformative, NOT evidence of "no colocalization".
#   5. Case-control outcomes need type="cc" AND s (case proportion).
# =============================================================================
suppressMessages(library(coloc))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript 20261009_coloc_run2.R <input_dir> <output_csv>")
DIR <- args[1]; OUT <- args[2]

read_reg <- function(f) {
  d <- read.csv(file.path(DIR, paste0(f, ".csv")), stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
  names(d)[1] <- sub("^\xef\xbb\xbf", "", names(d)[1])
  for (v in c("position", "eaf", "beta", "se", "p", "n")) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d <- d[!is.na(d$beta) & !is.na(d$se) & d$se > 0 & !is.na(d$rsid) & nzchar(d$rsid), ]
  d <- d[order(d$p), ]
  d <- d[!duplicated(d$rsid), ]
  d
}

# ---- MAF donor map (variant frequency is a property of the variant)
ALL <- sub("\\.csv$", "", list.files(DIR, pattern = "\\.csv$"))
ALL <- ALL[grepl("^(exp_|out_)", ALL)]
donor <- list()
for (f in ALL) {
  d <- tryCatch(read_reg(f), error = function(e) NULL)
  if (is.null(d) || !("eaf" %in% names(d))) next
  d <- d[!is.na(d$eaf), , drop = FALSE]
  for (i in seq_len(nrow(d))) {
    r <- d$rsid[i]
    if (is.null(donor[[r]])) donor[[r]] <- pmin(d$eaf[i], 1 - d$eaf[i])
  }
}
cat("MAF donor map:", length(donor), "variants from", length(ALL), "files\n\n")

align_pair <- function(expfile, outfile, maf_cut = 0.01) {
  e <- read_reg(expfile); o <- read_reg(outfile)
  m <- merge(e, o, by = "rsid", suffixes = c(".e", ".o"))
  if (nrow(m) == 0) return(list(m = m, note = "no overlapping rsID"))
  keep <- which((m$ea.e == m$ea.o & m$nea.e == m$nea.o) |
                  (m$ea.e == m$nea.o & m$nea.e == m$ea.o))
  m <- m[keep, ]
  if (nrow(m) == 0) return(list(m = m, note = "no allele-compatible rsID"))
  f <- which(m$ea.e == m$nea.o & m$nea.e == m$ea.o)
  if (length(f)) m$beta.o[f] <- -m$beta.o[f]
  maf.e <- pmin(m$eaf.e, 1 - m$eaf.e)
  maf.o <- pmin(m$eaf.o, 1 - m$eaf.o)
  na_e <- which(is.na(maf.e)); na_o <- which(is.na(maf.o))
  if (length(na_e)) maf.e[na_e] <- as.numeric(sapply(m$rsid[na_e], function(r) donor[[r]] %||% NA))
  if (length(na_o)) maf.o[na_o] <- as.numeric(sapply(m$rsid[na_o], function(r) donor[[r]] %||% NA))
  if (length(na_e) || length(na_o))
    cat(sprintf("   [MAF donated] exposure %d, outcome %d variants\n", length(na_e), length(na_o)))
  keep2 <- which(!is.na(maf.e) & !is.na(maf.o) & maf.e > maf_cut & maf.o > maf_cut)
  m <- m[keep2, ]; maf.e <- maf.e[keep2]; maf.o <- maf.o[keep2]
  m$maf.e <- maf.e; m$maf.o <- maf.o
  list(m = m, note = "ok")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

run_coloc <- function(expfile, outfile, tag, Nexp, Nout, p12 = 1e-5, maf_cut = 0.01,
                      sdY_forced = NA, otype = "quant", s_prop = NA,
                      win_center = NA, win_half = NA) {
  ap <- align_pair(expfile, outfile, maf_cut)
  m <- ap$m
  if (!is.na(win_center) && !is.na(win_half)) m <- m[which(abs(m$position.e - win_center) <= win_half), ]
  if (nrow(m) < 20) {
    cat(sprintf("%-50s shared=%d (%s) -> skipped\n", tag, nrow(m), ap$note)); return(NULL)
  }
  # * Do NOT pass position to coloc: coloc.abf internally does merge(df1, df2),
  #   and the merge key is **all columns shared by both sides** = snp + position.
  #   In this project the exposure side is GRCh37 and the FinnGen outcome side is GRCh38 -> position is bound to
  #   differ -> merge yields 0 rows and errors out with
  #   "dataset1 and dataset2 should contain the same snps in the same order..." (highly misleading).
  #   Therefore pass only snp on both sides (position is kept in m solely for the sliding window and reporting).
  d1 <- list(beta = m$beta.e, varbeta = m$se.e^2, snp = m$rsid,
             type = "quant", N = Nexp, MAF = m$maf.e)
  d2 <- list(beta = m$beta.o, varbeta = m$se.o^2, snp = m$rsid,
             type = otype, N = Nout, MAF = m$maf.o)
  if (otype == "cc" && !is.na(s_prop)) d2$s <- s_prop
  if (!is.na(sdY_forced)) { d1$sdY <- sdY_forced; if (otype == "quant") d2$sdY <- sdY_forced }
  r <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p12 = p12)),
                error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(NULL)
  sh <- r$summary
  topp <- function(h) {
    pp <- r$results[[paste0("SNP.PP.H", h)]]
    if (is.null(pp) || all(is.na(pp))) return("-")
    i <- which.max(pp); sprintf("%s(%.2f)", r$results$snp[i], pp[i])
  }
  lead_e <- m$rsid[which.min(m$p.e)]; lead_o <- m$rsid[which.min(m$p.o)]
  cat(sprintf("%-50s nSNP=%-5d H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.3f | leadE %s p=%.1e / leadO %s p=%.1e\n",
              tag, sh[["nsnps"]], sh[["PP.H0.abf"]], sh[["PP.H1.abf"]], sh[["PP.H2.abf"]],
              sh[["PP.H3.abf"]], sh[["PP.H4.abf"]], lead_e, min(m$p.e), lead_o, min(m$p.o)))
  data.frame(tag = tag, outcome_type = otype, p12 = p12, maf_cut = maf_cut, sdY_forced = sdY_forced,
             nsnp = sh[["nsnps"]], PP.H0 = sh[["PP.H0.abf"]], PP.H1 = sh[["PP.H1.abf"]],
             PP.H2 = sh[["PP.H2.abf"]], PP.H3 = sh[["PP.H3.abf"]], PP.H4 = sh[["PP.H4.abf"]],
             lead_exposure = lead_e, lead_outcome = lead_o,
             minP_exposure = min(m$p.e), minP_outcome = min(m$p.o),
             top_H4_snp = topp(4), stringsAsFactors = FALSE)
}

cat(strrep("=", 112), "\n", sep = "")
cat("Colocalization round 2: IL6R pQTL x MCH / MCV / bilirubin / neonatal jaundice   ",
    format(Sys.time()), "\n", sep = "")
cat(strrep("=", 112), "\n\n", sep = "")

N_MCH <- 447004; N_MCV <- 447040; N_BIL <- 190407
N_JAUND <- 500187; S_JAUND <- 499 / 500187   # case proportion

res <- list(); add <- function(x) if (!is.null(x)) res[[length(res) + 1]] <<- x

cat("-- A. positive control: IL6R x CRP (shared causal variant established)\n")
add(run_coloc("exp_IL6R_folkersen2018", "out_crp_said2022", "IL6R(F2018) x CRP [pos ctrl]", 3394, 575531))

cat("\n-- B. NEW: IL6R x MCH  (FinnGen R12 lab_values, n=447,004)  <-- the gap from section 15.9\n")
add(run_coloc("exp_IL6R_folkersen2018", "out_mch_finngen", "IL6R(F2018) x MCH FinnGen", 3394, N_MCH))
add(run_coloc("exp_IL6R_folkersen2020", "out_mch_finngen", "IL6R(F2020) x MCH FinnGen", 21758, N_MCH))
add(run_coloc("exp_IL6R_pietzner2020",  "out_mch_finngen", "IL6R(P2020) x MCH FinnGen", 10708, N_MCH))

cat("\n-- C. NEW: IL6R x MCV  (FinnGen R12 lab_values, n=447,040)\n")
add(run_coloc("exp_IL6R_folkersen2018", "out_mcv_finngen", "IL6R(F2018) x MCV FinnGen", 3394, N_MCV))
add(run_coloc("exp_IL6R_folkersen2020", "out_mcv_finngen", "IL6R(F2020) x MCV FinnGen", 21758, N_MCV))
add(run_coloc("exp_IL6R_pietzner2020",  "out_mcv_finngen", "IL6R(P2020) x MCV FinnGen", 10708, N_MCV))

cat("\n-- D. replication: IL6R x total bilirubin (FinnGen quantitative lab, n=190,407)\n")
add(run_coloc("exp_IL6R_folkersen2018", "out_bil_finngen_lab", "IL6R(F2018) x Bil FinnGen-lab", 3394, N_BIL))
add(run_coloc("exp_IL6R_folkersen2020", "out_bil_finngen_lab", "IL6R(F2020) x Bil FinnGen-lab", 21758, N_BIL))

cat("\n-- E. NEW: IL6R x neonatal jaundice (FinnGen P16, case-control, ncase=499)\n")
add(run_coloc("exp_IL6R_folkersen2018", "out_jaund_finngen", "IL6R(F2018) x Neonatal jaundice", 3394, N_JAUND,
              otype = "cc", s_prop = S_JAUND))
add(run_coloc("exp_IL6R_folkersen2020", "out_jaund_finngen", "IL6R(F2020) x Neonatal jaundice", 21758, N_JAUND,
              otype = "cc", s_prop = S_JAUND))
add(run_coloc("exp_IL6R_pietzner2020",  "out_jaund_finngen", "IL6R(P2020) x Neonatal jaundice", 10708, N_JAUND,
              otype = "cc", s_prop = S_JAUND))

cat("\n-- F. negative control: SAA1 x MCH (MR was null -> expect low H4)\n")
add(run_coloc("exp_SAA1_pietzner2020", "out_mch_finngen", "SAA1(P2020) x MCH FinnGen [neg ctrl]", 10708, N_MCH))

cat("\n-- G. sensitivity: p12 prior on the best-powered pair\n")
for (pp in c(1e-4, 1e-6)) {
  add(run_coloc("exp_IL6R_folkersen2018", "out_mch_finngen",
                sprintf("IL6R(F2018) x MCH FinnGen [p12=%.0e]", pp), 3394, N_MCH, p12 = pp))
}

cat("\n-- H. sensitivity: MAF > 5%\n")
add(run_coloc("exp_IL6R_folkersen2018", "out_mch_finngen", "IL6R(F2018) x MCH FinnGen [MAF>5%]", 3394, N_MCH, maf_cut = 0.05))
add(run_coloc("exp_IL6R_folkersen2018", "out_mcv_finngen", "IL6R(F2018) x MCV FinnGen [MAF>5%]", 3394, N_MCV, maf_cut = 0.05))

cat("\n-- I. sliding window centred on the pQTL lead variant (GRCh38 1:154,453,788 = hg19 154,426,264)\n")
CTR <- 154426264
for (hw in c(100000, 50000, 25000)) {
  add(run_coloc("exp_IL6R_folkersen2018", "out_mch_finngen",
                sprintf("IL6R(F2018) x MCH FinnGen [+-%dkb]", hw / 1000), 3394, N_MCH,
                win_center = CTR, win_half = hw))
}
cat("\n")

if (length(res) > 0) {
  out <- do.call(rbind, res)
  write.csv(out, OUT, row.names = FALSE, fileEncoding = "UTF-8")
  cat(strrep("=", 112), "\n", sep = "")
  cat("saved:", OUT, "\n")
  cat(sprintf("PP.H4 > 0.8 : %d / %d ;  PP.H3 > 0.8 : %d / %d ;  H4/H3 > 3 : %d\n",
              sum(out$PP.H4 > 0.8), nrow(out), sum(out$PP.H3 > 0.8), nrow(out),
              sum(out$PP.H4 / pmax(out$PP.H3, 1e-12) > 3)))
} else cat("!! no results\n")
