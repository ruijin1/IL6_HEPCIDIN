# =============================================================================
# Colocalization: IL6R / SAA locus x serum bilirubin      (coloc 5.2.3, ASCII source)
# Usage: Rscript coloc_run.R <input_dir> <output_csv>
#
# Design (controls mandatory):
#   main          IL6R pQTL (3 independent studies) x total bilirubin (Barton21 / UKB-Neale / BBJ)
#   positive ctrl IL6R pQTL x CRP        -> expected high PP.H4
#   negative ctrl SAA1 pQTL x bilirubin  -> expected low PP.H4 (MR was null)
#
# Data handling notes (all learned the hard way, keep them):
#  1. CSV may carry a UTF-8 BOM -> read with fileEncoding="UTF-8-BOM" and strip it,
#     otherwise merge(by="rsid") fails: "'by' must specify a uniquely valid column".
#  2. Some GWAS in OpenGWAS have EMPTY eaf (e.g. Said 2022 CRP). Never subset with a
#     logical index that can contain NA (m[cond,] inserts all-NA rows and coloc then
#     complains "beta contains missing values"). Use which(); donate MAF from an
#     independent European GWAS at the same rsID when eaf is missing.
#  3. Deduplicate rsIDs (keep the most significant) - some datasets repeat variants.
#  4. Report the number of shared variants per pair: with <200 SNPs coloc is underpowered
#     and a 50/50 H3/H4 split is uninformative, not evidence of "no colocalization".
#
# Limitation: coloc.abf assumes ONE causal variant per region; coloc.susie cannot be run
# because the OpenGWAS LD endpoints are down (SERVICES__PIPELINE = Error) ->
# declare in the manuscript.
# =============================================================================
suppressMessages(library(coloc))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript coloc_run.R <input_dir> <output_csv>")
DIR <- args[1]; OUT <- args[2]

read_reg <- function(f) {
  d <- read.csv(file.path(DIR, paste0(f, ".csv")), stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
  names(d)[1] <- sub("^\xef\xbb\xbf", "", names(d)[1])
  for (v in c("position", "eaf", "beta", "se", "p", "n")) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d <- d[!is.na(d$beta) & !is.na(d$se) & d$se > 0 & !is.na(d$rsid) & nzchar(d$rsid), ]
  d <- d[order(d$p), ]
  d <- d[!duplicated(d$rsid), ]          # keep the most significant record per rsID
  d
}

# ---- MAF donor: variant frequency is a property of the variant, so borrowing MAF for the
#      same rsID from an independent European GWAS is a reasonable approximation.
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
cat("MAF donor map built from", length(ALL), "files,", length(donor), "variants\n\n")

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
  d_e <- maf.e; d_o <- maf.o
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
                      sdY_forced = NA, win_center = NA, win_half = NA) {
  ap <- align_pair(expfile, outfile, maf_cut)
  m <- ap$m
  if (!is.na(win_center) && !is.na(win_half)) {
    m <- m[which(abs(m$position.e - win_center) <= win_half), ]
  }
  if (nrow(m) < 20) {
    cat(sprintf("%-46s shared=%d (%s) -> skipped\n", tag, nrow(m), ap$note)); return(NULL)
  }
  d1 <- list(beta = m$beta.e, varbeta = m$se.e^2, snp = m$rsid, position = m$position.e,
             type = "quant", N = Nexp, MAF = m$maf.e)
  d2 <- list(beta = m$beta.o, varbeta = m$se.o^2, snp = m$rsid, position = m$position.o,
             type = "quant", N = Nout, MAF = m$maf.o)
  if (!is.na(sdY_forced)) { d1$sdY <- sdY_forced; d2$sdY <- sdY_forced }
  r <- tryCatch(suppressWarnings(coloc.abf(d1, d2, p12 = p12)),
                error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(NULL)
  s <- r$summary
  topp <- function(h) {
    pp <- r$results[[paste0("SNP.PP.H", h)]]
    if (is.null(pp) || all(is.na(pp))) return("-")
    i <- which.max(pp); sprintf("%s(%.2f)", r$results$snp[i], pp[i])
  }
  # strongest signal for exposure/outcome in this region, to help judge "whether both sides carry a signal"
  lead_e <- m$rsid[which.min(m$p.e)]; lead_o <- m$rsid[which.min(m$p.o)]
  cat(sprintf("%-46s nSNP=%-5d H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.3f | leadE %s p=%.1e / leadO %s p=%.1e\n",
              tag, s[["nsnps"]], s[["PP.H0.abf"]], s[["PP.H1.abf"]], s[["PP.H2.abf"]],
              s[["PP.H3.abf"]], s[["PP.H4.abf"]], lead_e, min(m$p.e), lead_o, min(m$p.o)))
  data.frame(tag = tag, p12 = p12, maf_cut = maf_cut, sdY_forced = sdY_forced,
             nsnp = s[["nsnps"]], PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]],
             PP.H2 = s[["PP.H2.abf"]], PP.H3 = s[["PP.H3.abf"]], PP.H4 = s[["PP.H4.abf"]],
             lead_exposure = lead_e, lead_outcome = lead_o,
             minP_exposure = min(m$p.e), minP_outcome = min(m$p.o),
             top_H4_snp = topp(4), stringsAsFactors = FALSE)
}

cat(strrep("=", 108), "\n", sep = "")
cat("Colocalization: IL6R / SAA x bilirubin   ", format(Sys.time()), "\n", sep = "")
cat(strrep("=", 108), "\n\n", sep = "")

res <- list(); add <- function(x) if (!is.null(x)) res[[length(res) + 1]] <<- x

cat("-- A. positive control: IL6R x CRP (shared causal variant already established)", strrep("-", 21), "\n", sep = "")
add(run_coloc("exp_IL6R_folkersen2018", "out_crp_said2022", "IL6R(F2018) x CRP [pos ctrl]", 3394, 575531))
cat("\n-- B. main: IL6R x total bilirubin", strrep("-", 56), "\n", sep = "")
add(run_coloc("exp_IL6R_folkersen2018", "out_bil_ukb_neale",  "IL6R(F2018) x Bil UKB-Neale", 3394, 342829))
add(run_coloc("exp_IL6R_folkersen2020", "out_bil_ukb_neale",  "IL6R(F2020) x Bil UKB-Neale", 21758, 342829))
add(run_coloc("exp_IL6R_pietzner2020",  "out_bil_ukb_neale",  "IL6R(P2020) x Bil UKB-Neale", 10708, 342829))
add(run_coloc("exp_IL6R_folkersen2018", "out_bil_barton2021", "IL6R(F2018) x Bil Barton21", 3394, 436055))
add(run_coloc("exp_IL6R_folkersen2018", "out_bil_bbj",        "IL6R(F2018) x Bil BBJ EastAsia", 3394, 110207))
cat("\n-- C. negative control: SAA1 x bilirubin (MR null)", strrep("-", 40), "\n", sep = "")
add(run_coloc("exp_SAA1_pietzner2020", "out_bil_ukb_neale_SAAregion",  "SAA1(P2020) x Bil UKB-Neale", 10708, 342829))
add(run_coloc("exp_SAA1_pietzner2020", "out_bil_barton2021_SAAregion", "SAA1(P2020) x Bil Barton21", 10708, 436055))
cat("\n-- D. sensitivity: p12 prior (best-powered main pair)", strrep("-", 38), "\n", sep = "")
for (pp in c(1e-4, 1e-6, 1e-7)) {
  add(run_coloc("exp_IL6R_folkersen2018", "out_bil_ukb_neale",
                sprintf("IL6R(F2018) x UKB-Neale [p12=%.0e]", pp), 3394, 342829, p12 = pp))
}
cat("\n-- E. sensitivity: sdY fixed at 1 (standardised traits)", strrep("-", 35), "\n", sep = "")
add(run_coloc("exp_IL6R_folkersen2018", "out_bil_ukb_neale", "IL6R(F2018) x UKB-Neale [sdY=1]", 3394, 342829, sdY_forced = 1))
add(run_coloc("exp_IL6R_folkersen2018", "out_crp_said2022",  "IL6R(F2018) x CRP [sdY=1, pos ctrl]", 3394, 575531, sdY_forced = 1))
cat("\n-- F. sensitivity: MAF > 5%", strrep("-", 56), "\n", sep = "")
add(run_coloc("exp_IL6R_folkersen2018", "out_bil_ukb_neale", "IL6R(F2018) x UKB-Neale [MAF>5%]", 3394, 342829, maf_cut = 0.05))
cat("\n-- G. sliding window centred on the IL6R lead variant rs4129267 (chr1:154,426,264)", strrep("-", 12), "\n", sep = "")
CTR <- 154426264
for (hw in c(100000, 50000, 25000, 10000)) {
  add(run_coloc("exp_IL6R_folkersen2018", "out_bil_ukb_neale",
                sprintf("IL6R(F2018) x UKB-Neale [+-%dkb]", hw / 1000), 3394, 342829,
                win_center = CTR, win_half = hw))
}
for (hw in c(100000, 50000, 25000, 10000)) {
  add(run_coloc("exp_IL6R_folkersen2018", "out_crp_said2022",
                sprintf("IL6R(F2018) x CRP [+-%dkb, pos ctrl]", hw / 1000), 3394, 575531,
                win_center = CTR, win_half = hw))
}
cat("\n")

if (length(res) > 0) {
  out <- do.call(rbind, res)
  write.csv(out, OUT, row.names = FALSE, fileEncoding = "UTF-8")
  cat(strrep("=", 108), "\n", sep = "")
  cat("saved:", OUT, "\n")
  cat(sprintf("PP.H4 > 0.8 : %d / %d ;  PP.H3 > 0.8 : %d / %d ;  H4 vs H3 ratio > 3 : %d\n",
              sum(out$PP.H4 > 0.8), nrow(out), sum(out$PP.H3 > 0.8), nrow(out),
              sum(out$PP.H4 / pmax(out$PP.H3, 1e-12) > 3)))
  cat("\nInterpretation\n  high PP.H4 = shared causal variant (supports the MR causal reading)\n",
      " high PP.H3 = two distinct signals (pleiotropy, or an intermediate step)\n",
      " high PP.H1/H2 = signal in one trait only\n",
      " nSNP < 200 -> underpowered; a 50/50 H3/H4 split means 'uninformative', not 'no colocalization'\n", sep = "")
} else cat("!! no results\n")
