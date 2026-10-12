# =============================================================================
# coloc.susie:  IL6R pQTL x MCH / MCV / bilirubin / neonatal jaundice / CRP
#   usage: Rscript 20261009_coloc_susie.R <input_dir> <output_csv>
#
# WHY susie instead of abf
#   No LD-based method could be used in round 1 because the OpenGWAS LD endpoint is down.
#   In round 2, coloc.abf produced PP.H3 = 1.000 for EVERY pair -- including the positive
#   control IL6R x CRP, in which the shared causal variant (IL6R missense rs2228145) is
#   textbook. Diagnosis: the IL6R pQTL has p-values saturated at 1e-265 (z up to 38) and its
#   peak (chr1:154,426,264 GRCh37) sits ~20-30 kb from the CRP peak (chr1:154,395,212).
#   coloc.abf assumes ONE causal variant and does NOT use LD, so it cannot reconcile two
#   peaks in the same LD block and calls them "distinct" (H3).
#   => the abf results are UNINTERPRETABLE and must not be reported as "no colocalisation".
#   This script uses the FinnGen R12 LD reference panel (520,210 samples, tabix) to run
#   susie_rss + coloc.susie. If the positive control passes here, the negatives become usable.
#
# LD reference caveat (must be declared in the manuscript):
#   the LD long-table lists only part of all pairs, so the assembled matrix had to be
#   projected to the nearest positive-definite correlation matrix (min eigenvalue
#   -2.75 -> 9e-05, documented in the build script log). The panel is Finnish; exposure and
#   outcome GWAS are mostly European.
# =============================================================================
suppressMessages(library(coloc))
suppressMessages(library(susieR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript 20261009_coloc_susie.R <input_dir> <output_csv>")
DIR <- args[1]; OUT <- args[2]

## ---- LD matrix + variant map -------------------------------------------------
LD <- as.matrix(read.csv(file.path(DIR, "ld_matrix_il6r.csv"), row.names = 1,
                         check.names = FALSE))
cat("LD matrix:", nrow(LD), "x", ncol(LD), " min eig:", min(eigen(LD, only.values = TRUE)$values), "\n")
MAP <- read.csv(file.path(DIR, "ld_variant_map.csv"), stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")

# * Allele direction: the FinnGen LD id looks like chr1_154426264_T_C (a1_a2),
#   and the r from PLINK2 `--r-unphased` is computed with respect to **ALT (= the 2nd allele a2 of the id)**.
#   The r must be rotated onto the scale of "this dataset's effect allele ea", otherwise susie reports
#   "The estimated prior variance is unreasonably large ... mismatch between the summary
#    statistics and the LD matrix" (hit this in practice).
#   CONV: "a2" = the panel's counted allele is the 2nd allele of the id (default, PLINK2 semantics); "a1" = the 1st.
CONV <- Sys.getenv("LD_ALLELE_CONV", "a2")
sp <- strsplit(MAP$ld_id, "_", fixed = TRUE)
MAP$a1 <- vapply(sp, function(x) x[3], "")
MAP$a2 <- vapply(sp, function(x) paste(x[4:length(x)], collapse = "_"), "")
MAP$u <- if (CONV == "a2") ifelse(MAP$ea == MAP$a2, 1, -1) else ifelse(MAP$ea == MAP$a1, 1, -1)
cat("LD allele convention:", CONV, " (u=+1 for", sum(MAP$u == 1), "of", nrow(MAP), "records)\n")

read_reg <- function(f) {
  d <- read.csv(file.path(DIR, paste0(f, ".csv")), stringsAsFactors = FALSE,
                fileEncoding = "UTF-8-BOM")
  names(d)[1] <- sub("^\xef\xbb\xbf", "", names(d)[1])
  for (v in c("position", "eaf", "beta", "se", "p", "n")) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d <- d[!is.na(d$beta) & !is.na(d$se) & d$se > 0 & !is.na(d$rsid) & nzchar(d$rsid), ]
  d <- d[order(d$p), ]
  d[!duplicated(d$rsid), ]
}

## per-file: id -> (ld_id, sign)
orient <- function(d, f) {
  mp <- MAP[MAP$file == f, ]
  i <- match(d$rsid, mp$id)
  keep <- which(!is.na(i))
  d <- d[keep, ]; i <- i[keep]
  d$ld_id <- mp$ld_id[i]; d$u <- mp$u[i]
  d <- d[order(d$p), ]
  d[!duplicated(d$ld_id), ]
}

pair_data <- function(expf, outf) {
  e <- orient(read_reg(expf), expf)
  o <- orient(read_reg(outf), outf)
  merge(e, o, by = "ld_id", suffixes = c(".e", ".o"))
}

run_susie_pair <- function(expf, outf, tag, Nexp, Nout, otype = "quant", s_prop = NA,
                           L = 10, p12 = 1e-5) {
  m <- tryCatch(pair_data(expf, outf), error = function(e) { cat("  load error\n"); NULL })
  if (is.null(m) || nrow(m) < 30) {
    cat(sprintf("%-48s shared=%d -> skipped\n", tag, if (is.null(m)) 0 else nrow(m))); return(NULL)
  }
  sn <- m$ld_id
  R <- LD[sn, sn, drop = FALSE]
  # * Rotate the LD onto the scale of each dataset's "effect allele ea":
  #   R_ea[i,j] = R_panel[i,j] * u_i * u_j ; beta needs no change (it is per-ea by construction)
  uu <- outer(m$u.e, m$u.o)
  Rm <- R * uu
  dimnames(Rm) <- dimnames(R)
  # separate R for exposure and outcome sides (same u; identical if both sides use the same file convention)
  R1 <- R * outer(m$u.e, m$u.e); dimnames(R1) <- dimnames(R)
  R2 <- R * outer(m$u.o, m$u.o); dimnames(R2) <- dimnames(R)
  diag(R1) <- 1; diag(R2) <- 1
  d1 <- list(beta = m$beta.e, varbeta = m$se.e^2, snp = sn, N = Nexp, type = "quant",
             MAF = pmin(m$eaf.e, 1 - m$eaf.e), LD = R1)
  d2 <- list(beta = m$beta.o, varbeta = m$se.o^2, snp = sn, N = Nout, type = otype,
             MAF = pmin(m$eaf.o, 1 - m$eaf.o), LD = R2)
  if (otype == "cc" && !is.na(s_prop)) d2$s <- s_prop
  s1 <- tryCatch(suppressWarnings(suppressMessages(runsusie(d1, suffix = 1, L = L))),
                 error = function(e) { cat("  susie1 error:", conditionMessage(e), "\n"); NULL })
  if (is.null(s1)) return(NULL)
  s2 <- tryCatch(suppressWarnings(suppressMessages(runsusie(d2, suffix = 2, L = L))),
                 error = function(e) { cat("  susie2 error:", conditionMessage(e), "\n"); NULL })
  if (is.null(s2)) return(NULL)
  r <- tryCatch(suppressWarnings(coloc.susie(s1, s2)),
                error = function(e) { cat("  coloc.susie error:", conditionMessage(e), "\n"); NULL })
  if (is.null(r) || is.null(r$summary)) { cat(sprintf("%-48s coloc.susie returned NULL\n", tag)); return(NULL) }
  sm <- r$summary
  cat(sprintf("%-48s nSNP=%-5d n1=%-3d n2=%-3d H4=%-7s H3=%-7s | leadE %s p=%.1e / leadO %s p=%.1e\n",
              tag, nrow(m), length(s1$sets$cs), length(s2$sets$cs),
              paste(sprintf("%.3f", sm$PP.H4.abf), collapse = "/"),
              paste(sprintf("%.3f", sm$PP.H3.abf), collapse = "/"),
              m$rsid.e[which.min(m$p.e)], min(m$p.e), m$rsid.o[which.min(m$p.o)], min(m$p.o)))
  data.frame(tag = tag, nsnp = nrow(m), n_cs1 = length(s1$sets$cs), n_cs2 = length(s2$sets$cs),
             idx1 = paste(sm$idx1, collapse = ";"), idx2 = paste(sm$idx2, collapse = ";"),
             PP.H4 = paste(sprintf("%.4f", sm$PP.H4.abf), collapse = ";"),
             PP.H3 = paste(sprintf("%.4f", sm$PP.H3.abf), collapse = ";"),
             max_PP.H4 = max(sm$PP.H4.abf),
             lead_exposure = m$rsid.e[which.min(m$p.e)], lead_outcome = m$rsid.o[which.min(m$p.o)],
             stringsAsFactors = FALSE)
}

cat(strrep("=", 112), "\n", sep = "")
cat("coloc.susie  IL6R pQTL x MCH / MCV / bilirubin / neonatal jaundice / CRP   ",
    format(Sys.time()), "\n", sep = "")
cat(strrep("=", 112), "\n\n", sep = "")

res <- list(); add <- function(x) if (!is.null(x)) res[[length(res) + 1]] <<- x

cat("-- A. POSITIVE CONTROL: IL6R x CRP  (must show high PP.H4)\n")
add(run_susie_pair("exp_IL6R_folkersen2018", "out_crp_said2022", "IL6R(F2018) x CRP [pos ctrl]", 3394, 575531))

cat("\n-- B. IL6R x MCH (FinnGen lab, n=447,004)\n")
add(run_susie_pair("exp_IL6R_folkersen2018", "out_mch_finngen", "IL6R(F2018) x MCH", 3394, 447004))
add(run_susie_pair("exp_IL6R_folkersen2020", "out_mch_finngen", "IL6R(F2020) x MCH", 21758, 447004))
add(run_susie_pair("exp_IL6R_pietzner2020",  "out_mch_finngen", "IL6R(P2020) x MCH", 10708, 447004))

cat("\n-- C. IL6R x MCV (FinnGen lab, n=447,040)\n")
add(run_susie_pair("exp_IL6R_folkersen2018", "out_mcv_finngen", "IL6R(F2018) x MCV", 3394, 447040))
add(run_susie_pair("exp_IL6R_folkersen2020", "out_mcv_finngen", "IL6R(F2020) x MCV", 21758, 447040))
add(run_susie_pair("exp_IL6R_pietzner2020",  "out_mcv_finngen", "IL6R(P2020) x MCV", 10708, 447040))

cat("\n-- D. IL6R x total bilirubin (FinnGen quantitative lab, n=190,407)\n")
add(run_susie_pair("exp_IL6R_folkersen2018", "out_bil_finngen_lab", "IL6R(F2018) x Bil FinnGen-lab", 3394, 190407))
add(run_susie_pair("exp_IL6R_folkersen2020", "out_bil_finngen_lab", "IL6R(F2020) x Bil FinnGen-lab", 21758, 190407))

cat("\n-- E. IL6R x neonatal jaundice (FinnGen P16, case-control, ncase=499)\n")
add(run_susie_pair("exp_IL6R_folkersen2018", "out_jaund_finngen", "IL6R(F2018) x Neonatal jaundice", 3394, 500187,
                   otype = "cc", s_prop = 499 / 500187))
add(run_susie_pair("exp_IL6R_folkersen2020", "out_jaund_finngen", "IL6R(F2020) x Neonatal jaundice", 21758, 500187,
                   otype = "cc", s_prop = 499 / 500187))

cat("\n-- F. sensitivity: L = 5\n")
add(run_susie_pair("exp_IL6R_folkersen2018", "out_crp_said2022", "IL6R(F2018) x CRP [L=5, pos ctrl]", 3394, 575531, L = 5))
add(run_susie_pair("exp_IL6R_folkersen2018", "out_mch_finngen", "IL6R(F2018) x MCH [L=5]", 3394, 447004, L = 5))
cat("\n")

if (length(res) > 0) {
  out <- do.call(rbind, res)
  write.csv(out, OUT, row.names = FALSE, fileEncoding = "UTF-8")
  cat(strrep("=", 112), "\n", sep = "")
  cat("saved:", OUT, "\n")
  cat(sprintf("PP.H4 > 0.8 : %d / %d\n", sum(out$max_PP.H4 > 0.8), nrow(out)))
} else cat("!! no results\n")
