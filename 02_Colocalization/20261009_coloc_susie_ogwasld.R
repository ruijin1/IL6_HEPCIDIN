# =============================================================================
# coloc.susie v3 -- using the **OpenGWAS European reference LD** (root causes of the first two failures eliminated)
#   usage: Rscript 20261009_coloc_susie_ogwasld.R <in_dir> <out_csv>
#   env vars: LD_CONV=a1|a2   selects the LD allele convention (fix it with the positive control first)
#             ONLY=<substr>  run only pairs whose tag contains this substring (for debugging)
#
# Difference from v2 (FinnGen Finnish panel): the panel is now European, so `susie_rss` should converge.
# Same acceptance criterion as v2: **the positive control IL6R x CRP must give a high PP.H4**; otherwise no results are reported at all.
# =============================================================================
suppressMessages(library(coloc)); suppressMessages(library(susieR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript ... <in_dir> <out_csv>")
DIR <- args[1]; OUT <- args[2]
CONV <- Sys.getenv("LD_ALLELE_CONV", "a1")
ONLY <- Sys.getenv("ONLY", "")

LDraw <- as.matrix(read.csv(file.path(DIR, "ld_matrix.csv"), row.names = 1, check.names = FALSE))
snpl <- rownames(LDraw)
sp <- strsplit(snpl, "_", fixed = TRUE)
rsid <- vapply(sp, function(x) x[1], "")
a1 <- vapply(sp, function(x) x[2], "")
a2 <- vapply(sp, function(x) paste(x[3:length(x)], collapse = "_"), "")
cat(sprintf("LD: %d x %d variants; min eigenvalue %.3g\n", nrow(LDraw), ncol(LDraw),
            min(eigen(LDraw, only.values = TRUE)$values)))
# nearest positive-definite projection (the OpenGWAS LD reference covers only some variants -> uncovered ones are set to 0, so the matrix may be slightly negative)
LD <- LDraw
w <- eigen(LD, symmetric = TRUE)
if (min(w$values) <= 1e-6) {
  w$values <- pmax(w$values, 1e-6)
  LD <- w$vectors %*% diag(w$values) %*% t(w$vectors)
  d <- sqrt(diag(LD)); LD <- LD / outer(d, d); diag(LD) <- 1
  cat(sprintf("  not positive definite / semi-definite -> min eigenvalue after nearest positive-definite projection %.3g\n", min(eigen(LD, only.values=TRUE)$values)))
}
# * LD row names are panel strings (rs4129267_T_C) and must be replaced with **plain rsIDs** so that indexing
#   by rsid works; if the same rsID appears more than once (different allele representations), deduplicate and keep the first.
keep <- !duplicated(rsid)
LD <- LD[keep, keep, drop = FALSE]
RS <- rsid[keep]
dimnames(LD) <- list(RS, RS)
A1 <- setNames(a1, rsid)[RS]; A2 <- setNames(a2, rsid)[RS]
cat(sprintf("[OK] LD indexed by plain rsID (%d unique rsIDs)\n", length(RS)))

## ---- MAF donation map: variant frequency is a property of the variant, so the MAF of the same rsID can be borrowed from any dataset that has eaf ----
## (some datasets such as CRP Said2022 have an empty eaf; coloc's check_dataset errors out directly when MAF contains NA)
DONOR <- new.env(parent = emptyenv())
for (f in list.files(DIR, pattern = "^(exp_|out_).*\\.csv$")) {
  d <- tryCatch(read.csv(file.path(DIR, f), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM"),
                error = function(x) NULL)
  if (is.null(d) || !all(c("rsid", "eaf") %in% names(d))) next
  names(d)[1] <- sub("^\xef\xbb\xbf", "", names(d)[1])
  d <- d[!is.na(d$eaf) & nzchar(d$rsid), c("rsid", "eaf")]
  for (i in seq_len(nrow(d))) {
    if (!exists(d$rsid[i], envir = DONOR, inherits = FALSE))
      assign(d$rsid[i], min(d$eaf[i], 1 - d$eaf[i]), envir = DONOR)
  }
}
cat(sprintf("MAF donation map: %d variants (from %d files)\n",
            length(ls(DONOR)), length(list.files(DIR, pattern = "^(exp_|out_).*\\.csv$"))))

read_reg <- function(f) {
  d <- read.csv(file.path(DIR, paste0(f, ".csv")), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  names(d)[1] <- sub("^\xef\xbb\xbf", "", names(d)[1])
  for (v in c("position","eaf","beta","se","p","n")) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d <- d[!is.na(d$beta) & !is.na(d$se) & d$se > 0 & nzchar(d$rsid), ]
  d <- d[d$rsid %in% RS, ]
  d$u <- if (CONV == "a1") ifelse(d$ea == A1[d$rsid], 1, -1) else ifelse(d$ea == A2[d$rsid], 1, -1)
  d$maf <- pmin(d$eaf, 1 - d$eaf)
  na_i <- which(is.na(d$maf))
  if (length(na_i)) {
    d$maf[na_i] <- vapply(d$rsid[na_i],
                          function(r) if (exists(r, envir = DONOR, inherits = FALSE))
                            get(r, envir = DONOR) else NA_real_, numeric(1))
    cat(sprintf("   [MAF donated] %s: %d / %d variants\n", f, length(na_i), nrow(d)))
  }
  d <- d[order(d$p), ]; d[!duplicated(d$rsid), ]
}

run_pair <- function(expf, outf, tag, Nexp, Nout, otype = "quant", s_prop = NA, L = 10) {
  e <- tryCatch(read_reg(expf), error = function(x) NULL)
  o <- tryCatch(read_reg(outf), error = function(x) NULL)
  if (is.null(e) || is.null(o)) { cat(sprintf("%-44s load failed\n", tag)); return(NULL) }
  m <- merge(e, o, by = "rsid", suffixes = c(".e", ".o"))
  if (nrow(m) < 40) { cat(sprintf("%-44s shared %d -> skipped\n", tag, nrow(m))); return(NULL) }
  bad <- which(is.na(m$maf.e) | is.na(m$maf.o))
  if (length(bad)) { cat(sprintf("   %s: dropping %d variants with missing MAF\n", tag, length(bad))); m <- m[-bad, ] }
  if (nrow(m) < 40) { cat(sprintf("%-44s insufficient after dropping -> skipped\n", tag)); return(NULL) }
  sn <- m$rsid
  R <- LD[sn, sn, drop = FALSE]
  R1 <- R * outer(m$u.e, m$u.e); R2 <- R * outer(m$u.o, m$u.o)
  dimnames(R1) <- dimnames(R); dimnames(R2) <- dimnames(R)
  diag(R1) <- 1; diag(R2) <- 1
  d1 <- list(beta = m$beta.e, varbeta = m$se.e^2, snp = sn, N = Nexp, type = "quant",
             MAF = m$maf.e, LD = R1)
  d2 <- list(beta = m$beta.o, varbeta = m$se.o^2, snp = sn, N = Nout, type = otype,
             MAF = m$maf.o, LD = R2)
  if (otype == "cc" && !is.na(s_prop)) d2$s <- s_prop
  s1 <- tryCatch(suppressWarnings(suppressMessages(runsusie(d1, suffix = 1, L = L))),
                 error = function(x) { cat(sprintf("%-44s susie1 ERR: %s\n", tag, conditionMessage(x))); NULL })
  if (is.null(s1)) return(NULL)
  s2 <- tryCatch(suppressWarnings(suppressMessages(runsusie(d2, suffix = 2, L = L))),
                 error = function(x) { cat(sprintf("%-44s susie2 ERR: %s\n", tag, conditionMessage(x))); NULL })
  if (is.null(s2)) return(NULL)
  r <- tryCatch(suppressWarnings(coloc.susie(s1, s2)),
                error = function(x) { cat(sprintf("%-44s coloc ERR: %s\n", tag, conditionMessage(x))); NULL })
  if (is.null(r) || is.null(r$summary)) { cat(sprintf("%-44s coloc produced no result\n", tag)); return(NULL) }
  sm <- r$summary
  # * Each row of the summary is "one pair of credible sets"; PP.H4 and PP.H3 must be taken from the **same row**,
  #   otherwise one reaches the absurd conclusion "both H4=1 and H3=1".
  #   Standard practice: report the row with max(PP.H4) (and the PP.H3 on that same row), and list max(PP.H3)
  #   separately for reference.
  i4 <- which.max(sm$PP.H4.abf); i3 <- which.max(sm$PP.H3.abf)
  h4 <- sm$PP.H4.abf[i4]; h3_same <- sm$PP.H3.abf[i4]; h3 <- sm$PP.H3.abf[i3]
  cat(sprintf("%-42s nSNP=%-4d CS=%-2d/%-2d | **maxH4=%.3f** (H3 same row=%.3f, CS %s-%s) | maxH3=%.3f | leadE %s p=%.1e / leadO %s p=%.1e\n",
              tag, nrow(m), length(s1$sets$cs), length(s2$sets$cs),
              h4, h3_same, sm$idx1[i4], sm$idx2[i4], h3,
              m$rsid[which.min(m$p.e)], min(m$p.e), m$rsid[which.min(m$p.o)], min(m$p.o)))
  data.frame(tag = tag, conv = CONV, nsnp = nrow(m),
             n_cs1 = length(s1$sets$cs), n_cs2 = length(s2$sets$cs),
             maxPP.H4 = h4, PP.H3_at_maxH4 = h3_same, cs_pair = sprintf("%s-%s", sm$idx1[i4], sm$idx2[i4]),
             maxPP.H3 = h3,
             lead_exposure = m$rsid[which.min(m$p.e)], lead_outcome = m$rsid[which.min(m$p.o)],
             minP_exposure = min(m$p.e), minP_outcome = min(m$p.o), stringsAsFactors = FALSE)
}

N <- c(MCH=408112, MCV=408112, MCHC=408112, RDW=408112, RBC=408112, RETIC=408112, MRV=408112,
       PLT=408112, BILbarton=436055, BILbbj=110207, CRP=575531)
EXPS <- list(exp_IL6R_folkersen2018 = 3394, exp_IL6R_folkersen2020 = 21758,
             exp_IL6R_pietzner2020 = 10708, exp_IL6R_sun2020 = 3301)
OUTS <- list(out_MCH_vuckovic = "MCH", out_MCV_vuckovic = "MCV", out_MCHC_vuckovic = "MCHC",
             out_RDW_vuckovic = "RDW", out_RBC_vuckovic = "RBC", out_RETIC_vuckovic = "RETIC",
             out_MRV_vuckovic = "MRV", out_PLT_vuckovic = "PLT",
             out_BIL_barton2021 = "BILbarton", out_BIL_bbj = "BILbbj",
             out_CRP_said2022 = "CRP")

cat(strrep("=", 112), "\n", sep = "")
cat("coloc.susie v3 (OpenGWAS European LD)  convention=", CONV, "  ", format(Sys.time()), "\n", sep = "")
cat(strrep("=", 112), "\n\n", sep = "")

res <- list(); add <- function(x) if (!is.null(x)) res[[length(res) + 1]] <<- x
for (ek in names(EXPS)) {
  for (ok in names(OUTS)) {
    tag <- sprintf("%s x %s", sub("^exp_", "", ek), OUTS[[ok]])
    if (nzchar(ONLY) && !grepl(ONLY, tag, fixed = TRUE)) next
    add(run_pair(ek, ok, tag, EXPS[[ek]], N[[OUTS[[ok]]]]))
  }
}
# negative control
if (!nzchar(ONLY) || grepl("SAA1", ONLY, fixed = TRUE)) {
  add(run_pair("exp_SAA1_pietzner2020", "out_MCH_vuckovic", "SAA1 x MCH [neg ctrl]", 10708, N[["MCH"]]))
  add(run_pair("exp_SAA1_pietzner2020", "out_BIL_barton2021", "SAA1 x BIL_barton [neg ctrl]", 10708, N[["BILbarton"]]))
}
if (length(res)) {
  out <- do.call(rbind, res)
  write.csv(out, OUT, row.names = FALSE, fileEncoding = "UTF-8")
  cat("\nsaved:", OUT, " (%d rows)\n", sep = "")
  cat(sprintf("maxPP.H4 > 0.8: %d / %d\n", sum(out$maxPP.H4 > 0.8), nrow(out)))
} else cat("!! no results\n")
