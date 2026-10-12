## Final version of coloc.abf: allele harmonization + global MAF donation + narrow-window / full-region comparison
## Main analysis uses exp_IL6R_folkersen2018 (the only unsaturated one: min p=2.19e-265, number of variants tied at the minimum p = 1)
suppressMessages(library(coloc))
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]; OUT <- args[2]; MODE <- args[3]; EXPSEL <- args[4]

rd <- function(f) {
  p <- file.path(DIR, paste0(f, ".csv"))
  if (!file.exists(p)) return(NULL)
  d <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  names(d)[1] <- gsub("^[^A-Za-z]+", "", names(d)[1])
  for (c in c("position", "eaf", "beta", "se", "p"))
    if (c %in% names(d)) d[[c]] <- suppressWarnings(as.numeric(d[[c]]))
  d <- d[!is.na(d$beta) & !is.na(d$se) & d$p > 0 & d$se > 0 &
         !is.na(d$rsid) & nzchar(as.character(d$rsid)), ]
  d$maf <- pmin(d$eaf, 1 - d$eaf)
  d <- d[order(d$p), ]
  d[!duplicated(d$rsid), ]
}

## ---- global MAF donation: frequency is a property of a variant, so any record with the same rsID can lend it ----
GLOBAL <- new.env(parent = emptyenv())
for (f in list.files(DIR, pattern = "^(exp_|out_).*\\.csv$")) {
  d <- tryCatch(read.csv(file.path(DIR, f), stringsAsFactors = FALSE, check.names = FALSE),
                error = function(x) NULL)
  if (is.null(d) || !all(c("rsid", "eaf") %in% names(d))) next
  d$eaf <- suppressWarnings(as.numeric(d$eaf))
  for (i in which(!is.na(d$eaf))) {
    k <- as.character(d$rsid[i])
    if (!exists(k, envir = GLOBAL, inherits = FALSE))
      assign(k, min(d$eaf[i], 1 - d$eaf[i]), envir = GLOBAL)
  }
}
printf <- function(...) cat(sprintf(...), "\n")

EXPS <- list(exp_IL6R_folkersen2018 = 3394, exp_IL6R_folkersen2020 = 21758,
             exp_IL6R_pietzner2020 = 10708, exp_IL6R_sun2020 = 3301)
if (EXPSEL != "all") EXPS <- EXPS[names(EXPS) %in% strsplit(EXPSEL, ",")[[1]]]

OUTS <- list(out_MCH_vuckovic = 408112, out_MCV_vuckovic = 408112, out_MCHC_vuckovic = 408112,
             out_RDW_vuckovic = 408112, out_RBC_vuckovic = 408112, out_RETIC_vuckovic = 408112,
             out_MRV_vuckovic = 408112, out_PLT_vuckovic = 408112,
             out_BIL_barton2021 = 436055, out_CRP_said2022 = 575531)

cat(strrep("=", 126), "\n", sep = "")
cat("coloc.abf  MODE=", MODE, "  EXP=", paste(names(EXPS), collapse = ","), "\n", sep = "")
cat(strrep("=", 126), "\n", sep = "")

res <- list()
for (ek in names(EXPS)) {
  e <- rd(ek); if (is.null(e)) next
  for (ok in names(OUTS)) {
    o <- rd(ok); if (is.null(o)) next
    m <- merge(e, o, by = "rsid", suffixes = c(".e", ".o"))
    ## * allele harmonization: the betas on both sides must point to the same allele, otherwise H4 is completely distorted
    same <- (m$ea.e == m$ea.o) & (m$nea.e == m$nea.o)
    flip <- (m$ea.e == m$nea.o) & (m$nea.e == m$ea.o)
    n_conf <- sum(!(same | flip))
    m <- m[same | flip, ]
    if (any(flip[same | flip])) {
      f2 <- flip[same | flip]
      m$beta.o[f2] <- -m$beta.o[f2]
    }
    ## fill in missing MAF
    for (i in which(is.na(m$maf.o))) if (exists(m$rsid[i], envir = GLOBAL, inherits = FALSE)) m$maf.o[i] <- get(m$rsid[i], envir = GLOBAL)
    for (i in which(is.na(m$maf.e))) if (exists(m$rsid[i], envir = GLOBAL, inherits = FALSE)) m$maf.e[i] <- get(m$rsid[i], envir = GLOBAL)
    n_before <- nrow(m); m <- m[!is.na(m$maf.o) & !is.na(m$maf.e), ]
    printf("%-40s MAF dropped %d", paste(sub("^exp_", "", ek), sub("^out_", "", ok), sep = " x "), n_before - nrow(m))
    if (nrow(m) < 25) { printf("   -> skip (n=%d)", nrow(m)); next }
    d1 <- list(beta = m$beta.e, varbeta = m$se.e^2, snp = m$rsid, type = "quant",
               N = EXPS[[ek]], MAF = m$maf.e)
    d2 <- list(beta = m$beta.o, varbeta = m$se.o^2, snp = m$rsid, type = "quant",
               N = OUTS[[ok]], MAF = m$maf.o)
    r <- tryCatch(suppressWarnings(coloc.abf(d1, d2)),
                  error = function(x) { printf("   coloc error: %s", conditionMessage(x)); NULL })
    if (is.null(r)) next
    s <- r$summary
    printf("   n=%-4d flip=%-3d conf=%-3d | H3=%.3g  **H4=%.3f**  leadE=%-12s leadO=%-12s",
           nrow(m), sum(flip[same | flip]), n_conf,
           s[["PP.H3.abf"]], s[["PP.H4.abf"]],
           m$rsid[which.min(m$p.e)], m$rsid[which.min(m$p.o)])
    res[[length(res) + 1]] <- data.frame(exp = ek, out = ok, mode = MODE, nsnp = nrow(m),
      n_flip = sum(flip), n_conflict = n_conf,
      H0 = s[["PP.H0.abf"]], H1 = s[["PP.H1.abf"]], H2 = s[["PP.H2.abf"]],
      H3 = s[["PP.H3.abf"]], H4 = s[["PP.H4.abf"]],
      leadE = m$rsid[which.min(m$p.e)], leadO = m$rsid[which.min(m$p.o)], stringsAsFactors = FALSE)
  }
}
if (length(res)) {
  dd <- do.call(rbind, res); write.csv(dd, OUT, row.names = FALSE)
  printf("\nsaved %s (%d rows); H4>0.8: %d", OUT, nrow(dd), sum(dd$H4 > 0.8))
} else printf("!! no results")
