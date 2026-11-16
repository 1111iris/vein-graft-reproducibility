# 15_reproducibility_audit.R — teacher-mandated audit of concordance claims
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# (1) effect-size correlation + CI, (2) threshold sensitivity, (3) baselines,
# (4) regression slope, (5) clinical concordance CI + sensitivity analyses
suppressPackageStartupMessages({})  # inline Wilson CI used below
outdir <- "results"
d1 <- read.csv(file.path(outdir,"de_GSE310013_full.csv"))   # discovery (Augusta)
d2 <- read.csv(file.path(outdir,"de_GSE326710_full.csv"))   # replication (Lausanne)
occ <- read.csv(file.path(outdir,"de_GSE241205_fpkm_full.csv"))
core_all <- read.csv(file.path(outdir,"robust_core_priority_all.csv"))

m <- merge(d1[,c("gene","logFC","P.Value","adj.P.Val")],
           d2[,c("gene","logFC","P.Value","adj.P.Val")],
           by="gene", suffixes=c("_disc","_rep"))
cat("common genes tested in both cohorts:", nrow(m), "\n")

## (1) effect-size correlation across ALL common genes, with 95% CI
pt <- cor.test(m$logFC_disc, m$logFC_rep, method="pearson")
ps <- cor.test(m$logFC_disc, m$logFC_rep, method="spearman", exact=FALSE)
cat(sprintf("Pearson r = %.3f [%.3f, %.3f], p=%.3g\n",
    pt$estimate, pt$conf.int[1], pt$conf.int[2], pt$p.value))
cat(sprintf("Spearman rho = %.3f, p=%.3g\n", ps$estimate, ps$p.value))

## (4) regression slope (rep ~ disc) with CI
fit <- lm(logFC_rep ~ logFC_disc, data=m)
cf <- summary(fit)$coefficients
ci <- confint(fit)[2,]
cat(sprintf("slope = %.3f [%.3f, %.3f], intercept = %.3f, R2 = %.3f\n",
    cf[2,1], ci[1], ci[2], cf[1,1], summary(fit)$r.squared))

## (2) concordance across thresholds + (3) baselines

wilson <- function(k, n, z=1.96) {
  ph <- k/n; den <- 1 + z^2/n
  ctr <- (ph + z^2/(2*n))/den
  hw <- z*sqrt(ph*(1-ph)/n + z^2/(4*n^2))/den
  c(ctr - hw, ctr + hw)
}
conc <- function(x) mean(x > 0)
rows <- list()
add <- function(label, sel, denom) {
  if (is.null(sel)) return()
  k <- sum(sign(m$logFC_disc[sel]) == sign(m$logFC_rep[sel]))
  rows[[length(rows)+1]] <<- data.frame(
    criterion=label, n=sel_n <- sum(sel), concordant=k,
    rate=round(k/sel_n,4),
    wilson_lo=round(wilson(k, sel_n)[1],4),
    wilson_hi=round(wilson(k, sel_n)[2],4))
}
add("both FDR<0.05 (robust core)", m$adj.P.Val_disc<0.05 & m$adj.P.Val_rep<0.05)
add("discovery FDR<0.05", m$adj.P.Val_disc<0.05)
add("discovery FDR<0.01", m$adj.P.Val_disc<0.01)
add("discovery FDR<0.001", m$adj.P.Val_disc<0.001)
add("discovery p<0.05", m$P.Value_disc<0.05)
add("both p<0.05", m$P.Value_disc<0.05 & m$P.Value_rep<0.05)
add("discovery top 500 by FDR", rank(m$adj.P.Val_disc)<=500)
add("discovery top 2000 by FDR", rank(m$adj.P.Val_disc)<=2000)
add("ALL common genes (baseline)", rep(TRUE, nrow(m)))
add("non-significant in both (baseline)", m$adj.P.Val_disc>=0.05 & m$adj.P.Val_rep>=0.05)
sens <- do.call(rbind, rows)
# semantic assertions: every threshold analysis must have a non-empty, finite
# result. A script that "completes" with empty analyses must fail loudly, not
# pass silently (round-20 review: the nominal-p rows were silently empty).
stopifnot(all(sens$n > 0), all(is.finite(sens$rate)))
stopifnot(sens$n[sens$criterion == "both FDR<0.05 (robust core)"] == 4351,
          sens$n[sens$criterion == "discovery FDR<0.05"] == 9043)
cat("\n=== threshold sensitivity ===\n"); print(sens, row.names=FALSE)
write.csv(sens, file.path(outdir,"concordance_sensitivity.csv"), row.names=FALSE)

## heterogeneity: does the core differ in effect-size agreement from other significant genes?
core_set <- core_all$gene_symbol
m$core <- m$gene %in% core_set
h <- by(m, m$core, function(x) {
  dd <- abs(x$logFC_disc - x$logFC_rep)
  c(median_absdiff = median(dd), iqr = IQR(dd))
})
cat("\nabs(logFC_disc - logFC_rep): core vs non-core\n"); print(h)

## (5) clinical concordance audit
p <- core_all
p_occ <- p[!is.na(p$occ_logFC), ]
k <- sum(p_occ$occ_concordant); n <- nrow(p_occ)
w <- wilson(k, n)
bt <- binom.test(k, n, p=0.5)
cat(sprintf("\nclinical concordance: %d/%d = %.3f [Wilson %.3f, %.3f], vs 0.5: p=%.3g\n",
    k, n, k/n, w[1], w[2], bt$p.value))
# sensitivity: exclude immune-attributed genes
p_ni <- p_occ[which(p_occ$best_celltype != "Immune"), ]
k2 <- sum(p_ni$occ_concordant)
w2 <- wilson(k2, nrow(p_ni))
cat(sprintf("excluding immune-attributed: %d/%d = %.3f [%.3f, %.3f]\n",
    k2, nrow(p_ni), k2/nrow(p_ni), w2[1], w2[2]))
# sensitivity: exclude top 10% by |logFC_meta| (is it driven by big effects?)
thr <- quantile(abs(p_occ$logFC_meta), 0.9)
p_trim <- p_occ[abs(p_occ$logFC_meta) <= thr, ]
k3 <- sum(p_trim$occ_concordant)
w3 <- wilson(k3, nrow(p_trim))
cat(sprintf("excluding top 10%% by |logFC|: %d/%d = %.3f [%.3f, %.3f]\n",
    k3, nrow(p_trim), k3/nrow(p_trim), w3[1], w3[2]))
# baseline: concordance of ALL genes with occlusion data (not just core)
occ_lfc <- setNames(occ$logFC, occ$gene)
mm <- m[m$gene %in% names(occ_lfc), ]
kk <- sum(sign(mm$logFC_disc) == sign(occ_lfc[mm$gene]) & sign(mm$logFC_rep) == sign(occ_lfc[mm$gene]), na.rm=TRUE)
cat(sprintf("all-common-genes model-vs-occlusion concordance (both cohorts same sign as occ): %d/%d = %.3f\n",
    kk, nrow(mm), kk/nrow(mm)))
# simpler: discovery-significant genes vs occlusion
ds <- m[m$adj.P.Val_disc<0.05 & m$gene %in% names(occ_lfc), ]
kd <- sum(sign(ds$logFC_disc) == sign(occ_lfc[ds$gene]))
wd <- wilson(kd, nrow(ds))
cat(sprintf("discovery FDR<0.05 vs occlusion: %d/%d = %.3f [%.3f, %.3f]\n",
    kd, nrow(ds), kd/nrow(ds), wd[1], wd[2]))

## save audit table for manuscript
# core concordance taken dynamically from the threshold table (round-22 review:
# the summary previously hardcoded 0.953 / 4351, which could silently go stale)
core_row <- sens[sens$criterion == "both FDR<0.05 (robust core)", ]
stopifnot(core_row$n == 4351, core_row$concordant == 4148,
          abs(core_row$rate - round(core_row$concordant/core_row$n, 4)) < 1e-6)
audit_df <- data.frame(
  metric=c("pearson_r","pearson_ci_lo","pearson_ci_hi","spearman_rho","slope","slope_lo","slope_hi","R2",
           "core_conc","core_n","clin_conc","clin_n","clin_wilson_lo","clin_wilson_hi","clin_vs50_p"),
  value=c(pt$estimate, pt$conf.int, ps$estimate, cf[2,1], ci, summary(fit)$r.squared,
          core_row$rate, core_row$n, k/n, n, w[1], w[2], bt$p.value))
write.csv(audit_df, file.path(outdir,"reproducibility_audit_summary.csv"), row.names=FALSE)
cat("\naudit complete\n")
