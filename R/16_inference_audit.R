# 16_inference_audit.R — reviewer round 2: statistical inference-unit audit
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# (1) occlusion-free priority ranking (removes calibration-cohort circularity)
# (2) patient-level uncertainty for occlusion concordance (per-patient + leave-one-patient-out)
# (3) cell-attribution robustness (leave-one-marker-out + correlation-threshold sensitivity)
# (4) paired-design completeness check for both culture cohorts
suppressPackageStartupMessages({library(limma); library(Matrix)})
set.seed(20261029)
outdir <- "results"; d <- "data"

p <- read.csv(file.path(outdir,"robust_core_priority_all.csv"))
meta <- read.csv(file.path(outdir,"meta_analysis.csv"))

## ---------- (1) occlusion-free ranking ----------
# primary ranking uses ONLY culture-cohort meta statistics; occlusion never enters
p$score_culture_only <- with(p, abs(logFC_meta) * (-log10(pmax(FDR_meta, 1e-300))))
p <- p[order(-p$score_culture_only), ]
p$culture_rank <- seq_len(nrow(p))
top30_culture <- p$gene_symbol[1:30]
top30_old <- p$gene_symbol[order(-p$priority_score)][1:30]
rank_tab <- data.frame(
  rank = 1:30,
  top30_culture_only = top30_culture,
  top30_occlusion_penalized = top30_old,
  in_both = top30_culture %in% top30_old)
ov <- length(intersect(top30_culture, top30_old))
cat(sprintf("Top30 overlap culture-only vs occlusion-penalized: %d/30\n", ov))
write.csv(rank_tab, file.path(outdir,"top30_ranking_comparison.csv"), row.names=FALSE)
write.csv(p[, c("gene","gene_symbol","logFC_meta","FDR_meta","score_culture_only",
                "culture_rank","best_celltype","occ_logFC","occ_concordant","priority_score")],
          file.path(outdir,"robust_core_priority_culture_only.csv"), row.names=FALSE)

## ---------- (2) patient-level uncertainty for occlusion concordance ----------
x <- read.delim(gzfile(file.path(d,"GSE241205_gene_fpkm.txt.gz")), check.names=FALSE)
rownames(x) <- x$Symbol; x$Symbol <- NULL
x2 <- read.csv(gzfile(file.path(d,"GSE326710_Raw_counts_all_samples.csv.gz")), check.names=FALSE)
map <- x2[["Gene Symbol"]]; names(map) <- x2[["Gene ID"]]
sym <- unname(map[rownames(x)]); keepg <- !is.na(sym) & sym != ""
x <- x[keepg, ]; rownames(x) <- sym[keepg]
if (anyDuplicated(rownames(x))) {
  tot <- rowSums(x); keep_rn <- !duplicated(rownames(x))
  for (g in unique(rownames(x)[duplicated(rownames(x))])) {
    i <- which(rownames(x) == g); keep_rn[i[which.max(tot[i])]] <- TRUE
  }
  x <- x[keep_rn, , drop=FALSE]
}
mat <- log2(as.matrix(x) + 1)
pats <- c("1","2","3")
core_sign <- setNames(sign(p$logFC_meta), p$gene_symbol)

# per-patient paired logFC (T-C within patient), concordance with meta sign per patient
per_pat <- do.call(rbind, lapply(pats, function(i) {
  lfc_i <- mat[, paste0("T", i)] - mat[, paste0("C", i)]
  g <- intersect(p$gene_symbol, names(lfc_i))
  k <- sum(sign(lfc_i[g]) == core_sign[g], na.rm = TRUE)
  data.frame(patient = i, n_genes = length(g), concordant = k, prop = k / length(g))
}))
cat("\nper-patient direction concordance (core vs occlusion):\n"); print(per_pat)

# leave-one-patient-out cohort-level concordance (limma on 2 patients)
lopo <- do.call(rbind, lapply(pats, function(drop) {
  keep <- setdiff(pats, drop)
  cols <- c(paste0("C", keep), paste0("T", keep))
  pt <- factor(rep(keep, each = 2)); tm <- factor(rep(c("D0","D7"), length(keep)), levels = c("D0","D7"))
  des <- model.matrix(~ pt + tm)
  fit <- eBayes(lmFit(mat[, cols], des))
  lfc <- fit$coefficients[, "tmD7"]; names(lfc) <- rownames(fit)
  g <- intersect(p$gene_symbol, names(lfc))
  k <- sum(sign(lfc[g]) == core_sign[g], na.rm = TRUE)
  data.frame(held_out = drop, n_patients = 2, n_genes = length(g), concordant = k, prop = k / length(g))
}))
cat("\nleave-one-patient-out cohort concordance:\n"); print(lopo)
write.csv(per_pat, file.path(outdir,"occlusion_concordance_per_patient.csv"), row.names=FALSE)
write.csv(lopo, file.path(outdir,"occlusion_concordance_leave_one_out.csv"), row.names=FALSE)

## ---------- (3) cell-attribution robustness ----------
samples_d0 <- c("GSM9294393_1D0","GSM9294395_2D0","GSM9294397_3D0")
spot_list <- list()
for (f in samples_d0) {
  m <- readMM(gzfile(file.path(d, paste0(f, "_matrix.mtx.gz"))))
  feat <- read.delim(gzfile(file.path(d, paste0(f, "_features.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  bc <- read.delim(gzfile(file.path(d, paste0(f, "_barcodes.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  colnames(m) <- bc$V1
  lib <- Matrix::colSums(m); keep <- lib > 200
  Ln <- log2(t(t(m[,keep]) / lib[keep]) * 1e4 + 1)
  rownames(Ln) <- feat$V2
  spot_list[[f]] <- Ln
}
genes_common <- Reduce(intersect, lapply(spot_list, rownames))
Ln_all <- do.call(cbind, lapply(spot_list, function(L) L[genes_common, ]))
core_in <- intersect(p$gene_symbol, genes_common)

base_markers <- list(
  SMC        = c("ACTA2","MYH11","TAGLN","CNN1"),
  Endothelium= c("PECAM1","VWF","CDH5","CLDN5"),
  Fibroblast = c("DCN","LUM","COL1A2","GSN"),
  Immune     = c("PTPRC","LST1","S100A8","COTL1"),
  Adipocyte  = c("ADIPOQ","PLIN1","LPL")
)

assign_celltype <- function(markers, corr_thr = NULL) {
  sc <- sapply(names(markers), function(cn) {
    gg <- intersect(markers[[cn]], genes_common)
    if (length(gg) == 0) return(rep(NA_real_, ncol(Ln_all)))
    colMeans(Ln_all[gg, , drop=FALSE])
  })
  ct <- sapply(core_in, function(g) cor(Ln_all[g, ], sc[, "SMC"]))
  corr_mat <- sapply(names(markers), function(cn) {
    gg <- intersect(markers[[cn]], genes_common)
    if (length(gg) == 0) return(rep(NA_real_, length(core_in)))
    apply(Ln_all[core_in, ], 1, function(v) cor(v, sc[, cn]))
  })
  best <- names(markers)[max.col(corr_mat)]
  if (!is.null(corr_thr)) best[apply(corr_mat, 1, max) < corr_thr] <- "Unassigned"
  table(factor(best, levels = c(names(markers), "Unassigned")))
}

# (3a) leave-one-marker-out
lomo <- list()
for (cn in names(base_markers)) {
  for (mk in base_markers[[cn]]) {
    mk2 <- base_markers; mk2[[cn]] <- setdiff(base_markers[[cn]], mk)
    if (length(mk2[[cn]]) == 0) next
    tb <- assign_celltype(mk2)
    lomo[[paste(cn, mk, sep = "_drop")]] <- tb
  }
}
lomo_df <- as.data.frame(do.call(rbind, lomo))
lomo_summary <- data.frame(
  celltype = colnames(lomo_df),
  min = apply(lomo_df, 2, min), max = apply(lomo_df, 2, max))
cat("\nleave-one-marker-out attribution range (n variants =", nrow(lomo_df), "):\n"); print(lomo_summary)

# (3b) correlation-threshold sensitivity: unassigned if best |r| < thr
thr_tab <- do.call(rbind, lapply(c(0.2, 0.3, 0.4), function(th) {
  tb <- assign_celltype(base_markers, corr_thr = th)
  data.frame(threshold = th, as.list(tb))
}))
cat("\ncorrelation-threshold sensitivity:\n"); print(thr_tab)
write.csv(lomo_df, file.path(outdir,"attribution_leave_one_marker_out.csv"), row.names=FALSE)
write.csv(lomo_summary, file.path(outdir,"attribution_lomo_summary.csv"), row.names=FALSE)
write.csv(thr_tab, file.path(outdir,"attribution_threshold_sensitivity.csv"), row.names=FALSE)

## ---------- (4) paired-design completeness ----------
chk <- list()
c1 <- read.delim(gzfile(file.path(d,"GSE310013_raw_counts.txt.gz")), check.names=FALSE, nrows=1)
chk$GSE310013 <- colnames(c1)[-1]
c2 <- read.csv(gzfile(file.path(d,"GSE326710_Raw_counts_all_samples.csv.gz")), check.names=FALSE, nrows=1)
chk$GSE326710 <- colnames(c2)[-c(1,2)]
cat("\nGSE310013 samples:", paste(chk$GSE310013, collapse=", "), "\n")
cat("GSE326710 samples:", paste(chk$GSE326710, collapse=", "), "\n")

cat("\n=== 16_inference_audit.R done ===\n")
