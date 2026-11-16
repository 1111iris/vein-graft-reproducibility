# Round 9: stratified occlusion-direction calibration (descriptive).
# Run from the project root (the directory containing data/ and results/):
#   Rscript 24_stratified_occlusion_concordance.R
# Reads results/robust_core_priority_all.csv and results/meta_analysis.csv; writes
# results/stratified_occlusion_concordance.csv. No new statistical tests: all
# strata are descriptive concordance fractions (n = 3 calibration patients).

set.seed(20261029)

core <- read.csv("results/robust_core_priority_all.csv", stringsAsFactors = FALSE)
meta <- read.csv("results/meta_analysis.csv", stringsAsFactors = FALSE)

# occlusion logFC per gene (from the archived per-gene DE table)
occ <- read.csv("results/de_GSE241205_fpkm_full.csv", stringsAsFactors = FALSE)
occ_dir <- setNames(sign(occ$logFC), occ$gene)

# locked gene sets (Supplementary Material S1.1)
gene_sets <- list(
  SMC_CONTRACTILE = c("MYH11","ACTA2","TAGLN","CNN1","MYL9","MYLK","LMOD1","RGS5"),
  ECM_REMODELING  = c("MMP9","TIMP1","COL1A1","POSTN","MMP2","COL3A1","FN1","LOX"),
  INFLAMMATION    = c("IL6","CXCL2","CCL2","VCAM1","ICAM1","CXCL8","CCL20"),
  MS_CHANNEL      = c("PIEZO1","PIEZO2","TMC4","TMC5","TMC6","TMC8","TMEM63A","TMEM63B"),
  OSTEOGENIC      = c("RUNX2","SOX9","SPP1","BGLAP","BMP2")
)

# per-gene occlusion direction for all commonly detected genes
meta$occ_concordant_all <- ifelse(
  is.na(occ_dir[meta$gene]) | meta$logFC_meta == 0, NA,
  as.character(sign(meta$logFC_meta)) == as.character(occ_dir[meta$gene])
)

# core-gene concordance (recompute from occ_dir to stay self-contained)
core$occ_concordant2 <- ifelse(
  is.na(occ_dir[core$gene]) | core$logFC_meta == 0, NA,
  as.character(sign(core$logFC_meta)) == as.character(occ_dir[core$gene])
)
stopifnot(all(core$occ_concordant2 == core$occ_concordant, na.rm = TRUE))

frac <- function(x) {
  x <- x[!is.na(x)]
  c(n = length(x), concordant = sum(x), rate = round(100 * mean(x), 1))
}

rows <- list()
add_row <- function(stratum_type, stratum, cohort, stats) {
  rows[[length(rows) + 1]] <<- data.frame(
    stratum_type, stratum, cohort,
    n = stats["n"], concordant = stats["concordant"], rate_pct = stats["rate"]
  )
}

# 1. overall
add_row("overall", "all_commonly_detected_genes", "all_genes", frac(meta$occ_concordant_all))
add_row("overall", "robust_core", "core", frac(core$occ_concordant2))

# 2. by locked program (core members only)
for (s in names(gene_sets)) {
  idx <- core$gene_symbol %in% gene_sets[[s]]
  add_row("program", s, "core", frac(core$occ_concordant2[idx]))
}

# 3. by |logFC_meta| tertile, core vs same-band background
br <- quantile(abs(core$logFC_meta), probs = c(1/3, 2/3))
core$tertile <- cut(abs(core$logFC_meta), breaks = c(-Inf, br, Inf), labels = c("T1_low","T2_mid","T3_high"))
meta$tertile <- cut(abs(meta$logFC_meta), breaks = c(-Inf, br, Inf), labels = c("T1_low","T2_mid","T3_high"))
for (l in levels(core$tertile)) {
  add_row("effect_size_tertile", as.character(l), "core", frac(core$occ_concordant2[core$tertile == l]))
  add_row("effect_size_tertile", as.character(l), "all_genes_same_band", frac(meta$occ_concordant_all[meta$tertile == l]))
}

# 4. by spatial marker-associated category (forced assignment)
for (ct in sort(unique(core$best_celltype))) {
  idx <- core$best_celltype == ct
  add_row("cell_category", ct, "core", frac(core$occ_concordant2[idx]))
}

out <- do.call(rbind, rows)
write.csv(out, "results/stratified_occlusion_concordance.csv", row.names = FALSE)
cat("written: results/stratified_occlusion_concordance.csv\n")
print(out, row.names = FALSE)

# supplementary Table S10
tbl <- out
dir.create("results/supplementary_tables", showWarnings = FALSE, recursive = TRUE)
write.csv(tbl, "results/supplementary_tables/TableS10_stratified_occlusion_concordance.csv", row.names = FALSE)
cat("written: results/supplementary_tables/TableS10_stratified_occlusion_concordance.csv\n")
