# 25_leave_one_gene_out_attribution.R — advisor round 11, must-fix item 1
# Removes self-marker circularity from spatial cell attribution:
# for each core gene, its own expression is excluded from any marker score that
# contains it, before the gene-vs-score correlation is computed (leave-one-gene-out, LOGO).
# Also reports which core genes are themselves marker-set members and how the
# forced-assignment composition changes relative to the original (script 14) attribution.
# Run from the project root: Rscript 25_leave_one_gene_out_attribution.R
suppressPackageStartupMessages({library(Matrix)})
set.seed(20261029)
outdir <- "results"; d <- "data"

p <- read.csv(file.path(outdir, "robust_core_priority_all.csv"))

cell_markers <- list(
  SMC        = c("ACTA2","MYH11","TAGLN","CNN1"),
  Endothelium= c("PECAM1","VWF","CDH5","CLDN5"),
  Fibroblast = c("DCN","LUM","COL1A2","GSN"),
  Immune     = c("PTPRC","LST1","S100A8","COTL1"),
  Adipocyte  = c("ADIPOQ","PLIN1","LPL")
)

samples_d0 <- c("GSM9294393_1D0","GSM9294395_2D0","GSM9294397_3D0")
spot_list <- list()
for (f in samples_d0) {
  m <- readMM(gzfile(file.path(d, paste0(f, "_matrix.mtx.gz"))))
  feat <- read.delim(gzfile(file.path(d, paste0(f, "_features.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  bc <- read.delim(gzfile(file.path(d, paste0(f, "_barcodes.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  colnames(m) <- bc$V1
  lib <- Matrix::colSums(m); keep <- lib > 200
  Ln <- log2(t(t(m[, keep]) / lib[keep]) * 1e4 + 1)
  rownames(Ln) <- feat$V2
  spot_list[[f]] <- Ln
}
genes_common <- Reduce(intersect, lapply(spot_list, rownames))
Ln_all <- do.call(cbind, lapply(spot_list, function(L) L[genes_common, ]))
core_in <- intersect(p$gene_symbol, genes_common)
# same denominator as the manuscript: genes with a valid original attribution
# (excludes 38 genes absent from Visium and 91 genes with degenerate spatial variance)
orig_all <- setNames(p$best_celltype, p$gene_symbol)[core_in]
core_in <- core_in[!is.na(orig_all)]
orig <- orig_all[core_in]
cat("core genes with valid original attribution:", length(core_in), "\n")

# marker scores as in script 14 (full sets)
cell_sc <- sapply(cell_markers, function(g) {
  gg <- intersect(g, genes_common)
  colMeans(Ln_all[gg, , drop=FALSE])
})

# ---- (a) which core genes are themselves marker-set members? ----
member_tab <- do.call(rbind, lapply(names(cell_markers), function(cn) {
  m <- intersect(cell_markers[[cn]], core_in)
  if (length(m) == 0) return(NULL)
  data.frame(celltype = cn, marker_gene = m, stringsAsFactors = FALSE)
}))
cat("\ncore genes that are themselves marker-set members:\n"); print(member_tab)

# ---- (b) leave-one-gene-out correlations ----
# for gene g and cell type cn: if g is in marker set cn, recompute the score
# across spots excluding g; otherwise the score is unchanged.
logo_corr <- t(vapply(core_in, function(g) {
  vapply(names(cell_markers), function(cn) {
    gg <- intersect(cell_markers[[cn]], genes_common)
    gg <- setdiff(gg, g)          # exclude the gene itself when it is a marker member
    if (length(gg) == 0) return(NA_real_)
    cor(Ln_all[g, ], colMeans(Ln_all[gg, , drop=FALSE]))
  }, numeric(1))
}, numeric(length(cell_markers))))
colnames(logo_corr) <- names(cell_markers)

best_logo <- colnames(logo_corr)[max.col(logo_corr)]
best_logo[is.na(rowSums(logo_corr))] <- NA   # degenerate genes stay unassigned, as in script 14
best_logo_maxabs <- apply(abs(logo_corr), 1, max, na.rm = TRUE)  # largest absolute correlation of any score (diagnostic column logo_max_abs_r)
names(best_logo) <- core_in

# assignment rule check: highest signed correlation (primary, used throughout)
# versus largest absolute correlation regardless of sign
best_abs <- colnames(logo_corr)[max.col(abs(logo_corr))]
best_abs[is.na(rowSums(logo_corr))] <- NA
names(best_abs) <- core_in
rule_diff <- sum(!is.na(best_logo) & !is.na(best_abs) & best_logo != best_abs)
cat(sprintf("\nassignment-rule check: %d/%d genes change label under a max-|r| rule versus the primary max-signed-r rule\n",
            rule_diff, sum(!is.na(best_logo))))
comp_abs <- table(factor(best_abs, levels = c(names(cell_markers), "Unassigned")))
cat("composition under max-|r| rule:\n"); print(comp_abs)

agree <- mean(best_logo == orig, na.rm = TRUE)
cat(sprintf("\nLOGO vs original assignment agreement: %d/%d = %.1f%%\n",
            sum(best_logo == orig, na.rm = TRUE), sum(!is.na(best_logo) & !is.na(orig)), 100 * agree))

comp <- rbind(
  original = table(factor(orig, levels = c(names(cell_markers), "Unassigned"))),
  logo     = table(factor(best_logo, levels = c(names(cell_markers), "Unassigned"))))
pct <- round(100 * prop.table(comp, margin = 1), 1)
cat("\ncomposition (counts):\n"); print(comp)
cat("\ncomposition (% of 4,019):\n"); print(pct)

# ---- (c) marker-member genes: original vs LOGO classification ----
mm <- member_tab$marker_gene
mm_cmp <- data.frame(
  gene = mm,
  marker_set = member_tab$celltype,
  original = unname(orig[mm]),
  logo = unname(best_logo[mm]),
  stringsAsFactors = FALSE)
cat("\nmarker-member core genes, original vs LOGO classification:\n"); print(mm_cmp)

# ---- (d) per-gene output ----
out <- data.frame(
  gene_symbol = core_in,
  best_celltype_original = unname(orig),
  best_celltype_logo = unname(best_logo),
  unchanged = unname(orig == best_logo),
  logo_max_abs_r = unname(best_logo_maxabs),
  as.data.frame(logo_corr),
  row.names = NULL, check.names = FALSE)
write.csv(out, file.path(outdir, "attribution_leave_one_gene_out.csv"), row.names = FALSE)
write.csv(mm_cmp, file.path(outdir, "attribution_marker_member_genes.csv"), row.names = FALSE)

# composition table: celltype rows x rule columns (labels fixed, round-12 review)
ct_names <- colnames(comp)
out_comp <- data.frame(
  celltype = ct_names,
  original_n = as.integer(comp["original", ct_names]),
  logo_n = as.integer(comp["logo", ct_names]),
  original_pct = as.numeric(pct["original", ct_names]),
  logo_pct = as.numeric(pct["logo", ct_names]),
  stringsAsFactors = FALSE)
write.csv(out_comp, file.path(outdir, "attribution_logo_composition.csv"), row.names = FALSE)
write.csv(data.frame(rule = c("max_signed_r (primary)", "max_abs_r"),
                     n_changed_versus_primary = c(0, rule_diff)),
          file.path(outdir, "attribution_rule_comparison.csv"), row.names = FALSE)

cat("\n=== 25_leave_one_gene_out_attribution.R done ===\n")
