# 14_robust_core.R — Robust core signature: cell/spatial attribution + priority ranking
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Attribution sources: (a) Visium spot-level correlation with cell-type marker scores,
# (b) GeoMx MYH11pos/neg segment enrichment (media vs intima/adventitia bias).
suppressPackageStartupMessages({library(limma); library(Matrix); library(readxl)})
set.seed(20261029)
outdir <- "results"; d <- "data"

meta <- read.csv(file.path(outdir,"meta_analysis.csv"))
fdr1 <- read.csv(file.path(outdir,"de_GSE310013_full.csv"))
fdr2 <- read.csv(file.path(outdir,"de_GSE326710_full.csv"))
meta$FDR_310013 <- fdr1$adj.P.Val[match(meta$gene, fdr1$gene)]
meta$FDR_326710 <- fdr2$adj.P.Val[match(meta$gene, fdr2$gene)]
core <- subset(meta, FDR_310013 < 0.05 & FDR_326710 < 0.05 &
                 sign(logFC_310013) == sign(logFC_326710))
cat("robust core genes:", nrow(core), "\n")

# ---- (a) Visium spot-level cell attribution (pooled D0 sections, all 3 patients) ----
cell_markers <- list(
  SMC       = c("ACTA2","MYH11","TAGLN","CNN1"),
  Endothelium=c("PECAM1","VWF","CDH5","CLDN5"),
  Fibroblast= c("DCN","LUM","COL1A2","GSN"),
  Immune    = c("PTPRC","LST1","S100A8","COTL1"),
  Adipocyte = c("ADIPOQ","PLIN1","LPL")
)
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
cell_sc <- sapply(cell_markers, function(g) {
  gg <- intersect(g, genes_common)
  colMeans(Ln_all[gg, , drop=FALSE])
})
core_in <- intersect(core$symbol, genes_common)
cat("core genes with valid Visium assignments:", length(core_in), "\n")
# correlation of each core gene with each cell score across spots
ct_corr <- sapply(names(cell_markers), function(cn) {
  apply(Ln_all[core_in, ], 1, function(v) cor(v, cell_sc[, cn]))
})
ct_corr <- as.data.frame(ct_corr)
ct_corr$gene <- core_in
ct_best <- names(cell_markers)[max.col(as.matrix(ct_corr[, names(cell_markers)]))]
ct_corr$best_celltype <- ct_best

# ---- (a1) attribution denominator audit (explicit and asserted) ------------
# Manuscript denominators: 4,148 core genes; 4,110 represented in the Visium
# gene matrix; 91 of those lack valid correlations (constant expression), so
# 4,019 genes have valid original assignments. Asserted so the denominators
# cannot drift silently if gene mappings or data versions change.
den <- data.frame(
  step = c("robust_core", "in_visium_matrix", "not_in_visium_matrix",
           "in_matrix_no_valid_correlation", "valid_assignment"),
  n = c(nrow(core), length(core_in), nrow(core) - length(core_in),
        sum(is.na(ct_corr$SMC)), sum(!is.na(ct_corr$SMC))))
stopifnot(den$n[den$step == "robust_core"] == 4148,
          den$n[den$step == "in_visium_matrix"] == 4110,
          den$n[den$step == "not_in_visium_matrix"] == 38,
          den$n[den$step == "in_matrix_no_valid_correlation"] == 91,
          den$n[den$step == "valid_assignment"] == 4019)
write.csv(den, file.path(outdir, "attribution_denominator_audit.csv"), row.names = FALSE)
cat("attribution denominator audit: 4148 = 4110 in matrix (+38 not detected);",
    "4110 = 4019 valid assignments (+91 no valid correlation)\n")

# ---- (b) GeoMx MYH11pos/neg enrichment (media vs non-media bias) ----
gdir <- "geo_pages/gsm299685"
ann <- do.call(rbind, lapply(list.files(gdir, full.names=TRUE), function(f) {
  L <- readLines(f)
  t <- sub("^!Sample_title = ", "", L[grepl("^!Sample_title", L)][1])
  s <- sub("^!Sample_supplementary_file_1 = .*/", "", L[grepl("^!Sample_supplementary_file_1", L)][1])
  p <- strsplit(t, "_")[[1]]
  data.frame(dcc=sub("\\.gz$","",sub("^GSM\\d+_","",s)), patient=p[1], time=sub("ROI-\\d+_","",p[3]), seg=p[4])
}))
gx <- as.data.frame(read_excel(file.path(d,"GSE299685_Normalized_Counts.xlsx")))
rownames(gx) <- gx$gene; gx$gene <- NULL
ann <- ann[ann$dcc %in% colnames(gx), ]
gm <- as.matrix(gx[, ann$dcc])
ann$seg <- relevel(factor(ann$seg), ref="MYH11pos")
design_g <- model.matrix(~ patient + time + seg, data=ann)
fit_g <- eBayes(lmFit(gm, design_g))
tt_g <- topTable(fit_g, coef="segMYH11neg", number=Inf, sort.by="none")
seg_lfc <- setNames(tt_g$logFC, rownames(tt_g))  # positive = enriched in MYH11neg (intima/adventitia-biased)

# ---- (c) occlusion direction (GSE241205) ----
occ <- read.csv(file.path(outdir,"de_GSE241205_fpkm_full.csv"))
occ_lfc <- setNames(occ$logFC, occ$gene)

# ---- assemble priority table ----
p <- core
colnames(p)[colnames(p)=="symbol"] <- "gene_symbol"
p$ct_SMC <- ct_corr$SMC[match(p$gene_symbol, ct_corr$gene)]
p$ct_Endo <- ct_corr$Endothelium[match(p$gene_symbol, ct_corr$gene)]
p$ct_Fibro <- ct_corr$Fibroblast[match(p$gene_symbol, ct_corr$gene)]
p$ct_Immune <- ct_corr$Immune[match(p$gene_symbol, ct_corr$gene)]
p$best_celltype <- ct_best[match(p$gene_symbol, ct_corr$gene)]
p$geomx_seg_lfc <- seg_lfc[p$gene_symbol]
p$seg_region <- with(p, ifelse(is.na(geomx_seg_lfc), NA,
                        ifelse(geomx_seg_lfc > 0.5, "MYH11neg-biased",
                        ifelse(geomx_seg_lfc < -0.5, "MYH11pos-biased (media)", "diffuse"))))
p$occ_logFC <- occ_lfc[p$gene_symbol]
p$occ_concordant <- sign(p$logFC_meta) == sign(p$occ_logFC)
# priority score: effect size + attribution clarity + occlusion concordance
p$priority_score <- with(p, abs(logFC_meta) * (-log10(pmax(FDR_meta, 1e-300))) * ifelse(occ_concordant, 1, 0.5))
p <- p[order(-p$priority_score), ]
write.csv(p, file.path(outdir,"robust_core_priority_all.csv"), row.names=FALSE)

top30 <- p[1:30, c("gene_symbol","logFC_310013","logFC_326710","logFC_meta","FDR_meta",
                   "best_celltype","seg_region","occ_logFC","occ_concordant","priority_score")]
cat("\n=== TOP 30 robust core genes ===\n"); print(top30, digits=3)
write.csv(top30, file.path(outdir,"robust_core_top30.csv"), row.names=FALSE)
cat("\nCell-type composition of full core:\n"); print(table(p$best_celltype))
cat("Region composition:\n"); print(table(p$seg_region))
cat("Occlusion concordance rate:", mean(p$occ_concordant, na.rm=TRUE), "\n")
saveRDS(list(core=p, top30=top30), file.path(outdir,"robust_core_state.rds"))
