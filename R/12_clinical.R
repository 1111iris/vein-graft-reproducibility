# 12_clinical.R — GSE241205 (occluded vs spare, n=3, FPKM): direction calibration (EXPLORATORY)
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
suppressPackageStartupMessages({library(limma)})
set.seed(20261029)
outdir <- "results"

sets <- list(
  SMC_CONTRACTILE = c("MYH11","ACTA2","TAGLN","CNN1","MYL9","MYLK","LMOD1","RGS5"),
  ECM_REMODELING  = c("MMP9","TIMP1","COL1A1","POSTN","MMP2","COL3A1","FN1","LOX"),
  INFLAMMATION    = c("IL6","CXCL2","CCL2","VCAM1","ICAM1","CXCL8","CCL20"),
  MS_CHANNEL      = c("PIEZO1","PIEZO2","TMC4","TMC5","TMC6","TMC8","TMEM63A","TMEM63B"),
  OSTEOGENIC      = c("RUNX2","SOX9","SPP1","BGLAP","BMP2")
)

x <- read.delim(gzfile("data/GSE241205_gene_fpkm.txt.gz"), check.names=FALSE)
rownames(x) <- x$Symbol; x$Symbol <- NULL  # values are ENSG IDs (metadata sloppiness)
# map ENSG -> symbol using GSE326710 annotation
x2 <- read.csv(gzfile("data/GSE326710_Raw_counts_all_samples.csv.gz"), check.names=FALSE)
map <- x2[["Gene Symbol"]]; names(map) <- x2[["Gene ID"]]
sym <- unname(map[rownames(x)])
keepg <- !is.na(sym) & sym != ""
x <- x[keepg, ]; rownames(x) <- sym[keepg]
# collapse duplicates: retain the duplicate-symbol row with the highest
# summed FPKM across samples (a convention choice, not a uniqueness claim)
if (anyDuplicated(rownames(x))) {
  tot <- rowSums(x); keep_rn <- !duplicated(rownames(x))
  for (g in unique(rownames(x)[duplicated(rownames(x))])) {
    i <- which(rownames(x) == g); keep_rn[i[which.max(tot[i])]] <- TRUE
  }
  x <- x[keep_rn, , drop=FALSE]
}
time <- factor(ifelse(grepl("^T", colnames(x)), "D7", "D0"), levels=c("D0","D7"))  # T=occluded
patient <- factor(sub("^[CT]", "", colnames(x)))
mat <- log2(as.matrix(x) + 1)
design <- model.matrix(~ patient + time)
fit <- eBayes(lmFit(mat, design))
tt <- topTable(fit, coef="timeD7", number=Inf, sort.by="none")
tt$SE <- tt$logFC / tt$t; tt$gene <- rownames(tt)
write.csv(tt, file.path(outdir,"de_GSE241205_fpkm_full.csv"), row.names=FALSE)

# signature-level direction scores vs model meta signature
meta <- read.csv(file.path(outdir,"meta_analysis.csv"))
ranks <- setNames(meta$logFC_meta, meta$symbol)
sig_scores <- sapply(names(sets), function(nm) {
  g <- intersect(sets[[nm]], rownames(mat))
  if (length(g) < 2) return(NA)
  mean(tt[g,"logFC"])
})
# gene-set concordance: rank-based (fgsea) using meta signature ranks vs occlusion ranks
suppressPackageStartupMessages(library(fgsea))
ranks_occ <- setNames(tt$logFC, tt$gene)
fg <- fgsea(sets[sapply(sets, function(s) sum(s %in% names(ranks_occ))) > 1], stats=ranks_occ, eps=1e-50)
fg_out <- data.frame(set=fg$pathway, NES_occlusion=fg$NES, pval=fg$pval, padj=fg$padj, size=fg$size)
cat("=== GSE241205 signature mean logFC (occluded vs spare; positive = up in occlusion) ===\n")
print(round(sig_scores, 3))
cat("\n=== fgsea of locked sets on occlusion ranks ===\n"); print(fg_out)
write.csv(data.frame(set=names(sets), mean_logFC_occlusion=sig_scores),
          file.path(outdir,"gse241205_signature_direction.csv"), row.names=FALSE)
write.csv(fg_out, file.path(outdir,"fgsea_gse241205.csv"), row.names=FALSE)

# cross-signature concordance: model meta NES direction vs occlusion NES direction
fg_model <- read.csv(file.path(outdir,"fgsea_formal.csv"))
fm <- fg_model[fg_model$dataset=="GSE310013", ]
cmp <- merge(fm[,c("set","NES")], fg_out[,c("set","NES_occlusion")], by="set")
cmp$concordant <- sign(cmp$NES) == sign(cmp$NES_occlusion)
cat("\n=== Model (GSE310013 NES) vs occlusion NES ===\n"); print(cmp)
write.csv(cmp, file.path(outdir,"model_vs_occlusion_concordance.csv"), row.names=FALSE)

# PIEZO1 one-shot
pz <- tt["PIEZO1", c("logFC","SE","P.Value")]
cat("\nPIEZO1 occluded vs spare: logFC", round(pz$logFC,3), "P", signif(pz$P.Value,3), "\n")
write.csv(data.frame(dataset="GSE241205", logFC=pz$logFC, SE=pz$SE, P=pz$P.Value),
          file.path(outdir,"piezo1_gse241205.csv"), row.names=FALSE)
