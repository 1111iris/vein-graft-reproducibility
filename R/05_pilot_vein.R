# 05_pilot_vein.R — Phase-2 pilot: PIEZO1 detectability & predefined marker panels
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# EXPLORATORY ONLY: small-n pilot to inform feasibility audit, not confirmatory analysis.
# All gene panels below are standard markers fixed BEFORE looking at results
# (same MS_CHANNEL set as phase-1 audit for continuity).

suppressPackageStartupMessages({library(limma)})
set.seed(20261029)
d <- "data"

panels <- list(
  SMC_contractile = c("MYH11","ACTA2","TAGLN","CNN1","MYL9"),
  ECM_remodeling  = c("MMP9","TIMP1","COL1A1","POSTN","MMP2"),
  Inflammation    = c("IL6","CXCL2","CCL2","VCAM1","ICAM1"),
  MS_CHANNEL      = c("PIEZO1","PIEZO2","TMC4","TMC5","TMC6","TMC8","TMEM63A","TMEM63B")
)
genes_wanted <- unique(unlist(panels))

run_limma_paired <- function(mat, patient, time, label, logged=FALSE) {
  # mat: genes x samples; patient/time: vectors
  if (!logged) mat <- log2(mat + 1)
  design <- model.matrix(~ patient + time)
  fit <- eBayes(lmFit(mat, design))
  tt <- topTable(fit, coef="timeD7", number=Inf, sort.by="none")
  data.frame(dataset=label, gene=rownames(tt), logFC=tt$logFC, t=tt$t, P=tt$P.Value,
             stringsAsFactors=FALSE)
}

res <- list()

# ---- GSE310013: raw counts, 6 patients (R1-R3,S1-S3) paired D0/D7 ----
x <- read.delim(gzfile(file.path(d,"GSE310013_raw_counts.txt.gz")), check.names=FALSE)
rownames(x) <- x$Geneid; x$Geneid <- NULL
time <- factor(ifelse(grepl("D7", colnames(x)), "D7", "D0"), levels=c("D0","D7"))
patient <- sub("D[07]_001$", "", colnames(x))
stopifnot(length(unique(patient))==6, all(table(patient)==2))
g <- intersect(genes_wanted, rownames(x))
res$GSE310013 <- run_limma_paired(x[g,], patient, time, "GSE310013_bulk_n6")

# ---- GSE326710: raw counts, 4 patients (V1-V4) paired D0/D7 ----
x <- read.csv(gzfile(file.path(d,"GSE326710_Raw_counts_all_samples.csv.gz")), check.names=FALSE)
rn <- x[["Gene Symbol"]]; x[["Gene Symbol"]] <- NULL; x[["Type"]] <- NULL; x[["Gene ID"]] <- NULL
x <- as.matrix(x); rownames(x) <- rn
time <- factor(ifelse(grepl("D7", colnames(x)), "D7", "D0"), levels=c("D0","D7"))
patient <- sub("_D[07] Read Count$", "", colnames(x))
stopifnot(length(unique(patient))==4, all(table(patient)==2))
g <- intersect(genes_wanted, rownames(x))
res$GSE326710 <- run_limma_paired(x[g,], patient, time, "GSE326710_bulk_n4")

# ---- GSE241205: FPKM only, 3 patients paired (C=spare, T=occluded) ----
x <- read.delim(gzfile(file.path(d,"GSE241205_gene_fpkm.txt.gz")), check.names=FALSE)
rownames(x) <- x$Symbol; x$Symbol <- NULL  # header says Symbol, values are ENSG (metadata sloppiness)
# map ENSG -> symbol via the other datasets' gene info is not available here; use org.Hb annotation
# Simple approach: use biomaRt-free mapping from GSE310013 Geneid? Not available. Use limma on ENSG,
# then report PIEZO1 ENSG ENSG00000177556 directly.
time <- factor(ifelse(grepl("^T", colnames(x)), "D7", "D0"), levels=c("D0","D7"))  # T=occluded ~ 'late'
patient <- sub("^[CT]", "", colnames(x))
stopifnot(length(unique(patient))==3, all(table(patient)==2))
g <- intersect(c(genes_wanted, "ENSG00000177556"), rownames(x))
r <- run_limma_paired(x[g,], patient, time, "GSE241205_fpkm_n3")
r$gene[r$gene=="ENSG00000177556"] <- "PIEZO1(ENSG)"
res$GSE241205 <- r

# ---- GSE299685: GeoMx Q3-normalized AOIs, 4 patients, nested structure ----
library(readxl)
xl <- excel_sheets(file.path(d,"GSE299685_Normalized_Counts.xlsx"))
cat("GSE299685 xlsx sheets:", paste(xl, collapse=", "), "\n")
x <- as.data.frame(read_excel(file.path(d,"GSE299685_Normalized_Counts.xlsx"), sheet=xl[1]))
cat("dim:", dim(x), "; first cols:", paste(head(colnames(x),8), collapse=" | "), "\n")
# find gene id col and AOI annotation rows (GeoMx exports have annotation block)
print(head(x[,1:6]))

out <- do.call(rbind, res)
write.csv(out, "pilot_results_bulk.csv", row.names=FALSE)
cat("\n=== Pilot results (bulk, paired limma) ===\n")
print(subset(out, gene %in% c("PIEZO1","PIEZO1(ENSG)","MYH11","ACTA2","TAGLN","MMP9","IL6","CXCL2")))
