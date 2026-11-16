# 06_pilot_geomx.R — GSE299685 GeoMx AOI-level pilot (EXPLORATORY)
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# 36 AOIs nested in 4 patients; annotation rebuilt from GSM titles (GEO pages).
suppressPackageStartupMessages({library(limma); library(readxl)})
set.seed(20261029)

# --- annotation from GSM pages ---
gdir <- "geo_pages/gsm299685"
files <- list.files(gdir, full.names=TRUE)
ann <- do.call(rbind, lapply(files, function(f) {
  L <- readLines(f)
  t <- sub("^!Sample_title = ", "", L[grepl("^!Sample_title", L)][1])
  s <- sub("^!Sample_supplementary_file_1 = .*/", "", L[grepl("^!Sample_supplementary_file_1", L)][1])
  parts <- strsplit(t, "_")[[1]]
  data.frame(title=t, dcc=sub("\\.gz$","",s), patient=parts[1],
             roi=parts[2], time=sub("ROI-\\d+_","",parts[3]), seg=parts[4], stringsAsFactors=FALSE)
}))
ann$time[ann$time=="D0"] <- "D0"; ann$time[ann$time=="D7"] <- "D7"
cat("AOIs:", nrow(ann), "; patients:", paste(unique(ann$patient), collapse=","), "\n")
print(table(ann$patient, ann$time))

# --- expression matrix ---
x <- as.data.frame(read_excel("data/GSE299685_Normalized_Counts.xlsx"))
rownames(x) <- x$gene; x$gene <- NULL
m <- as.matrix(x)
ann$dcc <- sub("^GSM\\d+_", "", ann$dcc)  # xlsx columns lack GSM prefix
cat("matrix:", dim(m), "; all ann dcc in cols:", all(ann$dcc %in% colnames(m)), "\n")
missing <- setdiff(ann$dcc, colnames(m))
cat("DCC files absent from normalized matrix (QC-failed AOIs):", paste(missing, collapse=","), "\n")
ann <- ann[ann$dcc %in% colnames(m), ]
m <- m[, ann$dcc]

panels <- list(
  SMC_contractile = c("MYH11","ACTA2","TAGLN","CNN1","MYL9"),
  ECM_remodeling  = c("MMP9","TIMP1","COL1A1","POSTN","MMP2"),
  Inflammation    = c("IL6","CXCL2","CCL2","VCAM1","ICAM1"),
  MS_CHANNEL      = c("PIEZO1","PIEZO2","TMC4","TMC5","TMC6","TMC8","TMEM63A","TMEM63B")
)
genes_wanted <- unique(unlist(panels))
g <- intersect(genes_wanted, rownames(m))
cat("panel genes present:", length(g), "/", length(genes_wanted),
    "; missing:", paste(setdiff(genes_wanted, g), collapse=","), "\n")

mat <- log2(m[g, ] + 1)
ann$time <- factor(ann$time, levels=c("D0","D7"))
ann$seg <- factor(ann$seg)
design <- model.matrix(~ patient + seg + time, data=ann)
fit <- eBayes(lmFit(mat, design))
tt <- topTable(fit, coef="timeD7", number=Inf, sort.by="none")
out <- data.frame(dataset="GSE299685_geomx_36AOI_4pt", gene=rownames(tt),
                  logFC=tt$logFC, t=tt$t, P=tt$P.Value)
write.csv(out, "pilot_results_geomx.csv", row.names=FALSE)
cat("\n=== GeoMx pilot (time D7 vs D0, patient-blocked) ===\n")
print(out)
# PIEZO1 by segment (exploratory interaction)
design2 <- model.matrix(~ patient + time + seg + time:seg, data=ann)
fit2 <- eBayes(lmFit(mat, design2))
print(topTable(fit2, coef=grep("timeD7:seg", colnames(design2)), number=Inf, sort.by="none")["PIEZO1",])
