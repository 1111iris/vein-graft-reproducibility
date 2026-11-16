# 11_spatial.R — Visium (GSE310351) region-level analysis + GeoMx formal gene sets
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Locked rules: PREREGISTRATION_gene_sets.md. Spots aggregated to patient level (paired).
suppressPackageStartupMessages({library(limma); library(Matrix); library(readxl); library(ggplot2)})
set.seed(20261029)
outdir <- "results"; d <- "data"

sets <- list(
  SMC_CONTRACTILE = c("MYH11","ACTA2","TAGLN","CNN1","MYL9","MYLK","LMOD1","RGS5"),
  ECM_REMODELING  = c("MMP9","TIMP1","COL1A1","POSTN","MMP2","COL3A1","FN1","LOX"),
  INFLAMMATION    = c("IL6","CXCL2","CCL2","VCAM1","ICAM1","CXCL8","CCL20"),
  MS_CHANNEL      = c("PIEZO1","PIEZO2","TMC4","TMC5","TMC6","TMC8","TMEM63A","TMEM63B"),
  OSTEOGENIC      = c("RUNX2","SOX9","SPP1","BGLAP","BMP2")
)
samples <- data.frame(
  file=c("GSM9294393_1D0","GSM9294394_1D7","GSM9294395_2D0","GSM9294396_2D7","GSM9294397_3D0","GSM9294398_3D7"),
  patient=paste0("P", c(1,1,2,2,3,3)), time=c("D0","D7","D0","D7","D0","D7"), stringsAsFactors=FALSE)

score_sample <- function(f, sets) {
  m <- readMM(gzfile(file.path(d, paste0(f, "_matrix.mtx.gz"))))
  feat <- read.delim(gzfile(file.path(d, paste0(f, "_features.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  bc <- read.delim(gzfile(file.path(d, paste0(f, "_barcodes.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  pos <- read.csv(gzfile(file.path(d, paste0(f, "_tissue_positions.csv.gz"))), header=FALSE)
  if (ncol(pos) >= 8) colnames(pos) <- c("barcode","in_tissue","row","col","imagerow","imagecol","pxrow","pycol")
  else colnames(pos) <- c("barcode","in_tissue","row","col","imagerow","imagecol")
  colnames(m) <- bc$V1; m <- m[, bc$V1]
  lib <- Matrix::colSums(m)
  keep <- lib > 200
  m <- m[, keep]; lib <- lib[keep]; pos <- pos[match(colnames(m), pos$barcode), ]
  # log-normalize per spot
  Ln <- log2(t(t(m) / lib) * 1e4 + 1)
  sc <- sapply(names(sets), function(nm) {
    g <- intersect(sets[[nm]], feat$V2)
    if (length(g) < 2) return(rep(NA, ncol(Ln)))
    colMeans(Ln[match(g, feat$V2), , drop=FALSE])
  })
  list(scores=sc, pos=pos, n_spots=sum(keep), lib=lib)
}

res <- list(); spot_all <- list()
for (i in seq_len(nrow(samples))) {
  s <- samples[i, ]
  r <- score_sample(s$file, sets)
  res[[s$file]] <- data.frame(file=s$file, patient=s.patient <- s$patient, time=s$time,
                              n_spots=r$n_spots, med_lib=median(r$lib),
                              SMC=mean(r$scores[,"SMC_CONTRACTILE"], na.rm=TRUE),
                              ECM=mean(r$scores[,"ECM_REMODELING"], na.rm=TRUE),
                              INFL=mean(r$scores[,"INFLAMMATION"], na.rm=TRUE),
                              MS=mean(r$scores[,"MS_CHANNEL"], na.rm=TRUE),
                              OSTEO=mean(r$scores[,"OSTEOGENIC"], na.rm=TRUE),
                              PIEZO1=NA)
  # per-spot table for plots (PIEZO1 + SMC + INFL)
  df <- data.frame(file=s$file, patient=s$patient, time=s$time,
                   px=r$pos$row, py=r$pos$col,
                   SMC=as.numeric(r$scores[,"SMC_CONTRACTILE"]),
                   INFL=as.numeric(r$scores[,"INFLAMMATION"]))
  feat <- read.delim(gzfile(file.path(d, paste0(s$file, "_features.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  m <- readMM(gzfile(file.path(d, paste0(s$file, "_matrix.mtx.gz"))))
  bc <- read.delim(gzfile(file.path(d, paste0(s$file, "_barcodes.tsv.gz"))), header=FALSE, stringsAsFactors=FALSE)
  colnames(m) <- bc$V1; m <- m[, bc$V1]; lib <- Matrix::colSums(m); keep <- lib > 200
  Ln <- log2(t(t(m[,keep]) / lib[keep]) * 1e4 + 1)
  gi <- match("PIEZO1", feat$V2)
  df$PIEZO1 <- as.numeric(Ln[gi, ])
  spot_all[[s$file]] <- df
}
patient_summary <- do.call(rbind, res)
cat("=== Visium patient-level summary (mean spot scores) ===\n"); print(patient_summary)

tests <- lapply(c("SMC","ECM","INFL","MS","OSTEO","PIEZO1"), function(vv) {
  d0 <- patient_summary[patient_summary$time=="D0", vv]; d7 <- patient_summary[patient_summary$time=="D7", vv]
  data.frame(score=vv, mean_D0=mean(d0), mean_D7=mean(d7), direction=ifelse(mean(d7)>mean(d0),"Up","Down"),
             paired_t_p=tryCatch(t.test(d7, d0, paired=TRUE)$p.value, error=function(e) NA))
})
tests <- do.call(rbind, tests)
cat("\n=== Visium paired tests (n=3 patients, EXPLORATORY) ===\n"); print(tests)
write.csv(tests, file.path(outdir,"visium_paired_tests.csv"), row.names=FALSE)
write.csv(patient_summary, file.path(outdir,"visium_patient_summary.csv"), row.names=FALSE)

# paired patient-level tests (n=3, exploratory)
pv <- function(v) { t <- pairwise.t.test(v, paste0(patient_summary$patient, "_", patient_summary$time),
                                        paired=FALSE, pool.sd=FALSE)$p.value; t[1,2] }

spots <- do.call(rbind, spot_all)
# patient-level PIEZO1 from spots
pz <- aggregate(PIEZO1 ~ file, data=spots, FUN=mean)
patient_summary$PIEZO1 <- pz$PIEZO1[match(patient_summary$file, pz$file)]
cat("=== Visium patient-level summary (mean spot scores) ===\n"); print(patient_summary)

tests <- lapply(c("SMC","ECM","INFL","MS","OSTEO","PIEZO1"), function(vv) {
  d0 <- patient_summary[patient_summary$time=="D0", vv]; d7 <- patient_summary[patient_summary$time=="D7", vv]
  data.frame(score=vv, mean_D0=mean(d0), mean_D7=mean(d7), direction=ifelse(mean(d7)>mean(d0),"Up","Down"),
             paired_t_p=tryCatch(t.test(d7, d0, paired=TRUE)$p.value, error=function(e) NA))
})
tests <- do.call(rbind, tests)
cat("\n=== Visium paired tests (n=3 patients, EXPLORATORY) ===\n"); print(tests)
write.csv(tests, file.path(outdir,"visium_paired_tests.csv"), row.names=FALSE)
saveRDS(list(spots=spots, summary=patient_summary, tests=tests), file.path(outdir,"visium_state.rds"))

# ---- GeoMx formal: locked sets with FDR (patient-blocked) ----
gdir <- "geo_pages/gsm299685"
ann <- do.call(rbind, lapply(list.files(gdir, full.names=TRUE), function(f) {
  L <- readLines(f)
  t <- sub("^!Sample_title = ", "", L[grepl("^!Sample_title", L)][1])
  s <- sub("^!Sample_supplementary_file_1 = .*/", "", L[grepl("^!Sample_supplementary_file_1", L)][1])
  p <- strsplit(t, "_")[[1]]
  data.frame(title=t, dcc=sub("\\.gz$","",sub("^GSM\\d+_","",s)), patient=p[1],
             time=sub("ROI-\\d+_","",p[3]), seg=p[4], stringsAsFactors=FALSE)
}))
x <- as.data.frame(read_excel(file.path(d,"GSE299685_Normalized_Counts.xlsx")))
rownames(x) <- x$gene; x$gene <- NULL
ann <- ann[ann$dcc %in% colnames(x), ]
m <- as.matrix(x[, ann$dcc])
ann$time <- factor(ann$time, c("D0","D7")); ann$seg <- factor(ann$seg)
cam_out <- list()
for (nm in names(sets)) {
  idx <- na.omit(match(sets[[nm]], rownames(m)))
  if (length(idx) < 2) { cat(nm, ": <2 genes in GeoMx panel\n"); next }
  design <- model.matrix(~ patient + seg + time, data=ann)
  cam <- camera(m, setNames(list(idx), nm), design=design, contrast="timeD7")
  cam_out[[nm]] <- data.frame(set=nm, nGenes=length(idx), Direction=cam$Direction, P=cam$PValue)
}
cam_g <- do.call(rbind, cam_out)
cam_g$FDR <- p.adjust(cam_g$P, "BH")
cat("\n=== GeoMx camera (locked sets, patient-blocked) ===\n"); print(cam_g)
write.csv(cam_g, file.path(outdir,"camera_geomx_formal.csv"), row.names=FALSE)
