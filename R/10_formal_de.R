# 10_formal_de.R — Formal paired DE, cross-cohort concordance, meta, gene sets
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Locked rules: PREREGISTRATION_gene_sets.md (2026-10-29)
suppressPackageStartupMessages({library(limma); library(edgeR);  } )
suppressPackageStartupMessages({library(DESeq2); library(fgsea)})
set.seed(20261029)

outdir <- "results"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
d <- "data"

# ---- locked gene sets ----
sets <- list(
  SMC_CONTRACTILE = c("MYH11","ACTA2","TAGLN","CNN1","MYL9","MYLK","LMOD1","RGS5"),
  ECM_REMODELING  = c("MMP9","TIMP1","COL1A1","POSTN","MMP2","COL3A1","FN1","LOX"),
  INFLAMMATION    = c("IL6","CXCL2","CCL2","VCAM1","ICAM1","CXCL8","CCL20"),
  MS_CHANNEL      = c("PIEZO1","PIEZO2","TMC4","TMC5","TMC6","TMC8","TMEM63A","TMEM63B"),
  OSTEOGENIC      = c("RUNX2","SOX9","SPP1","BGLAP","BMP2")
)

# ---- load cohorts ----
load_cohort <- function(which) {
  if (which == "GSE310013") {
    x <- read.delim(gzfile(file.path(d,"GSE310013_raw_counts.txt.gz")), check.names=FALSE)
    rownames(x) <- x$Geneid; x$Geneid <- NULL
    time <- factor(ifelse(grepl("D7", colnames(x)), "D7", "D0"), levels=c("D0","D7"))
    patient <- factor(sub("D[07]_001$", "", colnames(x)))
    sym <- NULL
  } else {
    x <- read.csv(gzfile(file.path(d,"GSE326710_Raw_counts_all_samples.csv.gz")), check.names=FALSE)
    rn <- x[["Gene Symbol"]]; gid <- x[["Gene ID"]]; typ <- x[["Type"]]
    x[["Gene Symbol"]] <- NULL; x[["Type"]] <- NULL; x[["Gene ID"]] <- NULL
    x <- as.matrix(x); rownames(x) <- rn
    # collapse duplicated symbols: keep row with highest total counts
    if (anyDuplicated(rownames(x))) {
      tot <- rowSums(x)
      keep_rn <- !duplicated(rownames(x)) & !is.na(rownames(x)) & rownames(x) != ""
      for (g in unique(rownames(x)[duplicated(rownames(x)) & !is.na(rownames(x)) & rownames(x) != ""])) {
        i <- which(rownames(x) == g)
        keep_rn[i[which.max(tot[i])]] <- TRUE
      }
      x <- x[keep_rn & !is.na(rownames(x)) & rownames(x) != "", , drop=FALSE]
    }
    sym <- NULL  # rownames ARE symbols now
    time <- factor(ifelse(grepl("D7", colnames(x)), "D7", "D0"), levels=c("D0","D7"))
    patient <- factor(sub("_D[07] Read Count$", "", colnames(x)))
  }
  list(counts=x, time=time, patient=patient, sym=sym)
}

run_voom <- function(co) {
  y <- DGEList(counts=co$counts)
  keep <- filterByExpr(y, group=co$time)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ patient + time, data=data.frame(patient=co$patient, time=co$time))
  v <- voom(y, design, plot=FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef="timeD7", number=Inf, sort.by="none")
  tt$SE <- tt$logFC / tt$t
  tt$gene <- rownames(tt)
  tt$symbol <- tt$gene
  list(tt=tt, v=v, design=design, y=y)
}

c1 <- load_cohort("GSE310013"); r1 <- run_voom(c1)
c2 <- load_cohort("GSE326710"); r2 <- run_voom(c2)
cat("GSE310013 genes retained:", nrow(r1$tt), "; GSE326710:", nrow(r2$tt), "\n")

write.csv(r1$tt, file.path(outdir,"de_GSE310013_full.csv"), row.names=FALSE)
write.csv(r2$tt, file.path(outdir,"de_GSE326710_full.csv"), row.names=FALSE)

# ---- sensitivity: edgeR QL ----
run_edgeR <- function(co, label) {
  y <- DGEList(counts=co$counts, group=co$time)
  keep <- filterByExpr(y, group=co$time); y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ patient + time, data=data.frame(patient=co$patient, time=co$time))
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef="timeD7")
  tt <- topTags(qlf, n=Inf)$table
  data.frame(dataset=label, gene=rownames(tt), logFC=tt$logFC, P=tt$PValue, FDR=tt$FDR)
}
e1 <- run_edgeR(c1, "GSE310013"); e2 <- run_edgeR(c2, "GSE326710")

run_deseq <- function(co, label) {
  dds <- DESeqDataSetFromMatrix(round(as.matrix(co$counts)), data.frame(patient=co$patient, time=co$time),
                                ~ patient + time)
  dds <- DESeq(dds, quiet=TRUE)
  res <- as.data.frame(results(dds, contrast=c("time","D7","D0")))
  data.frame(dataset=label, gene=rownames(res), logFC=res$log2FoldChange, P=res$pvalue)
}
s1 <- run_deseq(c1, "GSE310013"); s2 <- run_deseq(c2, "GSE326710")

# ---- direction concordance ----
# NOTE (comment corrected in round 20): the "empty intersection" remark below
# dates from the early pilot era and referred to an interim DESeq2-based check;
# it is NOT true of the final primary voom model, whose both-cohort
# FDR<0.05 intersection is the 4,351-gene robust core. The primary metric per
# the preregistration spirit: discovery (GSE310013) FDR<0.05 genes -> direction
# agreement in replication (GSE326710); the strict both-FDR count is reported
# alongside it in direction_concordance.csv (n = 4,351).
conc <- function(a, b, sig_genes, label, rule) {
  if (length(sig_genes) == 0) return(data.frame(dataset=label, rule=rule, n_sig=0, concordant=NA, rate=NA, binomP=NA))
  k <- sum(sign(a[sig_genes,"logFC"]) == sign(b[sig_genes,"logFC"]))
  data.frame(dataset=label, rule=rule, n_sig=length(sig_genes), concordant=k, rate=k/length(sig_genes),
             binomP=binom.test(k, length(sig_genes), 0.5)$p.value)
}
fdr1 <- p.adjust(r1$tt$P.Value, "BH"); names(fdr1) <- rownames(r1$tt)
fdr2 <- p.adjust(r2$tt$P.Value, "BH"); names(fdr2) <- rownames(r2$tt)
disc <- names(fdr1)[fdr1 < 0.05]
res_conc <- rbind(
  conc(r1$tt, r2$tt, intersect(disc, rownames(r2$tt)), "voom_primary", "discFDR_rep_direction"),
  conc(r1$tt, r2$tt, intersect(disc, rownames(r2$tt))[fdr2[intersect(disc, rownames(r2$tt))] < 0.05], "voom_primary", "discFDR_repFDR_strict")
)
# sensitivity concordance
f1 <- p.adjust(e1$P, "BH"); names(f1) <- e1$gene; f2 <- p.adjust(e2$P, "BH"); names(f2) <- e2$gene
a1 <- setNames(e1$logFC, e1$gene); a2 <- setNames(e2$logFC, e2$gene)
disc_e <- names(f1)[f1 < 0.05]
res_conc <- rbind(res_conc, conc(data.frame(logFC=a1), data.frame(logFC=a2), intersect(disc_e, names(a2)), "edgeR_QL", "discFDR_rep_direction"))
p1 <- p.adjust(s1$P, "BH"); names(p1) <- s1$gene; p2 <- p.adjust(s2$P, "BH"); names(p2) <- s2$gene
b1 <- setNames(s1$logFC, s1$gene); b2 <- setNames(s2$logFC, s2$gene)
disc_s <- names(p1)[p1 < 0.05]
res_conc <- rbind(res_conc, conc(data.frame(logFC=b1), data.frame(logFC=b2), intersect(disc_s, names(b2)), "DESeq2", "discFDR_rep_direction"))
cat("\n=== Direction concordance ===\n"); print(res_conc)
write.csv(res_conc, file.path(outdir,"direction_concordance.csv"), row.names=FALSE)

# ---- inverse-variance meta on common genes (voom) ----
common <- intersect(rownames(r1$tt), rownames(r2$tt))
m1 <- r1$tt[common,"logFC"]; v1 <- r1$tt[common,"SE"]^2
m2 <- r2$tt[common,"logFC"]; v2 <- r2$tt[common,"SE"]^2
w1 <- 1/v1; w2 <- 1/v2
meta_logFC <- (m1*w1 + m2*w2)/(w1+w2)
meta_SE <- sqrt(1/(w1+w2))
meta_z <- meta_logFC/meta_SE
meta_P <- 2*pnorm(-abs(meta_z))
meta <- data.frame(gene=common, symbol=common,
                   logFC_310013=m1, logFC_326710=m2, logFC_meta=meta_logFC,
                   SE_meta=meta_SE, P_meta=meta_P, FDR_meta=p.adjust(meta_P,"BH"))
write.csv(meta, file.path(outdir,"meta_analysis.csv"), row.names=FALSE)

# ---- gene sets: camera (per cohort) + fgsea ----
cam_res <- list(); fg_res <- list()
for (nm in names(sets)) {
  for (i in 1:2) {
    v <- if (i==1) r1$v else r2$v; dsn <- if (i==1) r1$design else r2$design
    lab <- if (i==1) "GSE310013" else "GSE326710"
    idx <- na.omit(match(sets[[nm]], rownames(v$E)))
    if (length(idx) < 2) next
    cam <- camera(v$E, setNames(list(idx), nm), design=dsn, contrast="timeD7")
    cam_res[[paste(lab, nm, sep="_")]] <- data.frame(dataset=lab, set=nm, nGenes=length(idx),
                                                     Direction=cam$Direction, P=cam$PValue)
    ranks <- setNames(r1$tt[rownames(v$E),"logFC"], rownames(v$E)) # placeholder replaced below
    tt_use <- if (i==1) r1$tt else r2$tt
    ranks <- setNames(tt_use[rownames(v$E),"logFC"], rownames(v$E))
    fg <- fgsea(list(gs=sets[[nm]][sets[[nm]] %in% rownames(v$E)]), stats=ranks, eps=1e-50)
    fg_res[[paste(lab, nm, sep="_")]] <- data.frame(dataset=lab, set=nm, NES=fg$NES, pval=fg$pval,
                                                    padj=fg$padj, size=fg$size)
  }
}
cam_out <- do.call(rbind, cam_res)
cam_out$FDR <- unlist(lapply(c("GSE310013","GSE326710"), function(l)
  p.adjust(cam_out$P[cam_out$dataset==l], "BH")[cam_out$set[cam_out$dataset==l]])); fg_out <- do.call(rbind, fg_res)
cat("\n=== camera ===\n"); print(cam_out)
cat("\n=== fgsea ===\n"); print(fg_out)
write.csv(cam_out, file.path(outdir,"camera_formal.csv"), row.names=FALSE)
write.csv(fg_out, file.path(outdir,"fgsea_formal.csv"), row.names=FALSE)

# ---- PIEZO1 secondary (one-shot) ----
pie <- data.frame(
  dataset=c("GSE310013","GSE326710"),
  logFC=c(r1$tt["PIEZO1","logFC"], r2$tt["PIEZO1","logFC"]),
  SE=c(r1$tt["PIEZO1","SE"], r2$tt["PIEZO1","SE"]),
  P=c(r1$tt["PIEZO1","P.Value"], r2$tt["PIEZO1","P.Value"]))
pie$FDR <- p.adjust(pie$P, "BH")
cat("\n=== PIEZO1 ===\n"); print(pie)
write.csv(pie, file.path(outdir,"piezo1_secondary.csv"), row.names=FALSE)

saveRDS(list(r1=r1$tt, r2=r2$tt, meta=meta, conc=res_conc), "results/de_formal_state.rds")
cat("\nDONE\n")
