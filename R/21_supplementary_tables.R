# 21_supplementary_tables.R
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Compile Supplementary Tables S2-S9 from results/ result CSVs.
# Round-6 advisor revision. Also fixes a write-time indexing bug in
# camera_formal.csv (culture-cohort FDR was written as NA because
# p.adjust() output was indexed by set names; BH across the five sets
# within each cohort is recomputed here and the corrected file replaces
# camera_formal.csv).

indir  <- "results"
outdir <- "results/supplementary_tables"
dir.create(outdir, showWarnings = FALSE)

rd <- function(f) read.csv(file.path(indir, f), stringsAsFactors = FALSE)

# ---- S2: full DE results, both culture cohorts (merged by gene) ----
d1 <- rd("de_GSE310013_full.csv"); d2 <- rd("de_GSE326710_full.csv")
common <- rd("meta_analysis.csv")$gene            # 13,168-gene analysis universe
s2 <- merge(d1, d2[, c("gene","symbol","logFC","SE","P.Value","adj.P.Val")],
            by = c("gene","symbol"), all = FALSE, suffixes = c("_disc", "_rep"))
s2 <- s2[s2$gene %in% common, ]
s2 <- s2[order(s2$P.Value_disc), ]
names(s2)[names(s2) == "P.Value_disc"] <- "P_disc"
names(s2)[names(s2) == "adj.P.Val_disc"] <- "FDR_disc"
names(s2)[names(s2) == "P.Value_rep"]  <- "P_rep"
names(s2)[names(s2) == "adj.P.Val_rep"] <- "FDR_rep"
# 95% CIs for the reported logFCs (locked plan rule 3: logFC with 95% CIs).
# Cohort-level intervals use t quantiles at the paired-design RESIDUAL degrees
# of freedom (GSE310013: 6 patients, ~ patient + time -> df = 5; GSE326710:
# 4 patients -> df = 3). limma's moderated inference uses df.total =
# df.residual + df.prior >= df.residual, so these intervals are WIDER than
# model-exact intervals: they are a conservative approximation, not an exact
# match to limma's t-based inference (round-20 review).
s2$CI95_low_disc  <- s2$logFC_disc - qt(0.975, 5) * s2$SE_disc
s2$CI95_high_disc <- s2$logFC_disc + qt(0.975, 5) * s2$SE_disc
s2$CI95_low_rep   <- s2$logFC_rep  - qt(0.975, 3) * s2$SE_rep
s2$CI95_high_rep  <- s2$logFC_rep  + qt(0.975, 3) * s2$SE_rep
stopifnot(nrow(s2) == 13168)
write.csv(s2, file.path(outdir, "TableS2_full_DE_both_culture_cohorts.csv"), row.names = FALSE)

# ---- S3: 4,148-gene robust core, per-cohort and meta statistics ----
core <- rd("robust_core_priority_all.csv")   # 4,148 rows
stopifnot(nrow(core) == 4148)
s3 <- core[, c("gene","logFC_310013","FDR_310013","logFC_326710","FDR_326710",
               "logFC_meta","SE_meta","P_meta","FDR_meta",
               "best_celltype","ct_SMC","ct_Endo","ct_Fibro","ct_Immune",
               "geomx_seg_lfc","seg_region","occ_logFC","occ_concordant")]
# 95% CI for the meta-analytic logFC (locked plan rule 3). The normal
# quantile is kept here because the inverse-variance meta-analysis P values
# (script 10) are z-based (normal), so the interval matches the inference.
s3$CI95_low_meta  <- s3$logFC_meta - qnorm(0.975) * s3$SE_meta
s3$CI95_high_meta <- s3$logFC_meta + qnorm(0.975) * s3$SE_meta
write.csv(s3, file.path(outdir, "TableS3_robust_core_4148_genes.csv"), row.names = FALSE)

# ---- S4: fgsea + camera, all cohorts/platforms ----
fg_bulk  <- rd("fgsea_formal.csv")          # culture cohorts, per-set runs
fg_occ   <- rd("fgsea_gse241205.csv")       # occluded cohort, one run
cam_bulk <- rd("camera_formal.csv")
cam_gx   <- rd("camera_geomx_formal.csv")

# fix culture-cohort camera FDR (BH across the five sets within each cohort)
cam_bulk$FDR <- ave(cam_bulk$P, cam_bulk$dataset, FUN = function(p) p.adjust(p, "BH"))
write.csv(cam_bulk, file.path(indir, "camera_formal.csv"), row.names = FALSE)  # corrected archive copy

fg_occ$dataset <- "GSE241205"
fg <- rbind(
  cbind(platform="Bulk RNA-seq", method="fgsea", cohort=fg_bulk$dataset, set=fg_bulk$set,
        n_genes=fg_bulk$size, NES=fg_bulk$NES, P=fg_bulk$pval,
        P_adjusted=fg_bulk$padj, p_mode="nominal (per-set run)"),
  cbind(platform="Bulk RNA-seq", method="fgsea", cohort=fg_occ$dataset, set=fg_occ$set,
        n_genes=fg_occ$size, NES=fg_occ$NES_occlusion, P=fg_occ$pval,
        P_adjusted=fg_occ$padj, p_mode="BH across the five sets (one run)")
)
cam_gx$dataset <- "GSE299685 (GeoMx)"
cam_bulk$platform <- "Bulk RNA-seq"; cam_gx$platform <- "GeoMx"
fg2 <- rbind(
  cbind(platform="Bulk RNA-seq", method="fgsea", cohort=fg_bulk$dataset, set=fg_bulk$set,
        n_genes=fg_bulk$size, stat=fg_bulk$NES, P=fg_bulk$pval,
        P_adjusted=fg_bulk$padj, stat_or_direction="NES", p_mode="nominal (per-set run)"),
  cbind(platform="Bulk RNA-seq", method="fgsea", cohort=fg_occ$dataset, set=fg_occ$set,
        n_genes=fg_occ$size, stat=fg_occ$NES_occlusion, P=fg_occ$pval,
        P_adjusted=fg_occ$padj, stat_or_direction="NES", p_mode="BH across the five sets (one run)")
)
cam2 <- rbind(
  data.frame(platform="Bulk RNA-seq", method="camera", cohort=cam_bulk$dataset, set=cam_bulk$set,
             n_genes=cam_bulk$nGenes, stat=cam_bulk$Direction, P=cam_bulk$P,
             P_adjusted=cam_bulk$FDR, stat_or_direction="Direction", p_mode="BH across the five sets within cohort"),
  data.frame(platform="GeoMx", method="camera", cohort=cam_gx$dataset, set=cam_gx$set,
             n_genes=cam_gx$nGenes, stat=cam_gx$Direction, P=cam_gx$P,
             P_adjusted=cam_gx$FDR, stat_or_direction="Direction", p_mode="BH across tested sets (OSTEOGENIC not tested: <2 genes)")
)
s4 <- rbind(fg2, cam2)
s4 <- s4[order(s4$method, s4$cohort, s4$set), ]
write.csv(s4, file.path(outdir, "TableS4_fgsea_camera_all_cohorts.csv"), row.names = FALSE)

# ---- S5: DE-method sensitivity (voom vs edgeR-QL vs DESeq2) ----
# Binomial test on gene-level directional concordance (H0: p = 0.5).
# R underflows extremely small P values to 0 or 4.9e-324; these are
# displayed as "<1e-300" rather than as exact zeros. Genes are not
# independent observations, so these P values are descriptive.
s5 <- rd("direction_concordance.csv")
s5$binomP_display <- ifelse(s5$binomP < 1e-300, "<1e-300", formatC(s5$binomP, format = "e", digits = 1))
names(s5)[names(s5) == "binomP"] <- "binomP_gene_level_directional"
write.csv(s5, file.path(outdir, "TableS5_DE_method_sensitivity.csv"), row.names = FALSE)

# ---- S6: leave-one-marker-out attribution (19 variants) ----
base_markers <- list(
  SMC        = c("ACTA2","MYH11","TAGLN","CNN1"),
  Endothelium= c("PECAM1","VWF","CDH5","CLDN5"),
  Fibroblast = c("DCN","LUM","COL1A2","GSN"),
  Immune     = c("PTPRC","LST1","S100A8","COTL1"),
  Adipocyte  = c("ADIPOQ","PLIN1","LPL")
)
variant <- unlist(lapply(names(base_markers), function(cn)
  paste(cn, base_markers[[cn]], sep = " | drop ")))
lomo <- rd("attribution_leave_one_marker_out.csv")
stopifnot(nrow(lomo) == 19)
s6 <- cbind(variant = variant, lomo)
s6$prop_fibroblast <- round(s6$Fibroblast / 4019, 4)
write.csv(s6, file.path(outdir, "TableS6_leave_one_marker_out_attribution.csv"), row.names = FALSE)

# ---- S7: attribution correlation-threshold sensitivity ----
s7 <- rd("attribution_threshold_sensitivity.csv")
write.csv(s7, file.path(outdir, "TableS7_attribution_threshold_sensitivity.csv"), row.names = FALSE)

# ---- S8: occlusion concordance, per patient and leave-one-patient-out ----
pp  <- rd("occlusion_concordance_per_patient.csv");  pp$analysis  <- "per-patient"
loo <- rd("occlusion_concordance_leave_one_out.csv"); names(loo)[1] <- "patient"; loo$analysis <- "leave-one-patient-out (2-patient fits)"
s8 <- rbind(pp[, c("analysis","patient","n_genes","concordant","prop")],
            loo[, c("analysis","patient","n_genes","concordant","prop")])
write.csv(s8, file.path(outdir, "TableS8_occlusion_concordance_per_patient.csv"), row.names = FALSE)

# ---- S9: top-30 candidates (culture-only ranking) + ranking comparison ----
cult <- rd("robust_core_priority_culture_only.csv")
cult <- cult[order(cult$culture_rank), ]
s9a <- cult[cult$culture_rank <= 30,
            c("culture_rank","gene_symbol","logFC_meta","FDR_meta","best_celltype",
              "occ_logFC","occ_concordant")]
# display underflow FDR (numerical zero) as "<1e-300", matching the
# max(FDR, 1e-300) floor used in the ranking score (Methods 2.8)
s9a$FDR_meta <- ifelse(s9a$FDR_meta < 1e-300, "<1e-300", as.character(s9a$FDR_meta))
write.csv(s9a, file.path(outdir, "TableS9a_top30_culture_only_ranking.csv"), row.names = FALSE)
s9b <- rd("top30_ranking_comparison.csv")
write.csv(s9b, file.path(outdir, "TableS9b_top30_ranking_comparison.csv"), row.names = FALSE)

# ---- S11: leave-one-gene-out attribution (advisor round 11, must-fix item 1) ----
# removes self-marker circularity: each core gene's own expression is excluded
# from any marker score containing it before the gene-vs-score correlation
s11a <- rd("attribution_logo_composition.csv")
write.csv(s11a, file.path(outdir, "TableS11a_leave_one_gene_out_composition.csv"), row.names = FALSE)
s11b <- rd("attribution_marker_member_genes.csv")
write.csv(s11b, file.path(outdir, "TableS11b_marker_member_core_genes.csv"), row.names = FALSE)
s11c <- rd("attribution_rule_comparison.csv")  # max-signed-r (primary) vs max-|r| assignment rule
write.csv(s11c, file.path(outdir, "TableS11c_attribution_rule_comparison.csv"), row.names = FALSE)

cat("written to", outdir, ":\n"); print(list.files(outdir))
cat("\ncamera FDR fix check (culture cohorts):\n"); print(cam_bulk[, c("dataset","set","P","FDR")])
