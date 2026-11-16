# 22_finalize_audit_summary.R
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Rebuild reproducibility_audit_summary.csv from the FINAL result files,
# run after 16_inference_audit.R. Script 15 writes the core (15-metric)
# summary before the attribution analyses (16) exist; the attribution
# rows were previously appended by hand and had drifted onto the old
# 4,147 denominator. This script derives every attribution entry from
# the final CSVs on the 4,019-gene valid-assignment denominator, so
# the summary, Table S6/S7, and manuscript Section 3.5 share one source.

indir <- "results"
rd <- function(f) read.csv(file.path(indir, f), stringsAsFactors = FALSE)

base <- rd("reproducibility_audit_summary.csv")   # 15 metrics from script 15
stopifnot(nrow(base) == 15)

lomo  <- rd("attribution_leave_one_marker_out.csv")     # 19 variants, counts
thr   <- rd("attribution_threshold_sensitivity.csv")    # threshold sensitivity
top30 <- rd("top30_ranking_comparison.csv")
pp    <- rd("occlusion_concordance_per_patient.csv")
loo   <- rd("occlusion_concordance_leave_one_out.csv")

DENOM <- 4019  # core genes with valid spatial assignments

fib_props <- lomo$Fibroblast / DENOM
fib_range <- sprintf("%.1f%%-%.1f%% (largest in %d/19)",
                     100 * min(fib_props), 100 * max(fib_props), 19)

t02 <- thr[thr$threshold == 0.2, ]
assigned <- sum(t02[c("SMC","Endothelium","Fibroblast","Immune","Adipocyte")])
labeled  <- sprintf("%d/%d (fibroblast largest)", assigned, DENOM)
stopifnot(assigned + t02$Unassigned == DENOM)

extra <- data.frame(
  metric = c("top30_overlap_culture_only_vs_occlusion_penalized",
             "occlusion_concordance_per_patient_range",
             "occlusion_concordance_leave_one_out_2pt_range",
             "attribution_fibroblast_share_lomo_range",
             "attribution_labeled_at_abs_r_gt_0.2",
             "attribution_denominator"),
  value  = c(sprintf("%d/30", sum(top30$in_both == TRUE)),
             sprintf("%.1f%%-%.1f%%", 100 * min(pp$prop), 100 * max(pp$prop)),
             sprintf("%.1f%%-%.1f%%", 100 * min(loo$prop), 100 * max(loo$prop)),
             fib_range,
             labeled,
             DENOM))

out <- rbind(base, extra)
write.csv(out, file.path(indir, "reproducibility_audit_summary.csv"), row.names = FALSE)
cat("reproducibility_audit_summary.csv rebuilt:", nrow(out), "metrics\n")
print(out)
