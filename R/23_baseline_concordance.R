# 23_baseline_concordance.R
# Baseline directional-concordance definitions for Section 3.4, computed
# from the final DE result files. Run after 10_formal_de.R and
# 12_clinical.R, from the project root.
#
# Four definitions (advisor review, round 8):
#   (1) three-way   : discovery, replication, and occlusion all same direction
#   (2) meta vs occ : culture meta-analytic logFC vs occluded-graft logFC
#   (3) disc vs occ : discovery logFC vs occluded-graft logFC
#   (4) rep vs occ  : replication logFC vs occluded-graft logFC
# Each is reported genome-wide and for discovery-significant genes
# (discovery FDR < 0.05), the latter under both the discovery and the
# meta-analytic direction. Genes are not independent observations;
# these fractions are descriptive.

outdir <- "results"
m   <- read.csv(file.path(outdir, "meta_analysis.csv"))
occ <- read.csv(file.path(outdir, "de_GSE241205_fpkm_full.csv"))
d1  <- read.csv(file.path(outdir, "de_GSE310013_full.csv"))

mm <- merge(m, occ[, c("gene", "logFC")], by = "gene")
names(mm)[names(mm) == "logFC"] <- "logFC_occ"
mm <- merge(mm, d1[, c("gene", "adj.P.Val")], by = "gene")
names(mm)[names(mm) == "adj.P.Val"] <- "FDR_disc"

same <- function(a, b) sign(a) == sign(b)
n <- nrow(mm)

row <- function(def, vec, subset = rep(TRUE, n)) {
  # NA in vec (three-way definition: cohorts disagree) counts as discordant
  k <- sum(same(vec[subset], mm$logFC_occ[subset]), na.rm = TRUE)
  data.frame(definition = def, subset = if (all(subset)) "all common genes" else "discovery FDR<0.05",
             n = sum(subset), concordant = k, rate = k / sum(subset))
}
out <- rbind(
  row("three-way (disc, rep, occ same direction)",
      ifelse(same(mm$logFC_310013, mm$logFC_326710), mm$logFC_310013, NA) * 1),
  row("culture meta vs occlusion", mm$logFC_meta),
  row("discovery vs occlusion", mm$logFC_310013),
  row("replication vs occlusion", mm$logFC_326710)
)
ds <- mm$FDR_disc < 0.05
out <- rbind(out,
  row("discovery vs occlusion", mm$logFC_310013, ds),
  row("culture meta vs occlusion", mm$logFC_meta, ds))
out <- out[order(out$subset, out$definition), ]
write.csv(out, file.path(outdir, "baseline_concordance_definitions.csv"), row.names = FALSE)
cat("baseline concordance definitions:\n"); print(out, digits = 3)
