# 19_figure2_final.R — Manuscript Figure 2: locked gene programs across cohorts
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Revisions per advisor review: (A) NES dot plot for the three bulk cohorts with P values
# shown per point (nominal for culture cohorts, BH-adjusted across sets for occluded grafts);
# (B) GeoMx camera results as a separate panel (P and FDR), OSTEOGENIC marked not tested.
# No confidence intervals are shown, so the plot is a dot plot, not a forest plot.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
figdir <- "results/figures"; dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

f1 <- read.csv("results/fgsea_formal.csv")
f2 <- read.csv("results/fgsea_gse241205.csv")

set_order <- c("SMC_CONTRACTILE", "INFLAMMATION", "OSTEOGENIC", "ECM_REMODELING", "MS_CHANNEL")
coh <- c("GSE310013 (discovery)", "GSE326710 (replication)", "GSE241205 (occluded grafts)")

a <- f1[f1$dataset == "GSE310013", ]; a$cohort <- coh[1]; a$ptype <- "nominal"
b <- f1[f1$dataset == "GSE326710", ]; b$cohort <- coh[2]; b$ptype <- "nominal"
c2 <- data.frame(set = f2$set, NES = f2$NES_occlusion, cohort = coh[3], ptype = "adjusted",
                 pval = f2$padj)  # padj: BH across the five sets in one run
d <- rbind(a[, c("set","NES","cohort","ptype","pval")], b[, c("set","NES","cohort","ptype","pval")], c2)
d$set <- factor(d$set, set_order); d$cohort <- factor(d$cohort, coh)
d$y <- as.numeric(d$set) + 0.22 * (as.integer(d$cohort) - 2)   # manual dodge
fmt <- function(p) ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p))
d$plab <- ifelse(d$ptype == "nominal", paste0("P = ", fmt(d$pval)),
                                         paste0("adj P = ", fmt(d$pval)))

pA <- ggplot(d, aes(NES, y, color = cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_point(size = 2.7) +
  geom_text(aes(x = 2.75, label = plab), hjust = 0, size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = c("#0279EE", "#75A025", "#FF9400")) +
  scale_y_continuous(breaks = seq_along(set_order), labels = set_order) +
  scale_x_continuous(limits = c(-2.7, 4.6), breaks = -2:2) +
  labs(x = "fgsea NES (day 7 vs day 0; occluded vs spare)", y = NULL,
       color = NULL, title = "Bulk cohorts (fgsea)") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom", legend.margin = margin(0,0,0,0),
        plot.title = element_text(size = 11, face = "plain"))

cm <- read.csv("results/camera_geomx_formal.csv")
cm$set <- factor(cm$set, set_order)
cm$neglogP <- -log10(cm$P)
cm$plab <- paste0("P = ", fmt(cm$P), ", FDR = ", fmt(cm$FDR))
nt <- set_order[!set_order %in% cm$set]  # sets not tested (<2 genes on panel)

pB <- ggplot(cm, aes(neglogP, as.numeric(set))) +
  geom_point(size = 2.7, color = "#0279EE") +
  geom_text(aes(label = plab), hjust = -0.08, size = 2.8, color = "grey15") +
  geom_text(data = data.frame(set = factor(nt, set_order)),
            aes(x = 0, y = as.numeric(set)), label = "not tested\n(<2 genes on panel)",
            hjust = 0, size = 2.8, color = "grey45", inherit.aes = FALSE) +
  scale_y_continuous(breaks = seq_along(set_order), labels = set_order) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 8, 16, 24)) +
  labs(x = "camera -log10(P)", y = NULL, title = "GeoMx (camera, BH FDR across sets)") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 11, face = "plain"))

p <- (pA | pB) + plot_layout(nrow = 1, widths = c(1.85, 1), guides = "collect") &
  theme(legend.position = "bottom", text = element_text(family = "Liberation Sans"))
p <- p + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))
ggsave(file.path(figdir, "fig2_geneset_nes_camera_v2.svg"), p, width = 11.6, height = 4.6)
ggsave(file.path(figdir, "fig2_geneset_nes_camera_v2.png"), p, width = 11.6, height = 4.6, dpi = 300)
cat("Figure 2 written\n")
