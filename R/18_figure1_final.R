# 18_figure1_final.R — Manuscript Figure 1: cross-cohort effect-size agreement
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Revisions per advisor review: (1) robust core (4,148 same-direction) distinguished from
# the 203 direction-discordant doubly significant genes; (2) CIs removed (genes are not
# independent units; see Methods 2.6) — plain r, rho, slope, R2, n only; (3) horizontal
# layout, legend at bottom, subdued grid.
suppressPackageStartupMessages(library(ggplot2))
figdir <- "results/figures"; dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

meta <- read.csv("results/meta_analysis.csv")
f1 <- read.csv("results/de_GSE310013_full.csv")
f2 <- read.csv("results/de_GSE326710_full.csv")
meta$FDR_310013 <- f1$adj.P.Val[match(meta$gene, f1$gene)]
meta$FDR_326710 <- f2$adj.P.Val[match(meta$gene, f2$gene)]

sig1 <- meta$FDR_310013 < 0.05; sig2 <- meta$FDR_326710 < 0.05
same <- sign(meta$logFC_310013) == sign(meta$logFC_326710)
meta$cat <- ifelse(sig1 & sig2 & same, "Robust core (same direction, both FDR < 0.05)",
            ifelse(sig1 & sig2, "Doubly significant, direction-discordant",
                   "All other commonly detected genes"))
meta$cat <- factor(meta$cat, c("All other commonly detected genes",
                               "Robust core (same direction, both FDR < 0.05)",
                               "Doubly significant, direction-discordant"))
cat("\nCategory counts:\n"); print(table(meta$cat))

r  <- cor(meta$logFC_310013, meta$logFC_326710)
rho <- cor(meta$logFC_310013, meta$logFC_326710, method = "spearman")
fit <- lm(logFC_326710 ~ logFC_310013, data = meta)
sl  <- coef(summary(fit))["logFC_310013", "Estimate"]
r2  <- summary(fit)$r.squared
stat_lab <- sprintf("Pearson r = %.3f\nSpearman rho = %.3f\nSlope = %.3f\nR\u00b2 = %.3f\nn = %s genes",
                    r, rho, sl, r2, format(nrow(meta), big.mark = ","))

pd <- meta[order(meta$cat), ]  # draw core/discordant on top of grey
p <- ggplot(pd, aes(logFC_310013, logFC_326710, color = cat)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.55) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.7, se = FALSE) +
  scale_color_manual(values = c("All other commonly detected genes" = "grey72",
                                "Robust core (same direction, both FDR < 0.05)" = "#0279EE",
                                "Doubly significant, direction-discordant" = "#D55E00")) +
  annotate("text", x = -5.9, y = 8.6, label = stat_lab, hjust = 0, vjust = 1,
           size = 2.9, color = "grey15", lineheight = 0.95) +
  scale_x_continuous(breaks = -6:6) + scale_y_continuous(breaks = -6:9) +
  labs(x = "logFC, discovery cohort GSE310013 (day 7 vs day 0)",
       y = "logFC, replication cohort GSE326710 (day 7 vs day 0)",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 12, 6, 10))

ggsave(file.path(figdir, "fig1_cross_cohort_scatter_v3.svg"), p, width = 8.2, height = 6.0)
ggsave(file.path(figdir, "fig1_cross_cohort_scatter_v3.png"), p, width = 8.2, height = 6.0, dpi = 300)
cat(sprintf("\nStats: r=%.3f rho=%.3f slope=%.3f R2=%.3f n=%d\n", r, rho, sl, r2, nrow(meta)))
