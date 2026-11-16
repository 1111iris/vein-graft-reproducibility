# 17_figure3_final.R — Manuscript Figure 3: Visium program scores, paired day-0/day-7 per patient
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Data: results/visium_patient_summary.csv (patient-level means from 11_spatial.R)
# Stats: results/visium_paired_tests.csv (paired t-test, n=3, from 11_spatial.R)
suppressPackageStartupMessages(library(ggplot2))
figdir <- "results/figures"
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

ps <- read.csv("results/visium_patient_summary.csv")
pt <- read.csv("results/visium_paired_tests.csv")
pmap <- setNames(pt$paired_t_p, pt$score)

panels <- data.frame(
  set = c("SMC", "INFL", "OSTEO"),
  lab = c("Smooth muscle contractile", "Inflammatory", "Osteogenic"),
  stringsAsFactors = FALSE
)
panels$strip <- sprintf("%s\n(paired P = %.3f)", panels$lab, pmap[panels$set])

long <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  data.frame(program = panels$strip[i], patient = ps$patient, time = ps$time,
             score = ps[[panels$set[i]]])
}))
long$time <- factor(long$time, c("D0", "D7"))
long$program <- factor(long$program, levels = panels$strip)
stopifnot(identical(levels(long$program), panels$strip))  # guard panel order

p <- ggplot(long, aes(time, score, group = patient, color = patient)) +
  geom_point(size = 2.6) + geom_line(linewidth = 0.6) +
  facet_wrap(~program, scales = "free_y") +
  scale_color_manual(values = c(P1 = "#0279EE", P2 = "#FF9400", P3 = "#75A025")) +
  labs(x = NULL, y = "Patient-level program score", color = "Patient") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "bottom")

ggsave(file.path(figdir, "fig3_visium_program_scores.svg"), p, width = 8.2, height = 3.6)
ggsave(file.path(figdir, "fig3_visium_program_scores.png"), p, width = 8.2, height = 3.6, dpi = 300)
cat("Figure 3 written; panel order:\n"); print(levels(long$program))
