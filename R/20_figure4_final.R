# 20_figure4_final.R — Manuscript Figure 4: top thirty robust core candidates
# Run from the project root (the directory containing data/ and results/): Rscript <this script>.R
# Revisions per advisor review: legend title "Spatial marker-associated category"
# (not "Cell attribution"); genes ordered by the culture-only ranking; the four
# direction-discordant genes in occluded grafts highlighted on the open points.
suppressPackageStartupMessages(library(ggplot2))
figdir <- "results/figures"; dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

t30 <- read.csv("results/robust_core_priority_culture_only.csv")[1:30,]
stopifnot(all(diff(t30$culture_rank) > 0))  # culture-only ranking order

long_meta <- data.frame(gene = t30$gene_symbol, logFC = t30$logFC_meta,
                        celltype = t30$best_celltype)
long_occ  <- data.frame(gene = t30$gene_symbol, logFC = t30$occ_logFC,
                        status = ifelse(t30$occ_concordant, "Concordant", "Direction-discordant"))
long_meta$gene <- factor(long_meta$gene, levels = rev(t30$gene_symbol))  # rank 1 at top
long_occ$gene  <- factor(long_occ$gene,  levels = rev(t30$gene_symbol))

cell_col <- c(Fibroblast = "#0072B2", SMC = "#D55E00", Endothelium = "#009E73",
              Immune = "#CC79A7", Adipocyte = "#999999")

p <- ggplot() +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(data = long_occ, aes(logFC, gene, color = status), shape = 21,
             fill = "white", size = 2.1, stroke = 0.8) +
  geom_point(data = long_meta, aes(logFC, gene, fill = celltype), shape = 21,
             color = "grey20", size = 2.1) +
  scale_fill_manual(values = cell_col, name = "Spatial marker-associated category") +
  scale_color_manual(values = c("Concordant" = "grey55", "Direction-discordant" = "#FF9400"),
                     name = "Occluded grafts (open)") +
  labs(x = "logFC (day 7 vs day 0 in culture; occluded vs spare vein)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7.5),
        legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))

ggsave(file.path(figdir, "fig4_top30_dotplot_v2.svg"), p, width = 7.5, height = 7.6)
ggsave(file.path(figdir, "fig4_top30_dotplot_v2.png"), p, width = 7.5, height = 7.6, dpi = 300)
cat("Top30 cell types:\n"); print(table(t30$best_celltype))
cat("discordant in top30:", paste(t30$gene_symbol[!t30$occ_concordant], collapse = ", "), "\n")
