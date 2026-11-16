# 13_figures.R — publication figures (svg + png, Liberation Sans, editable svg text)
suppressPackageStartupMessages({library(ggplot2); library(readxl)})
matplotlib_like <- theme_bw(base_size = 11) +
  theme(text = element_text(family = "Liberation Sans"),
        panel.grid.minor = element_blank())
figdir <- "results/figures"; dir.create(figdir, showWarnings=FALSE, recursive=TRUE)
save_fig <- function(p, name, w=6, h=5) {
  ggsave(file.path(figdir, paste0(name, ".svg")), p, width=w, height=h)
  ggsave(file.path(figdir, paste0(name, ".png")), p, width=w, height=h, dpi=300)
}

# ---- Fig 1: cross-cohort logFC scatter ----
meta <- read.csv("results/meta_analysis.csv")
fdr1 <- read.csv("results/de_GSE310013_full.csv")
fdr2 <- read.csv("results/de_GSE326710_full.csv")
meta$sig310013 <- fdr1$adj.P.Val[match(meta$gene, fdr1$gene)] < 0.05
meta$sig326710 <- fdr2$adj.P.Val[match(meta$gene, fdr2$gene)] < 0.05
meta$cat <- with(meta, ifelse(sig310013 & sig326710, "Both FDR<0.05",
                       ifelse(sig310013, "Discovery only", ifelse(sig326710, "Replication only", "ns"))))
meta$cat <- factor(meta$cat, c("Both FDR<0.05","Discovery only","Replication only","ns"))
p1 <- ggplot(meta, aes(logFC_310013, logFC_326710, color=cat)) +
  geom_point(size=0.7, alpha=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey40") +
  scale_color_manual(values=c("Both FDR<0.05"="#0279EE","Discovery only"="#75A025",
                              "Replication only"="#FF9400","ns"="grey75"), name=NULL) +
  labs(x="logFC Day7 vs Day0 — GSE310013 (discovery, n=6)",
       y="logFC Day7 vs Day0 — GSE326710 (replication, n=4)") +
  matplotlib_like
save_fig(p1, "fig1_cross_cohort_scatter")

# ---- Fig 2: gene set NES forest across datasets ----
fg_bulk <- read.csv("results/fgsea_formal.csv")
fg_occ <- read.csv("results/fgsea_gse241205.csv")
nes <- rbind(
  data.frame(set=fg_bulk$set, NES=fg_bulk$NES, cohort=fg_bulk$dataset),
  data.frame(set=fg_occ$set, NES=fg_occ$NES, cohort="GSE241205_occlusion"))
nes$set <- factor(nes$set, c("SMC_CONTRACTILE","ECM_REMODELING","INFLAMMATION","OSTEOGENIC","MS_CHANNEL"))
nes$cohort <- factor(nes$cohort, c("GSE310013","GSE326710","GSE241205_occlusion"))
p2 <- ggplot(nes, aes(NES, set, color=cohort)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  geom_point(size=3, position=position_dodge(width=0.6)) +
  scale_color_manual(values=c("GSE310013"="#0279EE","GSE326710"="#75A025","GSE241205_occlusion"="#FF9400"), name=NULL) +
  labs(x="fgsea NES (Day7/D7-culture or occlusion vs baseline)", y=NULL) +
  matplotlib_like
save_fig(p2, "fig2_geneset_nes_forest")

# ---- Fig 3: Visium spatial maps, patient 1 ----
vs <- readRDS("results/visium_state.rds")
sp <- vs$spots; sp <- sp[sp$patient=="P1", ]
sp$time <- factor(sp$time, c("D0","D7"))
mk <- function(v, lab) {
  ggplot(sp, aes(px, py, color=.data[[v]])) + geom_point(size=0.4) +
    facet_wrap(~time, nrow=1) + scale_color_gradient2(midpoint=median(sp[[v]]),
      low="#0279EE", mid="#FAF9F3", high="#FF9400") +
    labs(x=NULL, y=NULL, color=lab, title=NULL) + coord_equal() +
    theme_bw(base_size=10) + theme(text=element_text(family="Liberation Sans"),
      strip.background=element_blank(), axis.ticks=element_blank(), axis.text=element_blank())
}
p3 <- mk("SMC", "SMC contractile score"); save_fig(p3, "fig3a_visium_P1_SMC", 7, 3.4)
p4 <- mk("INFL", "Inflammation score");   save_fig(p4, "fig3b_visium_P1_INFL", 7, 3.4)
p5 <- mk("PIEZO1", "PIEZO1 expression");  save_fig(p5, "fig3c_visium_P1_PIEZO1", 7, 3.4)

# ---- Fig 4: PIEZO1 across datasets (forest of logFC) ----
pz <- read.csv("results/piezo1_secondary.csv")
pz2 <- read.csv("results/piezo1_gse241205.csv")
vz <- aggregate(PIEZO1 ~ patient + time, data=vs$spots, FUN=mean)
vz_w <- reshape(vz, idvar="patient", timevar="time", direction="wide")
pz_v <- data.frame(dataset="GSE310351_Visium", logFC=mean(vz_w$PIEZO1.D7 - vz_w$PIEZO1.D0), SE=sd(vz_w$PIEZO1.D7 - vz_w$PIEZO1.D0)/sqrt(3), P=NA)
pzo <- read.csv("results/piezo1_gse241205.csv")
all_pz <- rbind(
  data.frame(dataset=c("GSE310013 (n=6)","GSE326710 (n=4)"), logFC=pz$logFC, SE=pz$SE, P=pz$P),
  data.frame(dataset="GSE241205 occlusion (n=3)", logFC=pzo$logFC, SE=pzo$SE, P=pzo$P),
  pz_v)
all_pz$up <- 1.96*all_pz$SE
p6 <- ggplot(all_pz, aes(logFC, dataset, color=dataset)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=logFC-up, xmax=logFC+up), height=0.15, linewidth=0.8) +
  geom_point(size=3.2) +
  scale_color_manual(values=c("#0279EE","#75A025","#FF9400","#FD9BED"), guide="none") +
  labs(x="PIEZO1 logFC (remodeled vs baseline)", y=NULL) +
  matplotlib_like
save_fig(p6, "fig4_piezo1_forest", 6, 3.2)
cat("figures done\n"); print(list.files(figdir))
