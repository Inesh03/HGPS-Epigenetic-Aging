# =============================================================================
# Script:  08_LAD_analysis.R
# Project: HGPS Epigenetic Aging
# Purpose: Map EPIC CpGs to Lamin A LADs and compare methylation HGPS vs Control
# Input:   data/processed/beta_fibroblasts_processed.RData
#          data/LAD/LaminA_consensus_LADs.bed (GSE54334 consensus)
# Output:  results/tables/LAD_methylation_summary.csv
#          results/tables/LAD_stats.csv
#          results/figures/Fig4_LAD_methylation.png/.pdf
# Reference: Köhler et al. 2020 (PMC7249329)
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(minfi)
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
})

setwd("~/HGPS_Project")

# ── Load data ─────────────────────────────────────────────────────────────────
load("data/processed/beta_fibroblasts_processed.RData")
lad_consensus <- import("data/LAD/LaminA_consensus_LADs.bed")

# ── 1. Set genome and filter to standard chromosomes ─────────────────────────
genome(lad_consensus) <- "hg19"
lad_consensus <- keepStandardChromosomes(lad_consensus, pruning.mode = "coarse")
cat("LAD regions on standard chromosomes:", length(lad_consensus), "\n")

# ── 2. Build CpG GRanges — drop probes with missing coordinates ──────────────
anno      <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_filt <- anno[rownames(beta_fibro_norm), ]

valid           <- !is.na(anno_filt$chr) & !is.na(anno_filt$pos)
cat("Probes with valid coordinates:", sum(valid), "/", nrow(anno_filt), "\n")
cat("Probes dropped (no coords):   ", sum(!valid), "\n")

anno_filt       <- anno_filt[valid, ]
beta_fibro_filt <- beta_fibro_norm[rownames(anno_filt), ]

cpg_gr <- GRanges(
  seqnames = anno_filt$chr,
  ranges   = IRanges(start = anno_filt$pos, width = 1),
  cpg      = rownames(anno_filt)
)
genome(cpg_gr) <- "hg19"
cat("CpG GRanges built:", length(cpg_gr), "probes\n")

# ── 3. Annotate CpGs as LAD or inter-LAD ─────────────────────────────────────
lad_hits <- findOverlaps(cpg_gr, lad_consensus)
lad_idx  <- unique(queryHits(lad_hits))
ilad_idx <- setdiff(seq_len(length(cpg_gr)), lad_idx)

cat("CpGs in LADs:      ", length(lad_idx),  "\n")
cat("CpGs in inter-LADs:", length(ilad_idx), "\n")
cat("LAD CpG fraction:  ", round(length(lad_idx)/length(cpg_gr)*100, 1), "%\n")

# ── 4. Mean methylation per sample per compartment ───────────────────────────
status <- ifelse(grepl("HGPS", sheet_fibro$disease.state.ch1), "HGPS", "Control")

lad_summary <- data.frame(
  Sample    = colnames(beta_fibro_filt),
  Status    = status,
  LAD_mean  = colMeans(beta_fibro_filt[lad_idx,  ], na.rm = TRUE),
  iLAD_mean = colMeans(beta_fibro_filt[ilad_idx, ], na.rm = TRUE)
) |> mutate(LAD_iLAD_diff = LAD_mean - iLAD_mean)

print(lad_summary)

# ── 5. Wilcoxon tests ─────────────────────────────────────────────────────────
lad_stats <- lad_summary |>
  pivot_longer(c(LAD_mean, iLAD_mean, LAD_iLAD_diff),
               names_to  = "Compartment",
               values_to = "Methylation") |>
  group_by(Compartment) |>
  wilcox_test(Methylation ~ Status) |>
  add_significance()

print(lad_stats)
write.csv(lad_stats,   "results/tables/LAD_stats.csv",              row.names = FALSE)
write.csv(lad_summary, "results/tables/LAD_methylation_summary.csv", row.names = FALSE)

# ── 6. Fig 4 — LAD vs inter-LAD methylation ──────────────────────────────────
plot_df <- lad_summary |>
  pivot_longer(c(LAD_mean, iLAD_mean),
               names_to  = "Compartment",
               values_to = "Methylation") |>
  mutate(Compartment = recode(Compartment,
                              LAD_mean  = "LAD (Lamin A)",
                              iLAD_mean = "inter-LAD"
  ))

stat_pos <- lad_stats |>
  filter(Compartment %in% c("LAD_mean", "iLAD_mean")) |>
  mutate(
    Compartment = recode(Compartment,
                         LAD_mean  = "LAD (Lamin A)",
                         iLAD_mean = "inter-LAD"
    ),
    y.position = c(
      max(lad_summary$LAD_mean)  * 1.05,
      max(lad_summary$iLAD_mean) * 1.05
    )
  )

p4 <- ggplot(plot_df, aes(x = Status, y = Methylation, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(colour = Status), width = 0.15, size = 2.5) +
  stat_pvalue_manual(stat_pos, label = "p = {p}", tip.length = 0.01) +
  facet_wrap(~ Compartment) +
  scale_fill_manual(values  = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  scale_color_manual(values = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  labs(
    title    = "DNA Methylation in Lamin A LAD vs inter-LAD Regions",
    subtitle = "HGPS vs Control Skin Fibroblasts | GSE54334 consensus LADs",
    x = "", y = "Mean Beta Value"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("results/figures/Fig4_LAD_methylation.pdf", p4, width = 7, height = 5)
ggsave("results/figures/Fig4_LAD_methylation.png", p4, width = 7, height = 5, dpi = 300)
cat("Fig 4 saved\n")

# Add LAD_iLAD_diff panel to Fig 4
diff_stat <- lad_stats |>
  filter(Compartment == "LAD_iLAD_diff") |>
  mutate(y.position = max(lad_summary$LAD_iLAD_diff) * 0.85)

p4b <- ggplot(lad_summary, aes(x = Status, y = LAD_iLAD_diff, fill = Status)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(colour = Status), width = 0.15, size = 2.5) +
  stat_pvalue_manual(diff_stat, label = "p = {p}", tip.length = 0.01) +
  scale_fill_manual(values  = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  scale_color_manual(values = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  labs(
    title    = "Loss of LAD Methylation Compartmentalisation in HGPS",
    subtitle = "LAD - inter-LAD methylation difference | p = 0.008",
    x = "", y = "LAD mean beta - inter-LAD mean beta"
  ) +
  theme_bw(base_size = 13) + theme(legend.position = "none")

ggsave("results/figures/Fig4b_LAD_diff.pdf", p4b, width = 5, height = 5)
ggsave("results/figures/Fig4b_LAD_diff.png", p4b, width = 5, height = 5, dpi = 300)
cat("Fig 4b saved\n")

