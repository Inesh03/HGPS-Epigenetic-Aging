# =============================================================================
# Script:  09_EAA_LAD_correlation.R
# Project: HGPS Epigenetic Aging
# Purpose: Correlate LAD methylation disruption with epigenetic age acceleration
# Input:   results/tables/DNAmAge_EAA_results.csv
#          results/tables/LAD_methylation_summary.csv
# Output:  results/figures/Fig5_EAA_LAD_correlation.png/.pdf
#          results/tables/EAA_LAD_correlation.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
})

setwd("~/HGPS_Project")

# ── Load inputs (both from previous scripts — no raw data needed) ─────────────
dnam     <- read.csv("results/tables/DNAmAge_EAA_results.csv")
lad_summ <- read.csv("results/tables/LAD_methylation_summary.csv")

# ── Compute residual EAA for fibroblasts ─────────────────────────────────────
dnam_fibro <- dnam |>
  filter(Group %in% c("HGPS", "Control"), Tissue == "Skin Fibroblast") |>
  mutate(EAA_resid = residuals(lm(DNAmAge ~ chron_age)))

# ── Clean sample IDs to match (strip array suffix from LAD summary) ───────────
lad_clean <- lad_summ |>
  mutate(Sample = gsub("_\\d{12}_R\\d{2}C\\d{2}$", "", Sample)) |>
  select(Sample, LAD_iLAD_diff)

# ── Merge ─────────────────────────────────────────────────────────────────────
corr_df <- left_join(dnam_fibro, lad_clean, by = "Sample")

cat("Merged rows:", nrow(corr_df), "| NAs in LAD_iLAD_diff:", 
    sum(is.na(corr_df$LAD_iLAD_diff)), "\n")

# ── Spearman correlation ──────────────────────────────────────────────────────
cor_test <- cor.test(corr_df$LAD_iLAD_diff, corr_df$EAA_resid, 
                     method = "spearman")
cat("Spearman rho:", round(cor_test$estimate, 3), 
    "| p-value:", round(cor_test$p.value, 4), "\n")

write.csv(corr_df, "results/tables/EAA_LAD_correlation.csv", row.names = FALSE)

# ── Fig 5 — correlation plot ──────────────────────────────────────────────────
p5 <- ggplot(corr_df, aes(x = LAD_iLAD_diff, y = EAA_resid, colour = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 3.5) +
  stat_cor(aes(label = after_stat(paste0("rho == ", ..r.., "~`,`~p == ", ..p..))),
           method = "spearman", label.x.npc = 0.05, label.y.npc = 0.95,
           colour = "black", size = 4, parse = TRUE) +
  scale_color_manual(values = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  labs(
    title    = "LAD Methylation Disruption vs Epigenetic Age Acceleration",
    subtitle = "Skin Fibroblasts | Spearman correlation",
    x        = "LAD - inter-LAD methylation difference",
    y        = "Residual EAA (years)",
    colour   = "Group"
  ) +
  theme_bw(base_size = 13)

ggsave("results/figures/Fig5_EAA_LAD_correlation.pdf", p5, width = 6, height = 5)
ggsave("results/figures/Fig5_EAA_LAD_correlation.png", p5, width = 6, height = 5, dpi = 300)
cat("Fig 5 saved\n")
