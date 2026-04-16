# =============================================================================
# Script:  
# Purpose: EAA statistics — Wilcoxon tests, per-tissue, with plots
# Input:   results/tables/DNAmAge_EAA_results.csv
# Output:  results/tables/EAA_stats.csv
#          results/figures/Fig2_EAA_boxplot.png
#          results/figures/Fig3_DNAmAge_scatter.png
# =============================================================================

library(tidyverse)
library(ggpubr)
library(rstatix)

setwd("~/HGPS_Project")
dnam <- read.csv("results/tables/DNAmAge_EAA_results.csv")

# ── Keep only classical HGPS vs Control (exclude non-classical) ───────────────
dnam_core <- dnam |> filter(Group %in% c("HGPS", "Control"))

# ── 1. Wilcoxon tests per tissue ──────────────────────────────────────────────
stats_results <- dnam_core |>
  group_by(Tissue) |>
  wilcox_test(EAA ~ Group) |>
  add_significance()

print(stats_results)
write.csv(stats_results, "results/tables/EAA_stats.csv", row.names = FALSE)

# ── 2. EAA boxplot per tissue (Fig 2) ─────────────────────────────────────────
p2 <- ggplot(dnam_core, aes(x = Group, y = EAA, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(colour = Group), width = 0.15, size = 2.5) +
  stat_pvalue_manual(
    stats_results, label = "p = {p}", y.position = max(dnam_core$EAA) * 1.1
  ) +
  facet_wrap(~ Tissue, scales = "free_y") +
  scale_fill_manual(values  = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  scale_color_manual(values = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  labs(title    = "Epigenetic Age Acceleration (EAA)",
       subtitle = "Horvath Skin & Blood Clock | HGPS vs Control",
       x = "", y = "EAA (DNAmAge - Chronological Age)") +
  theme_bw(base_size = 13) + theme(legend.position = "none")

ggsave("results/figures/Fig2_EAA_boxplot.pdf", p2, width = 8, height = 5)
ggsave("results/figures/Fig2_EAA_boxplot.png", p2, width = 8, height = 5, dpi = 300)
cat("Fig 2 saved\n")

# ── 3. DNAmAge ~ ChronAge scatter per tissue (Fig 3) ──────────────────────────
p3 <- ggplot(dnam_core, aes(x = chron_age, y = DNAmAge, colour = Group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  facet_wrap(~ Tissue, scales = "free") +
  scale_color_manual(values = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  labs(title    = "DNAm Age vs Chronological Age",
       subtitle = "Dashed line = perfect age prediction",
       x = "Chronological Age (years)", y = "DNAm Age (years)") +
  theme_bw(base_size = 13)

ggsave("results/figures/Fig3_DNAmAge_scatter.pdf", p3, width = 9, height = 5)
ggsave("results/figures/Fig3_DNAmAge_scatter.png", p3, width = 9, height = 5, dpi = 300)
cat("Fig 3 saved\n")
