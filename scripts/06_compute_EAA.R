library(tidyverse)
library(GEOquery)

setwd("~/HGPS_Project")

# ── 1. Blood chronological ages ──────────────────────────────────────────────
pheno_blood$chron_age <- as.numeric(pheno_blood$`age:ch1`)
pheno_blood$disease_state <- gsub("disease state: ", "", pheno_blood$characteristics_ch1.2)

# ── 2. Fibroblast chronological ages ─────────────────────────────────────────
fibro_chron_ages <- c(
  GSM4518921 = 1.17, GSM4518922 = 1.25, GSM4518923 = 2.25,
  GSM4518924 = 8.42, GSM4518925 = 8.83, GSM4518926 = 37.00,
  GSM4518927 = 40.00, GSM4518928 = 1.17, GSM4518929 = 2.00,
  GSM4518930 = 11.00, GSM4518931 = 11.00, GSM4518932 = 4.67,
  GSM4518933 = 5.00,  GSM4518934 = 6.92,  GSM4518935 = 8.50
)

# ── 3. Load DNAmAge and merge ─────────────────────────────────────────────────
dnam <- read.csv("results/tables/DNAmAge_results.csv")

dnam$chron_age <- NA

# Fill fibroblast ages
fibro_idx <- dnam$Sample %in% names(fibro_chron_ages)
dnam$chron_age[fibro_idx] <- fibro_chron_ages[dnam$Sample[fibro_idx]]

# Fill blood ages
blood_idx <- dnam$Sample %in% rownames(pheno_blood)
dnam$chron_age[blood_idx] <- pheno_blood[dnam$Sample[blood_idx], "chron_age"]

# Add granular group label for non-classical
blood_nonclassical <- rownames(pheno_blood)[
  grepl("non-classical", pheno_blood$disease_state, ignore.case = TRUE)]
dnam$Group <- dnam$Status
dnam$Group[dnam$Sample %in% blood_nonclassical] <- "Non-classical"

# ── 4. Compute EAA ────────────────────────────────────────────────────────────
dnam$EAA <- dnam$DNAmAge - dnam$chron_age

# ── 5. Verify ─────────────────────────────────────────────────────────────────
stopifnot(!any(is.na(dnam$chron_age)))
cat("Complete cases:", nrow(dnam), "\n")
print(dnam[, c("Sample", "Tissue", "Group", "chron_age", "DNAmAge", "EAA")])

# ── 6. Save ───────────────────────────────────────────────────────────────────
write.csv(dnam, "results/tables/DNAmAge_EAA_results.csv", row.names = FALSE)
cat("Saved → results/tables/DNAmAge_EAA_results.csv\n")


# ── Regression-based EAA (per tissue) ────────────────────────────────────────
dnam_core <- dnam |> filter(Group %in% c("HGPS", "Control"))

dnam_core <- dnam_core |>
  group_by(Tissue) |>
  mutate(
    lm_resid = residuals(lm(DNAmAge ~ chron_age)),
    EAA_resid = lm_resid
  ) |>
  ungroup()

# Re-run Wilcoxon on residual EAA
stats_resid <- dnam_core |>
  group_by(Tissue) |>
  wilcox_test(EAA_resid ~ Group) |>
  add_significance()

print(stats_resid)
# See complete stats
print(stats_resid, n = Inf, width = Inf)

library(ggpubr); library(rstatix)

# Position p-value brackets correctly
stat_pos <- stats_resid |>
  mutate(y.position = max(dnam_core$EAA_resid) * 1.15)

p2b <- ggplot(dnam_core, aes(x = Group, y = EAA_resid, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(colour = Group), width = 0.15, size = 2.5) +
  stat_pvalue_manual(stat_pos, label = "p = {p}", tip.length = 0.01) +
  facet_wrap(~ Tissue, scales = "free_y") +
  scale_fill_manual(values  = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  scale_color_manual(values = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  labs(title    = "Epigenetic Age Acceleration (Residual EAA)",
       subtitle = "Horvath Skin & Blood Clock | Regression-corrected",
       x = "", y = "Residual EAA (years)") +
  theme_bw(base_size = 13) + theme(legend.position = "none")

ggsave("results/figures/Fig2b_EAA_resid_boxplot.pdf", p2b, width = 8, height = 5)
ggsave("results/figures/Fig2b_EAA_resid_boxplot.png", p2b, width = 8, height = 5, dpi = 300)
cat("Fig 2b saved\n")

# Save updated stats
write.csv(stats_resid, "results/tables/EAA_resid_stats.csv", row.names = FALSE)


