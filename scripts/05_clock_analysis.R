library(methylclock)
library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("~/HGPS_Project")
load("data/processed/beta_fibroblasts_processed.RData")
load("data/processed/beta_blood_processed.RData")


ages_fibro <- DNAmAge(beta_fibro_norm, clocks = "skinHorvath")
ages_blood  <- DNAmAge(beta_blood_norm,  clocks = "skinHorvath")

cat("Fibroblast DNAm ages computed:", nrow(ages_fibro), "\n")
cat("Blood DNAm ages computed:",      nrow(ages_blood),  "\n")


fibro_df <- data.frame(
  Sample    = sheet_fibro$Sample_Name,
  DNAmAge   = ages_fibro$skinHorvath,
  Status    = ifelse(grepl("HGPS", sheet_fibro$disease.state.ch1), "HGPS", "Control"),
  Tissue    = "Skin Fibroblast"
)

blood_df <- data.frame(
  Sample    = sheet_blood$Sample_Name,
  DNAmAge   = ages_blood$skinHorvath,
  ChronAge  = as.numeric(sheet_blood$age.ch1),
  Status    = ifelse(grepl("HGPS", sheet_blood$disease.state.ch1), "HGPS", "Control"),
  Tissue    = "Whole Blood"
)

cat("Fibroblast age range:", range(fibro_df$DNAmAge, na.rm = TRUE), "\n")
cat("Blood age range:",      range(blood_df$DNAmAge, na.rm = TRUE), "\n")


combined_df <- bind_rows(
  fibro_df,
  blood_df[, c("Sample", "DNAmAge", "Status", "Tissue")]
)

p1 <- ggplot(combined_df, aes(x = Tissue, y = DNAmAge, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(color = Status), width = 0.15, size = 2.5) +
  stat_compare_means(aes(group = Status), method = "wilcox.test",
                     label = "p.format", label.y.npc = 0.95) +
  scale_fill_manual(values  = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  scale_color_manual(values = c("HGPS" = "#D62728", "Control" = "#1F77B4")) +
  labs(
    title    = "Epigenetic Age: HGPS vs Control",
    subtitle = "Horvath Skin & Blood Clock",
    x        = "",
    y        = "DNAm Age (years)"
  ) +
  theme_bw(base_size = 13)

ggsave("results/figures/Fig1_DNAmAge_boxplot.pdf", p1, width = 7, height = 5)
ggsave("results/figures/Fig1_DNAmAge_boxplot.png", p1, width = 7, height = 5, dpi = 300)
cat("Figure 1 saved\n")


write.csv(combined_df, "results/tables/DNAmAge_results.csv", row.names = FALSE)
cat("Results saved to results/tables/DNAmAge_results.csv\n")
 