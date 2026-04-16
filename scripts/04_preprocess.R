library(minfi)
library(wateRmelon)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)

setwd("~/HGPS_Project")
options(timeout = 1200)

sheet_fibro <- read.csv("data/raw/GSE149960/sample_sheet.csv")
sheet_blood  <- read.csv("data/raw/GSE182991/sample_sheet.csv")

sheet_fibro <- sheet_fibro[sheet_fibro$disease.state.ch1 %in%
                             c("Hutchinson-Gilford Progeria Syndrome (HGPS)", "normal"), ]

sheet_blood  <- sheet_blood[sheet_blood$disease.state.ch1 %in%
                              c("classical Hutchinson Gilford Progeria Syndrome (HGPS)", "Control"), ]

cat("Fibroblast samples:", nrow(sheet_fibro), "\n")
cat("Blood samples:", nrow(sheet_blood), "\n")


RGset_fibro <- read.metharray(sheet_fibro$Basename, force = TRUE)
RGset_blood  <- read.metharray(sheet_blood$Basename,  force = TRUE)

cat("Fibroblast RGset:", dim(RGset_fibro), "\n")
cat("Blood RGset:",      dim(RGset_blood),  "\n")


detP_fibro <- detectionP(RGset_fibro)
detP_blood  <- detectionP(RGset_blood)

failed_fibro <- rowMeans(detP_fibro > 0.01) > 0.1
failed_blood  <- rowMeans(detP_blood  > 0.01) > 0.1

cat("Failed probes — Fibroblast:", sum(failed_fibro), "\n")
cat("Failed probes — Blood:",      sum(failed_blood),  "\n")


MSet_fibro <- preprocessNoob(RGset_fibro)
MSet_blood  <- preprocessNoob(RGset_blood)

beta_fibro <- getBeta(MSet_fibro)
beta_blood  <- getBeta(MSet_blood)

beta_fibro <- beta_fibro[!failed_fibro, ]
beta_blood  <- beta_blood[!failed_blood,  ]

cat("Fibroblast beta matrix:", dim(beta_fibro), "\n")
cat("Blood beta matrix:",      dim(beta_blood),  "\n")


library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

anno_EPIC   <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
design_epic <- ifelse(anno_EPIC[rownames(beta_fibro), "Type"] == "I", 1, 2)
design_epic[is.na(design_epic)] <- 2

beta_fibro_bmiq <- BMIQ(beta_fibro, design.v = design_epic)


design_blood <- ifelse(anno_EPIC[rownames(beta_blood), "Type"] == "I", 1, 2)
design_blood[is.na(design_blood)] <- 2

beta_blood_bmiq <- BMIQ(beta_blood, design.v = design_blood)


save(beta_fibro_bmiq, sheet_fibro,
     file = "data/processed/beta_fibroblasts_processed.RData")
save(beta_blood_bmiq, sheet_blood,
     file = "data/processed/beta_blood_processed.RData")

cat("Fibroblast BMIQ matrix:", dim(beta_fibro_bmiq), "\n")
cat("Blood BMIQ matrix:",      dim(beta_blood_bmiq),  "\n")
cat("Saved to data/processed/\n")


cat("Saved processed beta matrices to data/processed/\n")


class(beta_fibro_bmiq)
names(beta_fibro_bmiq)


dim(beta_fibro_bmiq$nbeta)
dim(beta_blood_bmiq$nbeta)


beta_fibro_norm <- beta_fibro_bmiq$nbeta
beta_blood_norm  <- beta_blood_bmiq$nbeta

cat("Fibroblast normalised matrix:", dim(beta_fibro_norm), "\n")
cat("Blood normalised matrix:",      dim(beta_blood_norm),  "\n")


save(beta_fibro_norm, sheet_fibro,
     file = "data/processed/beta_fibroblasts_processed.RData")
save(beta_blood_norm, sheet_blood,
     file = "data/processed/beta_blood_processed.RData")

cat("Resaved with correct beta matrices\n")


cat("Beta value range — Fibroblast:", range(beta_fibro_norm, na.rm = TRUE), "\n")
cat("Beta value range — Blood:",      range(beta_blood_norm,  na.rm = TRUE), "\n")
cat("NA count — Fibroblast:", sum(is.na(beta_fibro_norm)), "\n")
cat("NA count — Blood:",      sum(is.na(beta_blood_norm)),  "\n")

