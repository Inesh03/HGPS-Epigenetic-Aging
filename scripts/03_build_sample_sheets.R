library(minfi)
library(dplyr)
setwd("~/HGPS_Project")

build_sample_sheet <- function(idat_dir, pheno_df, disease_col, age_col = NULL) {
  files        <- list.files(idat_dir, pattern = "_Grn.idat.gz", full.names = FALSE)
  gsm_ids      <- sub("_.*", "", files)
  sentrix_ids  <- sub(".*_(\\d+)_R(\\d+)C(\\d+)_Grn.*", "\\1", files)
  sentrix_pos  <- sub(".*_(\\d+)_(R\\d+C\\d+)_Grn.*", "\\2", files)
  sheet <- data.frame(
    Sample_Name      = gsm_ids,
    Sentrix_ID       = sentrix_ids,
    Sentrix_Position = sentrix_pos,
    Basename         = file.path(idat_dir, sub("_Grn.idat.gz", "", files)),
    stringsAsFactors = FALSE
  )
  pheno_df$Sample_Name <- rownames(pheno_df)
  sheet <- left_join(sheet,
                     pheno_df[, c("Sample_Name", disease_col,
                                  if (!is.null(age_col)) age_col)],
                     by = "Sample_Name")
  return(sheet)
}

sheet_fibro <- build_sample_sheet(
  idat_dir    = "data/raw/GSE149960/GSE149960",
  pheno_df    = pheno_149960,
  disease_col = "disease state:ch1",
  age_col     = "characteristics_ch1.1"
)
write.csv(sheet_fibro, "data/raw/GSE149960/sample_sheet.csv", row.names = FALSE)
head(sheet_fibro)

sheet_blood <- build_sample_sheet(
  idat_dir    = "data/raw/GSE182991/GSE182991",
  pheno_df    = pheno_blood,
  disease_col = "disease state:ch1",
  age_col     = "age:ch1"
)
write.csv(sheet_blood, "data/raw/GSE182991/sample_sheet.csv", row.names = FALSE)
head(sheet_blood)

save(pheno_149960, sheet_fibro, file = "data/processed/pheno_fibroblasts.RData")
save(pheno_blood,  sheet_blood,  file = "data/processed/pheno_blood.RData")
save(pheno_ref,                  file = "data/processed/pheno_reference.RData")


