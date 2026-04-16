# 🧬 HGPS Epigenetic Aging: Systemic vs. Structural Tissues

<div align="center">

**Analyzing the Mechanistic Drivers of Tissue-Specific Epigenetic Aging in Hutchinson-Gilford Progeria Syndrome**

[![R](https://img.shields.io/badge/R-4.x-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.22-brightgreen?style=flat-square)](https://bioconductor.org/)
[![Platform](https://img.shields.io/badge/Platform-Illumina%20EPIC%20850k-orange?style=flat-square)](https://www.illumina.com/)
[![License](https://img.shields.io/badge/License-MIT-blue?style=flat-square)](LICENSE)
[![VIT Chennai](https://img.shields.io/badge/VIT-Chennai-red?style=flat-square)](https://chennai.vit.ac.in/)

*CSE4067 Bioinformatics · J-Component Project · Winter Semester 2025–26*

**Akash R (22MIA1065) · Inesh ST (22MIA1108)**

*M.Tech Business Analytics (5-Year Integrated) · SCOPE, VIT Chennai*

*Submitted to Dr. Softya Sebastian, Assistant Professor, SCOPE, VIT Chennai*

</div>

---

## 📋 Table of Contents

- [Background](#-background)
- [Project Objectives](#-project-objectives)
- [Datasets](#-datasets)
- [Repository Structure](#-repository-structure)
- [Methodology](#-methodology)
- [Key Results](#-key-results)
- [Requirements](#-requirements)
- [Reproducing the Analysis](#-reproducing-the-analysis)
- [Figures](#-figures)
- [Acknowledgements](#-acknowledgements)

---

## 🔬 Background

**Hutchinson-Gilford Progeria Syndrome (HGPS)** is a rare fatal accelerated-aging disorder caused by a *de novo* point mutation (`c.1824C>T`) in the *LMNA* gene. This mutation produces a truncated, permanently farnesylated Lamin A protein — **Progerin** — which progressively destroys the nuclear lamina, causing premature aging symptoms within the first two years of life. Median age of death is ~13–14 years, predominantly from cardiovascular disease.

A key function of the nuclear lamina is to anchor large blocks of heterochromatin at the nuclear periphery through **Lamina-Associated Domains (LADs)**. When Progerin disrupts the lamina, LAD-heterochromatin contacts are lost and the distinctive DNA methylation patterns within these domains erode.

This project tests a quantitative hypothesis:

> *Does the degree of LAD methylation compartmentalization loss in individual HGPS fibroblasts directly predict the magnitude of their epigenetic clock dysregulation?*

---

## 🎯 Project Objectives

| # | Objective | Output |
|---|---|---|
| 1 | Compute DNAm age using the **Horvath Skin & Blood Clock** across HGPS skin fibroblast and whole blood samples; calculate raw and residual Epigenetic Age Acceleration (EAA) | Figs 1–3 |
| 2 | Map all EPIC 850k CpG probes onto Lamin A consensus LAD coordinates and quantify methylation compartmentalization loss per sample | Figs 4–4b |
| 3 | Test whether the degree of LAD disruption in individual samples directly correlates with their residual EAA — a sample-level mechanistic link | Fig 5 |

---

## 📦 Datasets

> **Note:** Raw IDAT files and processed `.RData` matrices are **not stored in this repository** due to size. Download directly from GEO using the accession numbers below.

| GEO Accession | Description | Samples | Role |
|---|---|---|---|
| [GSE149960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149960) | HGPS + Control skin fibroblasts · EPIC 850k · Köhler et al. (2020) | 9 HGPS · 6 Control | Primary analysis dataset |
| [GSE182991](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182991) | HGPS + Control whole blood · EPIC 850k · Bejaoui et al. (2022) | 8 classical HGPS · 12 Control *(7 non-classical excluded)* | Systemic tissue comparator |
| [GSE54334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54334) | Lamin A ChIP-seq consensus LADs · hg19 · Normal human fibroblasts | — | LAD genomic reference |

**Download raw data:**
```bash
# Using GEOquery in R
library(GEOquery)
getGEOSuppFiles("GSE149960", baseDir = "data/raw/")
getGEOSuppFiles("GSE182991", baseDir = "data/raw/")
```

---

## 📁 Repository Structure

```
HGPS-Epigenetic-Aging/
│
├── 📂 scripts/                         # Analysis pipeline (run in order)
│   ├── 03_build_sample_sheets.R        # Step 1 — Parse IDAT filenames + GEO metadata
│   ├── 04_preprocess.R                 # Step 2 — Noob correction + BMIQ normalization
│   ├── 05_clock_analysis.R             # Step 3 — Horvath Skin & Blood Clock
│   ├── 06_compute_EAA.R                # Step 4 — Raw + residual EAA computation
│   ├── 07_statistical_analysis.R       # Step 5 — Wilcoxon tests + scatter plots
│   ├── 08_LAD_analysis.R               # Step 6 — LAD/inter-LAD methylation analysis
│   └── 09_EAA_LAD_correlation.R        # Step 7 — Spearman EAA ~ LAD correlation
│
├── 📂 results/
│   ├── figures/                        # All output figures (PNG + PDF)
│   │   ├── Fig1_DNAmAge_boxplot.*      # DNAm age: HGPS vs Control, by tissue
│   │   ├── Fig2_EAA_boxplot.*          # Raw EAA: HGPS vs Control, by tissue
│   │   ├── Fig2b_EAA_resid_boxplot.*   # Residual EAA (regression-corrected)
│   │   ├── Fig3_DNAmAge_scatter.*      # DNAm age vs chronological age scatter
│   │   ├── Fig4_LAD_methylation.*      # LAD vs inter-LAD mean beta values
│   │   ├── Fig4b_LAD_diff.*            # LAD compartmentalization difference score
│   │   └── Fig5_EAA_LAD_correlation.*  # Spearman correlation: EAA ~ LAD disruption
│   │
│   └── tables/                         # All output CSVs
│       ├── DNAmAge_results.csv         # Per-sample DNAm ages
│       ├── DNAmAge_EAA_results.csv     # DNAm ages + chronological ages + EAA
│       ├── EAA_stats.csv               # Wilcoxon results — raw EAA
│       ├── EAA_resid_stats.csv         # Wilcoxon results — residual EAA
│       ├── LAD_methylation_summary.csv # Per-sample LAD/inter-LAD means
│       ├── LAD_stats.csv               # Wilcoxon results — LAD compartments
│       └── EAA_LAD_correlation.csv     # Merged data for Spearman correlation
│
├── 📂 data/
│   └── LAD/
│       └── LaminA_consensus_LADs.bed   # ✅ Included — hg19 LAD reference (GSE54334)
│
└── README.md
```

---

## ⚙️ Methodology

The pipeline runs across seven sequential R scripts. Below is a summary of each step.

```
┌──────────────────┐    ┌──────────────────┐    ┌──────────────────┐
│  03_build_sample │    │  04_preprocess   │    │ 05_clock_analysis│
│    _sheets.R     │───▶│       .R         │───▶│       .R         │
│                  │    │                  │    │                  │
│ Parse IDAT names │    │  Noob + BMIQ     │    │ skinHorvath clock│
│ Join GEO pheno   │    │  normalization   │    │ (391 CpGs)       │
└──────────────────┘    └──────────────────┘    └────────┬─────────┘
                                                          │
                        ┌─────────────────────────────────┘
                        ▼
┌──────────────────┐    ┌──────────────────┐    ┌──────────────────┐
│ 06_compute_EAA   │    │ 07_statistical   │    │ 08_LAD_analysis  │
│       .R         │───▶│  _analysis.R     │    │       .R         │
│                  │    │                  │    │                  │
│ Raw EAA          │    │ Wilcoxon tests   │    │ findOverlaps()   │
│ Residual EAA     │    │ Scatter plots    │    │ LAD compartments │
└──────────────────┘    └──────────────────┘    └────────┬─────────┘
                                                          │
                        ┌─────────────────────────────────┘
                        ▼
                ┌──────────────────┐
                │ 09_EAA_LAD_corr  │
                │   elation.R      │
                │                  │
                │ Spearman ρ       │
                │ EAA ~ LAD diff   │
                └──────────────────┘
```

### Normalization
- **Noob** (Normal-exponential out-of-band) background correction via `minfi::preprocessNoob`
- **BMIQ** (Beta-Mixture Quantile Normalization) via `wateRmelon` to harmonize Type I / Type II probe distributions

### Epigenetic Clock
- **Horvath Skin & Blood Clock** (`skinHorvath`, 391 CpGs) via `methylclock::DNAmAge()`
- Specifically calibrated for fibroblast and blood cell types; validated on HGPS fibroblast lines

### EAA Computation
- **Raw EAA** = `DNAmAge − ChronologicalAge`
- **Residual EAA** = residuals from `lm(DNAmAge ~ ChronologicalAge)` fitted per tissue across all samples (HGPS + Control combined) — removes chronological age confounding

### LAD Analysis
- Lamin A consensus LAD coordinates imported from `LaminA_consensus_LADs.bed` (hg19)
- CpG → LAD mapping via `GenomicRanges::findOverlaps()`
- **Compartmentalization score** = `mean_LAD_beta − mean_iLAD_beta` per sample
- A score closer to 0 indicates greater loss of methylation compartmentalization

---

## 📊 Key Results

| Comparison | Tissue | Metric | p-value | Significance |
|---|---|---|---|---|
| HGPS vs Control DNAm age | Skin Fibroblast | Raw DNAmAge (Wilcoxon) | **0.0076** | * |
| HGPS vs Control DNAm age | Whole Blood | Raw DNAmAge (Wilcoxon) | 0.624 | ns |
| HGPS vs Control EAA | Skin Fibroblast | Raw EAA | 0.181 | ns |
| HGPS vs Control EAA | Whole Blood | Raw EAA | 0.208 | ns |
| HGPS vs Control EAA | Skin Fibroblast | Residual EAA | **0.0663** | (trend) |
| HGPS vs Control EAA | Whole Blood | Residual EAA | 0.384 | ns |
| LAD methylation | Skin Fibroblast | Mean LAD beta | 0.0663 | (trend) |
| inter-LAD methylation | Skin Fibroblast | Mean iLAD beta | 0.113 | ns |
| **LAD compartmentalization** | **Skin Fibroblast** | **LAD − iLAD difference** | **0.00759** | **\*\*** |
| **EAA ~ LAD disruption** | **Skin Fibroblast** | **Spearman ρ = −0.85** | **6 × 10⁻⁵** | **\*\*\*\*** |

> **Core finding:** A strong negative Spearman correlation (ρ = −0.85, p = 6×10⁻⁵) links LAD methylation compartmentalization loss to reduced residual EAA in skin fibroblasts — providing, to our knowledge, the first sample-level quantitative mechanistic evidence connecting Progerin-induced nuclear structural damage to epigenetic clock dysregulation.

---

## 🛠️ Requirements

### R Packages

```r
# Bioconductor packages
BiocManager::install(c(
  "minfi",
  "wateRmelon",
  "methylclock",
  "GenomicRanges",
  "GenomeInfoDb",
  "rtracklayer",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICmanifest"
))

# CRAN packages
install.packages(c(
  "tidyverse",
  "ggplot2",
  "ggpubr",
  "rstatix"
))
```

### System

| Requirement | Version |
|---|---|
| R | ≥ 4.0 |
| Bioconductor | 3.22 |
| RAM | ≥ 16 GB recommended (IDAT loading) |
| Disk | ≥ 10 GB (raw IDAT files) |
| OS | macOS / Linux (tested on macOS arm64) |

---

## ▶️ Reproducing the Analysis

```bash
# 1. Clone the repository
git clone https://github.com/Inesh03/HGPS-Epigenetic-Aging.git
cd HGPS-Epigenetic-Aging

# 2. Create the required data directories
mkdir -p data/raw data/processed

# 3. Download raw data from GEO (in R)
# See the Datasets section above for GEOquery commands

# 4. Copy LaminA_consensus_LADs.bed into data/LAD/
#    (already included in this repo under data/LAD/)

# 5. Run scripts in order from within R/RStudio
#    Set working directory to the project root: setwd("~/HGPS-Epigenetic-Aging")
source("scripts/03_build_sample_sheets.R")
source("scripts/04_preprocess.R")
source("scripts/05_clock_analysis.R")
source("scripts/06_compute_EAA.R")
source("scripts/07_statistical_analysis.R")
source("scripts/08_LAD_analysis.R")
source("scripts/09_EAA_LAD_correlation.R")
```

All output figures and tables are written to `results/figures/` and `results/tables/` respectively.

---

## 🖼️ Figures

| Figure | Description |
|---|---|
| **Fig 1** | DNAm age (Horvath Skin & Blood Clock) in HGPS vs Control, skin fibroblast and whole blood |
| **Fig 2** | Raw Epigenetic Age Acceleration (EAA = DNAmAge − ChronAge), by tissue and disease status |
| **Fig 2b** | Regression-corrected residual EAA, by tissue and disease status |
| **Fig 3** | DNAm age vs chronological age scatter plot with group-specific regression lines |
| **Fig 4** | Mean beta values in LAD vs inter-LAD regions, HGPS vs Control fibroblasts |
| **Fig 4b** | LAD − inter-LAD methylation difference (compartmentalization score) per sample |
| **Fig 5** | Spearman correlation: LAD disruption score vs residual EAA (ρ = −0.85, p = 6×10⁻⁵) |

---

## 🙏 Acknowledgements

- **Köhler et al. (2020)** — GSE149960 fibroblast dataset and LAD epigenetic analysis framework
- **Bejaoui et al. (2022)** — GSE182991 whole blood dataset
- **Horvath et al. (2018)** — Skin & Blood Clock development and HGPS validation
- **Bioconductor** community for `minfi`, `methylclock`, `GenomicRanges`, and `wateRmelon`

---

<div align="center">

*VIT Chennai · SCOPE · CSE4067 Bioinformatics · Winter 2025–26*

</div>
