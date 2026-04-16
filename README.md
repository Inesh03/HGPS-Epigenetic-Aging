# HGPS Epigenetic Aging — Systemic vs. Structural Tissues

**Bioinformatics J-Component | CSE4067 | VIT Chennai | Winter 2025-26**

Analyzes tissue-specific epigenetic aging in Hutchinson-Gilford Progeria Syndrome (HGPS)
using DNA methylation clock analysis and Lamin A LAD disruption mapping.

## Datasets (download from GEO — not stored in this repo)

| GEO ID | Description | Use |
|---|---|---|
| GSE149960 | HGPS + Control skin fibroblasts (EPIC 850k) | Primary analysis |
| GSE182991 | HGPS + Control whole blood (EPIC 850k) | Tissue comparison |
| GSE54334 | Lamin A ChIP-seq consensus LADs (hg19) | LAD annotation |

## Repository Structure
scripts/ — R analysis pipeline (03–09)
results/
figures/ — All output figures (PNG + PDF)
tables/ — All output CSVs
data/LAD/ — LaminA_consensus_LADs.bed (reference annotation)


## Pipeline

Run scripts in order from `scripts/`:
03_build_sample_sheets.R → Sample metadata
04_preprocess.R → Noob + BMIQ normalization
05_clock_analysis.R → Horvath Skin & Blood Clock
06_compute_EAA.R → EAA + residual EAA
07_statistical_analysis.R → Wilcoxon tests + scatter plots
08_LAD_analysis.R → LAD methylation analysis
09_EAA_LAD_correlation.R → Spearman correlation


## Requirements

R 4.x + Bioconductor 3.22. Key packages: `minfi`, `wateRmelon`,
`methylclock`, `GenomicRanges`, `rtracklayer`, `tidyverse`, `ggpubr`, `rstatix`
