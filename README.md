# Sepsis longitudinal transcriptomics

Code accompanying the manuscript:

**Cross-sectional and longitudinal covariation between an HLA-II antigen-presentation transcriptional programme and a T-cell dysfunction-associated transcriptional signature in sepsis: a secondary analysis of two longitudinal cohorts**

## Overview

The analysis uses two longitudinal whole-blood transcriptomic cohorts:

- E-MTAB-5273 (primary cohort)
- GSE236713 (validation cohort)

A 19-gene HLA-II antigen-presentation transcriptional programme and a 14-gene T-cell dysfunction-associated transcriptional signature are quantified with `singscore`. The analysis includes fixed-timepoint Spearman correlations, paired within-patient change correlations, linear mixed-effects models, and MCPcounter Immune8 composition-adjusted analyses.

## Repository structure

- `analysis.R` — data preparation and statistical analyses
- `run.R` — runs the complete analysis and generates tables and figures
- `scripts/tables/` — Supplementary Tables S1-S3
- `scripts/figures/` — Figures 2-4 and Supplementary Figures S1-S4
- `sessionInfo.txt` — software environment used for the verified full run

Figure 1 is a study-design schematic and is not generated programmatically.

## Data

The source datasets are publicly available under accessions E-MTAB-5273 and GSE236713. Source data are not redistributed in this repository.

The scripts use local files when present and otherwise attempt retrieval from the original repositories. Recognised local files include:

- `Burnham_sepsis_discovery_normalised_231.txt`
- `E-MTAB-5273.sdrf.txt`
- `GSE236713_series_matrix.txt.gz`
- GPL17077 annotation files matching `GPL17077*.txt` or `GPL17077*.gz`

Local input files may be placed in the repository root. Alternatively, set the environment variables `EMTAB_DIR` and `GEO_LOCAL_DIR` to directories containing the E-MTAB-5273 and GEO files, respectively.

## Requirements

The analysis requires R and the following packages.

CRAN packages:

```r
install.packages(c(
    "lme4", "lmerTest", "emmeans", "dplyr", "stringr",
    "ggplot2", "patchwork"
))
```

Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c(
    "GEOquery", "Biobase", "singscore", "MCPcounter",
    "illuminaHumanv4.db", "AnnotationDbi"
))
```

The verified software environment is recorded in `sessionInfo.txt`.

## Running the analysis

From the repository root:

```r
source("run.R")
```

or from a shell:

```bash
Rscript run.R
```

`run.R` first executes `analysis.R`, then generates Supplementary Tables S1-S3 and Figures 2-4 and S1-S4.

## Outputs

Numerical analysis outputs are written to `outputs_publication/`. Publication tables are written to `tables/`, and figures to `figures/`.

Figures are generated as PNG and TIFF files. PDF files are also written when Cairo graphics support is available.

For Figures 2 and 3, least-squares lines are visual guides only. Reported association estimates and confidence intervals are based on Spearman rank correlations.

## Reproducibility

Bootstrap and permutation procedures use fixed random seeds defined in `analysis.R`. The analysis checks gene-set completeness, sample matching, duplicate patient-timepoint records, MCPcounter population availability, and overlap between MCPcounter markers and the study constructs.

The complete pipeline was successfully reproduced under R 4.5.0 on Windows 11. Package versions used for the verified run are listed in `sessionInfo.txt`.

## License

This code is released under the MIT License. See `LICENSE`.
