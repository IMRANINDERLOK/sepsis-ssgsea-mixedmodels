# sepsis-ssgsea-mixedmodels
Reproducible R workflow for longitudinal sepsis whole-blood transcriptomics (GSE54514) with ssGSEA/GSVA program scoring, mixed-effects models, and external validation in GSE65682.

<img width="3200" height="2500" alt="Graphical Abstract" src="https://github.com/user-attachments/assets/888b4c3e-19cc-4a8c-a25b-8919d5180ef3" />


# An Interpretable Program-Level Framework for Longitudinal Sepsis Gene Expression

This repository contains the analysis code and key result tables/figures for a program-level, longitudinal analysis of whole-blood sepsis transcriptomics. The workflow integrates:

- Day-1 baseline differential expression (Survivor vs NonSurvivor)
- Program-level scoring using ssGSEA/GSVA (MSigDB Hallmark gene sets)
- Longitudinal mixed-effects modelling across Day 1–Day 5
- Program-linked gene prioritisation (DE at Day 1 + correlation with program scores)
- External validation in the MARS cohort (GSE65682) using 28-day mortality and endotype context

## Study datasets
**Discovery cohort:** GSE54514 (GPL6947; Illumina HumanHT-12 v4)  
**External validation cohort:** GSE65682 (GPL13667; Affymetrix HG-U219)

> Note: This is a computational secondary analysis of public datasets. Results should be interpreted within the scope of the available sampling design and metadata.

---

## Repository structure

- `scripts/`  
  R scripts for download, preprocessing, scoring, modelling, validation, and figure generation.

- `04_results_tables/`  
  Main output tables generated from the discovery cohort (CSV).

- `05_figures_main/`  
  Main figures (PDF/PNG).

- `06_figures_supp/`  
  Supplementary figures (PDF/PNG).

- `validation/`  
  Validation cohort processing, tables, and figures:
  - `validation/02_data_processed/`
  - `validation/04_results_tables/`
  - `validation/05_figures_main/`
  - `validation/06_figures_supp/`

---

## Key outputs (high level)

### Discovery cohort (GSE54514)
- Day-1 DE results (Survivor vs NonSurvivor)
- ssGSEA program scores (IFNG_RESPONSE, INFLAMMATORY_RESPONSE, OXPHOS)
- Mixed-effects model summaries (main model, covariate-adjusted model, Day1–3 sensitivity)
- Program-linked candidate gene lists and overlap summaries

### Validation cohort (GSE65682)
- ssGSEA program scores by 28-day mortality
- Wilcoxon tests with FDR correction across programs
- Endotype (MARS) stratification of inflammatory program scores

---

## Requirements

This project was developed in **R (>= 4.2 recommended)**.

### R packages
Core packages used include:
- GEOquery
- limma
- dplyr, tidyr
- GSVA, GSEABase
- lme4, lmerTest
- ggplot2
- pROC (optional; for ROC experiments if included)

Install packages (example):
```r
install.packages(c("dplyr","tidyr","ggplot2","pROC"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","limma","GSVA","GSEABase"))
install.packages(c("lme4","lmerTest"))
