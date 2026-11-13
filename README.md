# BLUP using R and ASReml

Genetic analysis of **slaughter class** in Norwegian Angus cattle using a pedigree-based animal model, fitted both with **ASReml-R** and by solving the **mixed model equations (MME)** directly in R.

Course project for HFA350 – *From Phenotypes to Breeding Values* (December 2025).

## Model

y = Xβ + Z₁u + Z₂h + e

- **Fixed effects:** slaughter age, birth month, slaughter year
- **Random effects:** animal (additive genetic, with pedigree relationship matrix **A**), herd × birth year
- u ~ N(0, Aσ²_A), h ~ N(0, Iσ²_H), e ~ N(0, Iσ²_E)

## Scripts

Run them in this order:

| Script | What it does |
|---|---|
| `01_filter_animals.R` | Keeps Angus males, splits herd/year codes, samples 10 herds × 100 animals, derives ADG and AMG |
| `02_prune_pedigree.R` | Strict pedigree pruning (`ggroups::pruneped`), then keeps animals + parents + grandparents |
| `03_asreml_and_model_comparison.R` | Linear-regression screening of fixed effects, ASReml animal model, heritability with SE (delta method), convergence plots, ASReml vs MME EBV comparison |
| `04_mme_blup.R` | Builds X, Z, G, R and solves the MME by hand; EBVs, PEV, SE, accuracy, and model-fit checks |

## Key results

| Component | Estimate | SE |
|---|---|---|
| Additive genetic variance | 0.273 | 0.100 |
| Herd-year variance | 0.447 | 0.085 |
| Residual variance | 0.484 | 0.077 |
| Heritability (h²) | 0.36 | 0.092 |

- EBV correlation between ASReml and MME: **r = 0.927** (n = 2516 animals)
- Correlation of EBV with fixed-effect-corrected phenotype: 0.851
- Correlation of observed with predicted phenotype: 0.896

## Packages

`dplyr`, `readr`, `tidyverse`, `ggplot2`, `ggroups`, `kinship2`, `pedigreemm`, `Matrix`, `flextable`, `officer`, `openxlsx`, and `asreml` (needs a licence).

## Data

The phenotype and pedigree data came from Tyr (Norway) and are **not included** in this repository.
