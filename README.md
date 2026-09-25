# Charting Brain Structure in 22q11.2 Deletion Syndrome with Clinical Neuroimaging

Analysis code for:

> Jung B, et al. Charting brain structure in 22q11.2 deletion syndrome with clinical neuroimaging.


## Overview

This repository contains the Quarto notebooks used to generate the normative brain deviation scores, case-control analyses, imaging-transcriptomic analyses, and syndrome-specific growth charts reported in the paper.

## Data availability

The individual-level data analyzed in this study are **not included** in this repository and cannot be shared. The primary clinical cohort consists of protected health information from pediatric patients whose clinical MRIs were used for secondary research under an IRB-exempt protocol at the Children's Hospital of Philadelphia, and these data cannot be shared or transferred outside the institution.

Summary statistics underlying all analyses and figures are provided as Additional file 2 of the paper.

Third-party data used in this study can be obtained from their maintainers:

- **ENIGMA-22q:** access is governed by the ENIGMA 22q11.2 Deletion Syndrome Working Group (https://enigma.ini.usc.edu/ongoing/enigma-22q-working-group/).
- **Allen Human Brain Atlas:** publicly available (https://human.brain-map.org), accessed here via the `abagen` and `neuromaps` Python packages.
- **Lifespan Brain Chart Consortium growth charts:** available at https://brainchart.io/.

The notebooks are shared so that the analytic methods can be inspected and reused. They will not run end-to-end without the underlying data.

## Notebooks

Run in numerical order. Notebooks 01–03 prepare the analysis dataset; 04–10 produce the reported results.

| Notebook | Description | Outputs |
|---|---|---|
| `01_enigma_normative_centiles.qmd` | Generates population centiles for the ENIGMA-22q cohort using the Lifespan Brain Chart models | Processed ENIGMA dataset |
| `02_chop_normative_centiles.qmd` | Processes population centiles for the CHOP clinical cohort | Processed CHOP dataset |
| `03_merge_analysis_dataset.qmd` | Merges the CHOP and ENIGMA-22q datasets and defines analysis features | Combined analysis dataset |
| `04_case_control_effect_sizes.qmd` | Case-control effect sizes and sensitivity analyses | Fig. 1; Figs. S3, S6, S17–S23; Tables S3–S6, S16–S25 |
| `05_age_by_diagnosis_effects.qmd` | Age-by-diagnosis interaction effects | Fig. S4 |
| `06_sex_by_diagnosis_effects.qmd` | Sex-by-diagnosis interaction effects | Fig. S5 |
| `07_extreme_deviation_rates.qmd` | Rates of extreme deviations | Fig. 2; Figs. S7, S8; Tables S7, S8 |
| `08_gene_expression_correlation.qmd` | Spatial correlation between brain deviations and AHBA gene expression | Fig. 3; Fig. S9; Tables S9, S10 |
| `09_syndrome_specific_growth_charts.qmd` | Syndrome-specific (LMSz) growth charts | Fig. 4; Figs. S10, S16; Table S11 |
| `10_clinical_outcome_associations.qmd` | Associations with cognitive, language, and clinical features | Figs. S12, S13 |

