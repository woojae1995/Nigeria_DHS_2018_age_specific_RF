# 2018 Nigeria DHS Cluster-Aware Multistage Selection Analysis

This repository contains the author-generated R code used to reproduce the
statistical analysis, tables, and supplementary figures for a study using the
2018 Nigeria Demographic and Health Survey (NDHS).

The analysis applies a cluster-aware multistage selection (CMS) procedure to
identify individual-, household-, and community-level factors associated with
malaria parasitemia among children aged 6-59 months.

## Repository contents

- `Full_code_journal_submission.R`: complete analysis used for journal
  submission and public reproducibility.

The DHS microdata are not included in this repository because the DHS data-use
conditions do not permit users to redistribute downloaded datasets.

## Data access

This analysis uses de-identified data from the **2018 Nigeria Demographic and
Health Survey (NDHS)**. Researchers can request access through the
[DHS Program data portal](https://dhsprogram.com/data/available-datasets.cfm).
The dataset is identified by its survey title, year, and DHS filenames; no DOI
or accession number is assigned in the DHS data catalog.

To request the data:

1. Register for a DHS Program account.
2. Submit a research project description.
3. Request access to the 2018 Nigeria DHS.
4. Accept the applicable DHS data-use conditions.
5. After approval, download the required recode files.

The analysis requires these original DHS filenames:

| File | DHS recode type |
|---|---|
| `NGKR7BFL.DTA` | Children's Recode (KR) |
| `NGPR7BFL.DTA` | Household Member Recode (PR) |
| `NGHR7BFL.DTA` | Household Recode (HR) |

Relevant DHS documentation includes:
*Demographic and Health Surveys Standard Recode Manual for DHS7* 
(ICF, 2018; publication identifier DHSG4) 

*Guide to DHS Statistics, DHS-7, version 2* 


These documents describe the standard recode variables and indicator definitions
They are not part of the analytic code for this study.

## Software requirements

The code was tested using R 4.5.1. It requires the following R packages:

```r
install.packages(c(
  "dplyr", "ggplot2", "glmnet", "gridExtra", "haven", "labelled",
  "naniar", "openxlsx", "patchwork", "purrr", "sjlabelled",
  "stringr", "survey", "tibble", "tidyr"
))
```

## Full analysis settings

The journal-submission script uses the prespecified full settings:

- 200 Step 1 bootstrap iterations
- 200 Step 2 bootstrap iterations
- 10 cross-validation folds
- Full 144-combination tuning grid

Runtime depends on the computer and may be substantial. 

## Outputs

The script creates the following principal outputs under `outputs/`:

- `cms_tables_and_figures.xlsx`: manuscript and supplementary tables and
  figures.
- `figures/Figure_1S_model_stability.png`
- `figures/Figure_2S_residual_calibration.png`
- `results/cms_tuning_results_light.rds`
- `results/cms_tuning_results_full.rds`
- `results/cms_final_results.RData`
- `session_info.txt`: R version, platform, and package versions used for the
  run.

## Reproducibility notes

- Access approval must be obtained independently from the DHS Program.
- DHS filenames must remain unchanged unless the `input_files` mapping in the
  script is updated.
- The full analysis settings should be retained when reproducing reported
  results.
- Results may vary across materially different R or package versions despite
  use of a fixed seed. Consult `session_info.txt` when comparing runs.

## Citation

When using this code, please cite the associated article:

TBD


