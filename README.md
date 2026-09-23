# Urbanization-Ethnicity-and-Diet-Shape-Gut-Microbial-Structural-Variations
Gut Microbial Structural Variations Code
# Gut Microbiome Structural Variation Analysis

# Gut Microbiome Structural Variation Analysis

This repository contains the R scripts used to investigate microbial structural variation (SV) patterns and their associations with host ethnicity, residential environment, dietary factors, and other host characteristics.

Two types of structural variations are analyzed:

- **dSVs** — deletion structural variations
- **vSVs** — variable structural variations

The analytical workflow follows five major stages:

1. Data filtering
2. PERMANOVA
3. PCoA
4. Association models
5. Mediation analysis

---

# Analysis workflow

```text
Raw SV data and metadata
          |
          v
   Data filtering
          |
          v
      PERMANOVA
          |
          v
         PCoA
          |
          v
   Association models
          |
    +-----+-----+
    |     |     |
    v     v     v
Residency Ethnicity Diet
 Model 1   Model 2  Model 3
    \       |       /
     \      |      /
      +-----+-----+
            |
            v
    Mediation analysis
```

The rationale for this order is:

- **Data filtering** prepares the SV and metadata datasets.
- **PERMANOVA** first evaluates the overall contribution of host and environmental factors to dSV and vSV variation.
- **PCoA** then visualizes the global SV structure and group-level separation.
- **Model 1–3 analyses** identify specific SV associations with residency, ethnicity, and dietary factors.
- **Mediation analysis** integrates significant associations from the previous models to investigate potential `Ethnicity -> Diet -> SV` pathways.

---

# Repository structure

## 1. Data preprocessing

### `Data_filtering.R`

Preprocesses the SV and dietary datasets before downstream analyses.

Main steps include:

- matching samples across datasets;
- filtering dSV and vSV features according to ethnicity-specific prevalence;
- retaining SVs with a maximum ethnicity-specific prevalence of at least 10%;
- filtering dietary variables according to ethnicity-specific prevalence;
- generating final filtered datasets for downstream analyses.

Main outputs include:

```text
final_filtered_dsv.csv
final_filtered_vsv.csv
final_filtered_diet.csv
```

---

# 2. PERMANOVA

### `PERMANOVA.R`

PERMANOVA is performed immediately after data filtering to evaluate how much of the overall dSV and vSV variation can be explained by host, environmental, and dietary factors.

Species abundance is incorporated as a covariate using the first three abundance-derived PCoA axes:

```text
Abun_PCoA1
Abun_PCoA2
Abun_PCoA3
```

### dSV distance

dSV profiles are analyzed using:

```text
Binary Jaccard distance
```

### vSV distance

vSV values are min-max normalized and analyzed using:

```text
Euclidean distance
```

The PERMANOVA workflow includes:

- full multivariable models;
- species-abundance adjustment;
- group-level environmental effects;
- leave-one-group-out marginal R²;
- individual factor-level marginal R²;
- Benjamini-Hochberg correction;
- separate dSV and vSV analyses;
- combined dSV/vSV results.

The number of permutations is:

```r
permutations = 999
```

## Group-level PERMANOVA

Environmental variables are evaluated by predefined groups.

For each group, the marginal contribution is estimated by comparing:

```text
Full model
vs
Reduced model excluding that group
```

## Factor-level PERMANOVA

Individual factors are evaluated using marginal R²:

```text
Marginal R² =
Full model R² - Reduced model R²
```

## Diet classification-level PERMANOVA

Dietary variables are additionally grouped into dietary classifications.

Each classification is treated as a block, and variables within the same classification are jointly permuted using the same permutation order.

This preserves the correlation structure among dietary variables within each classification.

Representative outputs include:

```text
permanova_vsv_dsv_bray_results.csv

permanova_diet_factors_dsv.csv
permanova_diet_factors_vsv.csv
permanova_diet_factors_dSV_vSV_with_classification.csv

PERMANOVA_Diet_Classification_dSV.csv
PERMANOVA_Diet_Classification_vSV.csv
PERMANOVA_Diet_Classification_dSV_vSV_merged.csv
```

---

# 3. PCoA

### `PCoA.R`

PCoA is performed after PERMANOVA to visualize overall sample-level structural variation patterns and group separation.

The script uses:

```text
vsv.csv
dsv.csv
group.csv
```

Samples are matched across the three datasets before analysis.

## vSV preprocessing

vSV values are min-max scaled feature by feature:

```text
(x - min) / (max - min)
```

Missing values are replaced with zero after scaling.

## dSV preprocessing

Missing dSV values are replaced with zero.

Samples containing no variation are removed before distance calculation.

## Distance calculation

The following distance measures are used:

```text
vSV: Canberra distance
dSV: Jaccard distance
```

A combined distance matrix is then calculated as:

```text
Combined distance =
(vSV distance + dSV distance) / 2
```

## PCoA calculation

PCoA is performed using:

```r
capscale(d_comb ~ 1, add = TRUE)
```

The first two PCoA axes are extracted for visualization.

The percentage of explained variance for PCoA1 and PCoA2 is calculated from the eigenvalues.

## Group comparison

A PERMANOVA test is also displayed in the PCoA visualization:

```r
adonis2(
  d_comb ~ group,
  data = groups,
  permutations = 999
)
```

The figure contains:

- sample-level PCoA coordinates;
- group-specific colors;
- 95% confidence ellipses;
- PCoA1 marginal boxplots;
- PCoA2 marginal boxplots;
- PERMANOVA R²;
- PERMANOVA P value.

Main output:

```text
PCoA_result.pdf
```

---

# 4. Association models

After the global PERMANOVA and PCoA analyses, specific SV associations are evaluated using three groups of statistical models.

---

## Model 1 — Residency-associated SVs

### `Model1-dSVs.R`

Tests associations between residency and dSVs across the full cohort.

Main comparison:

```text
Rural vs Urban
```

The analysis includes:

- chi-squared or Fisher's exact screening;
- multivariable logistic regression;
- species abundance;
- sequencing coverage;
- age;
- sex;
- BMI;
- ethnicity;
- diet;
- medication;
- urbanization;
- general metadata;
- Benjamini-Hochberg FDR correction;
- cross-validation;
- marginal R² analysis.

---

### `Model1-vSVs.R`

Tests associations between residency and vSVs.

Two versions are included:

```text
Version A:
Coverage included as a covariate
without a coverage cutoff.

Version B:
Sensitivity analysis restricted to
species coverage >= 5X.
```

The analysis uses:

- Spearman screening;
- multivariable linear regression;
- species abundance and sequencing coverage;
- demographic and environmental covariates;
- FDR correction;
- cross-validation;
- marginal R² estimation.

---

### `Model1-dSVs within ethnicity.R`

Performs Rural vs Urban dSV analyses separately within each ethnicity.

The analysis includes:

- ethnicity-stratified screening;
- chi-squared or Fisher's exact testing;
- logistic regression;
- species abundance and coverage adjustment;
- demographic and environmental covariates;
- FDR correction.

---

### `Model1-vSVs within ethnicity.R`

Performs Rural vs Urban vSV analyses separately within each ethnicity.

The analysis includes:

- ethnicity-stratified Spearman screening;
- multivariable linear regression;
- species abundance;
- species coverage;
- host and environmental covariates;
- FDR correction.

---

# Model 2 — Ethnicity-associated SVs

### `Model2-dSVs.R`

Tests associations between ethnicity and dSVs.

Each ethnicity is compared with the remaining cohort:

```text
Target ethnicity vs all other ethnicities
```

Two versions are evaluated:

```text
Version A:
Coverage included as a covariate
without a coverage cutoff.

Version B:
Sensitivity analysis restricted to
coverage >= 5X.
```

The multivariable logistic model includes:

```text
Ethnicity
Species coverage
Species abundance
Age
Sex
BMI
Residency
Diet
Medication
Urbanization
General metadata
```

Additional analyses include:

- screening;
- FDR correction;
- 10-fold cross-validation;
- variable importance;
- marginal R² estimation.

---

### `Model2-vSVs.R`

Tests associations between ethnicity and vSVs.

The workflow includes:

- Spearman screening;
- multivariable linear regression;
- species abundance;
- sequencing coverage;
- host and environmental covariates;
- FDR correction;
- cross-validation;
- marginal R² analysis.

Both no-coverage-cutoff and ≥5X sensitivity analyses are included.

---

# Model 3 — Diet-associated SVs

### `Model3-dSVs.R`

Tests associations between individual dietary factors and dSVs.

The analysis follows two stages.

### Stage 1 — Raw screening

Each dietary variable is tested against each dSV using logistic regression.

### Stage 2 — Adjusted model

Candidate associations are tested using multivariable logistic regression.

Covariates include:

```text
Target dietary factor
Species coverage
Species abundance
Age
Sex
BMI
Ethnicity
Residency
Medication
Urbanization
General metadata
```

---

### `Model3-vSVs.R`

Tests associations between dietary factors and vSVs.

The workflow also follows two stages.

### Stage 1 — Raw screening

Each diet-vSV association is evaluated using linear regression.

### Stage 2 — Adjusted model

Candidate associations are tested using multivariable linear regression.

Covariates include:

```text
Target dietary factor
Species coverage
Species abundance
Age
Sex
BMI
Ethnicity
Residency
Medication
Urbanization
General metadata
```

---

# 5. Mediation analysis

### `Mediation analysis.R`

Mediation analysis is performed after the association models.

Candidate pathways follow:

```text
Ethnicity
    |
    v
   Diet
    |
    v
    SV
```

Candidate pathways are constructed by integrating significant associations from:

```text
Ethnicity-Diet
Ethnicity-dSV
Ethnicity-vSV
Diet-dSV
Diet-vSV
```

For each candidate pathway, the following paths are estimated:

```text
A path:
Ethnicity -> Diet

B path:
Diet -> SV

C path:
Ethnicity -> SV
(total effect)

C' path:
Ethnicity -> SV
(direct effect after adjustment for Diet)
```

The mediation analysis also estimates:

```text
ACME
Average Causal Mediation Effect

ADE
Average Direct Effect

Total Effect

Proportion Mediated
```

## Primary analysis

The primary model includes sequencing coverage and host/environmental covariates.

## Sensitivity analysis

The sensitivity model additionally includes:

```text
Species abundance
```

## SV-specific models

For dSVs:

```text
Logistic regression
```

For vSVs:

```text
Linear regression
```

Mediation inference uses:

```r
N_SIMS <- 5000
```

Results are written to:

```text
mediation_results/
```

---

# Recommended running order

The intended analysis order for this repository is:

```text
1. Data_filtering.R

2. PERMANOVA.R

3. PCoA.R

4. Model1-dSVs.R
5. Model1-vSVs.R
6. Model1-dSVs within ethnicity.R
7. Model1-vSVs within ethnicity.R

8. Model2-dSVs.R
9. Model2-vSVs.R

10. Model3-dSVs.R
11. Model3-vSVs.R

12. Mediation analysis.R
```

In summary:

```text
Filtering
   ↓
PERMANOVA
   ↓
PCoA
   ↓
Residency / Ethnicity / Diet models
   ↓
Mediation analysis
```

---

# Required R packages

The scripts use packages including:

```text
tidyverse
dplyr
tidyr
readr
readxl
stringr
purrr
ggplot2
patchwork
broom
openxlsx
vegan
brglm2
lmtest
caret
pROC
mediation
```

Example installation:

```r
install.packages(c(
  "tidyverse",
  "readxl",
  "ggplot2",
  "patchwork",
  "broom",
  "openxlsx",
  "vegan",
  "brglm2",
  "lmtest",
  "caret",
  "pROC",
  "mediation"
))
```

---

# Reproducibility

Most statistical analyses use:

```r
set.seed(12345)
```

The PCoA group-comparison analysis uses:

```r
set.seed(123)
```

Before running the scripts, verify that all required input files are available and that file paths match the local directory structure.

Some scripts contain local absolute paths that should be changed before running the workflow on another computer.

---

# Main output directories

Outputs are written to directories including:

```text
results/
results_vsv/
results_ethnicity/
mediation_results/
```

PCoA visualization:

```text
PCoA_result.pdf
```

PERMANOVA results are primarily stored as `.csv` files, while regression and mediation outputs are mainly stored as `.xlsx` files.





