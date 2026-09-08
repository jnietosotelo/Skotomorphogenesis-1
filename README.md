# Skotomorphogenesis-1
# Evolution of skotomorphogenic development in seedlings of Zea and Tripsacum

This repository contains the data and scripts used to analyze skotomorphogenic development across different groups (**Maize, Teosinte, and Tripsacum**). It includes statistical analyses, dimensionality reduction, clustering, and figure generation.

---

## Project structure

```
├── data/
│   ├── EvolEscoDB.csv
│   ├── SR_primordium.csv
│   ├── Photoperiod.csv
│   ├── EvolEscoMean.csv
│   ├── EvolEscoMeanAlt.csv
│   ├── EvolEscoDBHCJ.csv
│
├── results/
│   ├── Tables (CSV)
│   ├── Figures (TIFF)
│
├── scripts/
│   ├── GroupAnalysis.R
│   ├── CorreGroups.R
│   ├── SR_primordium.R
│   ├── Photoperiod.R
│   ├── AltitudeCorr.R
│   ├── AltitudeCorrTeo.R
│   ├── MultAnalysis.R
│   ├── SkotoDays.R
|
└── README.md
```

---

##  Requirements

The analysis was performed using **R (version ≥ 4.4.3)**.

Install required packages:

```r
install.packages(c(
  "dplyr", "ggplot2", "tidyr", "cowplot",
  "agricolae", "ggfortify", "MASS",
  "patchwork", "ape"
))
```

---

##  How to reproduce the analysis

1. Clone this repository:

```bash
git clone <your-repo-url>
```

2. Open R or RStudio

3. Set working directory:

```r
setwd("path/to/repository")
```

4. Run the select script:

```r
source("scripts/GroupAnalysis.R")
```

```r
source("scripts/CorreGroups.R")
```

```r
source("scripts/SR_primordium.R")
```

```r
source("scripts/Photoperiod.R")
```

```r
source("scripts/AltitudeCorr.R")
```

```r
source("scripts/AltitudeCorrTeo.R")
```

```r
source("scripts/MultAnalysis.R")
```

```r
source("scripts/SkotoDays.R")
```

All outputs (tables and figures) will be saved in the `results/` folder.

---

## Scripts description

### `GroupAnalysis.R`

**Input:**
- `EvolEscoDB.csv`  
  Each row represents an individual seedling.

#### Outputs

| File | Description |
|------|------------|
| Table_S10_Summary_by_Group.csv | Summary statistics (Min, Max, Mean, SD, CV) by group |
| Table_S10_ANOVA_Tukey.csv | ANOVA followed by Tukey's HSD test among groups |
| Table_S679.csv | Mean trait values per accession |

#### Figures

| Figure | Description |
|--------|------------|
| Figure_1.tiff | Violin and bubble plots for all traits across Maize, Teosinte, and Tripsacum ||
| Figure_S5.tiff | Breakage percentage per accession for Maize and Teocinte|
| Figure_S6_ARCN_Maize.tiff | ARCN frequency distribution |
| Figure_S7_SR_Maize.tiff | SR frequency distribution |
| Figure_S8.tiff | Coleoptile breakage (%) comparison across Maize, Teosinte, and Tripsacum ||

---

### `CorreGroups.R`

**Input:**
- `EvolEscoDB.csv`  
  Each row represents an individual seedling.

#### Outputs

| File | Description |
|------|------------|
| Table_Correlation_Group_Full_Matrix.csv | Full correlation matrix by group, r and p value |
| Table_Correlation_Maize.csv | Correlation matrix for maize, r value |
| Table_Correlation_Teosinte.csv | Correlation matrix for teosinte, r value |
| Table_Correlation_Tripsacum.csv | Correlation matrix for tripsacum, r value |

---

### `SR_primordium.R`

**Input:**
- `PrimordiumSR.csv`  
  Each row represents an individual seedling.

#### Outputs

| File | Description |
|------|------------|
| Table_ANOVA_Tukey_Primordium.csv | ANOVA followed by Tukey's HSD test among varieties |

#### Figures

| Figure | Description |
|--------|------------|
| Figure_2.tiff | Bubble plots of seminal roots in maize varieties |

---

### `Photoperiod.R`

**Input:**
- `Photoperiod.csv`  
  Each row represents an individual seedling.

#### Outputs

| File | Description |
|------|------------|
| Table_S8_Summary_Photoperiod.csv | Summary statistics (Min, Max, Mean, SD, CV) by treatment |
| Table_S8_ANOVA_Tukey_Photoperiod.csv | ANOVA followed by Tukey's HSD test among treatment |

#### Figures

| Figure | Description |
|--------|------------|
| Figure_3.tiff | Violin and bubble plots for all traits for maize |
| Figure_4.tiff | Violin and bubble plots for all traits for teocinte|

---

### `AltitudeCorr.R`

**Input:**
- `EvolEscoMeanAlt.csv`  
  Each row represents mean values per accession.

#### Outputs

| File | Description |
|------|------------|
| Table_S12_Correlation_Altitude_Maize.csv | Correlation between altitude and maize traits (r and p-values) |
| Table_S12_Correlation_Altitude_Tripsacum.csv | Correlation between altitude and Tripsacum traits (r and p-values) |

#### Figures

| Figure | Description |
|--------|------------|
| Figure_6.tiff | Scatter plots (X–Y) for traits significantly correlated with altitude and maize |

---

### `AltitudeCorrTeo.R`

**Input:**
- `EvolEscoMeanAlt.csv`  
  Each row represents mean values per accession.

#### Outputs

| File | Description |
|------|------------|
| Table_S12_Correlation_Altitude_Teosinte.csv | Correlation between altitude and teosinte traits (r and p-values) |

#### Figures

| Figure | Description |
|--------|------------|
| Figure_S13.tiff | Scatter plots (X–Y) for traits significantly correlated with altitude in teosinte |

---

### `MultAnalysis.R`

**Input:**

For PCA and LDA:
- `EvolEscoDB.csv`  
  Each row represents an individual seedling.

For hierarchical clustering:
- `EvolEscoMean.csv`  
  Each row represents mean values per accession.
 - `EvolEscoDBHCJ.csv`  
  Each row represents an individual seedling assigned to a cluster. 

#### Outputs

| File | Description |
|------|------------|
| Table_S13_PCA_loadings.csv | PCA loadings (variable contributions) |
| Table_S13_PCA_variance_explained.csv | Variance explained by principal components |
| Table_S14_LDA_coefficients.csv | Linear discriminant coefficients |
| Table_Clusters_Assignation.csv | Cluster assignment for each accession |
| Table_S15_Summary_Clusters.csv | Summary statistics by cluster |
| Table_S15_ANOVA_Tukey.csv |  ANOVA followed by Tukey's HSD test among cluster |

#### Figures

| Figure | Description |
|--------|------------|
| Figure_6.tiff | PCA scatter plot (first two components) |
| Figure_7.tiff | LDA scatter plot first two discriminat functions |
| Figure_8_Cluster_circular.tiff | Circular dendrogram showing cluster relationships |
| Figure_9.tiff | Trait variation across clusters |
| Figure_S14.tiff | LDA first two discriminant functions |

### `SkotoDays.R`

**Input:**
- `SDDBD.csv`  
  Each row represents mean values per seedling.

#### Outputs

| File | Description |
|------|------------|
| Table_SkotoDays_Summary.csv | Summary statistics by Group |
| Table_Sko_ANOVA_Tukey.csv | NOVA followed by Tukey's HSD test among Group |

#### Figures

| Figure | Description |
|--------|------------|
| Figure_S15.tiff | Violin and bubble plots for all traits across Maize seven days, Teosinte seven days, and Tripsacum 17 days ||

---

---
##  Notes
- Each script runs independently.
- Scripts include automatic cleaning of the R environment, console, and plots before execution.
- You can execute only the script relevant to your analysis.
- All outputs (tables and figures) are automatically saved in the `results/` folder.
- All figures and tables are generated from raw data.
- Missing values are handled within scripts.
- The scripts reproduce all statistical analyses and figures.
- Some manuscript tables required manual formatting and post-processing
  for publication purposes and are therefore not included as generated outputs.

---

## 📄 Citation

If you use this repository, please cite the corresponding manuscript.
