# Pan-Cancer CRISPR Dependency & Biomarker Analysis of PLK1

## Overview
This repository contains an end-to-end, automated computational oncology pipeline analyzing **PLK1 (Polo-like Kinase 1)** CRISPR-Cas9 dependency across human cancer lineages using Broad Institute DepMap (Chronos) data, along with co-dependency evaluations against **CDK7** and **CDH1**.

---

## Key Results & Figures

### 1. Lineage Essentiality Profile
Across all evaluated lineages, *PLK1* dependency scores tightly cluster around the canonical pan-essential benchmark line (**Score = -1.0**), confirming *PLK1* as an indispensable mitotic kinase across both solid and hematologic malignancies.

![PLK1 Lineage Boxplot](results/figures/plk1_lineage_boxplot.png)

### 2. Lineage Differential Dependency (Volcano Analysis)
Lineage-level differential sensitivity shifts remain constrained within a narrow range ($\Delta \approx -0.07$ to $+0.13$). Following Benjamini-Hochberg FDR correction, no single lineage demonstrates statistically significant outlier vulnerability ($q > 0.05$), reflecting universal mitotic essentiality.

![PLK1 Volcano Plot](results/figures/plk1_volcano_plot.png)

### 3. Co-Dependency & Biomarker Correlation
* **CDK7 (Positive Co-dependency Control)**: Demonstrates strong, statistically significant positive co-dependency ($r \approx 0.42$, $p < 10^{-15}$), validating that rapidly cycling models sensitive to mitotic arrest (*PLK1*) simultaneously rely heavily on transcriptional/cell-cycle progression (*CDK7*).
* **CDH1 (Negative Control)**: Demonstrates near-zero correlation ($r \approx -0.02$, $p > 0.05$), confirming orthogonal phenotypic independence.

![PLK1 Co-Dependency Scatter](results/figures/plk1_cdk7_cdh1_correlation.png)

---

## Statistical Summary

| Analysis | Contrast / Gene | Metric | Value | Biological Invariant |
| :--- | :--- | :--- | :--- | :--- |
| **Essentiality Benchmark** | *PLK1* Pan-Cancer | Median Chronos | $-0.91$ to $-0.94$ | Universal Essential Kinase |
| **Differential Sensitivity** | Pan-Cancer Lineages | Effect Size ($\Delta$) | $-0.07$ to $+0.13$ | Uniform dependency profile |
| **Lineage Significance** | Pan-Cancer Cohort | Adjusted FDR | $q > 0.05$ | Absence of tissue-restricted bias |
| **Co-dependency Control 1** | *CDK7* | Pearson $r$ | $\sim 0.42$ ($p < 10^{-15}$) | Synchronous cell-cycle vulnerability |
| **Co-dependency Control 2** | *CDH1* | Pearson $r$ | $\sim -0.02$ ($p > 0.05$) | Phenotypic independence |

---

## Project Structure
```text
├── data/
│   ├── raw/                # Extracted DepMap raw benchmarks & metadata
│   └── processed/          # Merged and cleaned multi-gene datasets
├── results/
│   ├── figures/            # High-resolution publication plots
│   └── tables/             # ANOVA, Welch t-test & correlation tables
├── src/
│   ├── download_depmap.py          # Data ingestion and Chronos benchmarking
│   ├── clean_lineage_data.py       # Metadata mapping and lineage cleaning
│   ├── analyze_plk1.py             # Differential statistics & plotting
│   └── correlate_biomarkers.py     # Multi-gene co-dependency module
└── README.md