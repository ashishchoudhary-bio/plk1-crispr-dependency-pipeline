# Systematic Pan-Cancer CRISPR Dependency and Functional Co-Dependency Profiling of *PLK1*

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DepMap](https://img.shields.io/badge/Data-DepMap%20Achilles%20Chronos-red.svg)](https://depmap.org/portal/)

---

## Abstract
**Background:** Polo-like Kinase 1 (*PLK1*) is an essential serine/threonine protein kinase required for centrosome maturation, bipolar spindle assembly, and mitotic exit. While frequently proposed as an oncogenic therapeutic vulnerability, systematic cross-lineage benchmarking is essential to delineate whether *PLK1* sensitivity represents a targetable cancer-specific vulnerability or an invariant pan-cancer essentiality.  
**Methods:** Utilizing genome-wide CRISPR-Cas9 loss-of-function screening data from the Broad Institute Cancer Dependency Map (DepMap Project Achilles, Chronos scoring framework), we engineered an automated computational oncology pipeline. We evaluated differential lineage dependencies across diverse cancer lineages using Welch's two-tailed *t*-tests and one-way ANOVA with Benjamini-Hochberg False Discovery Rate (FDR) control. Furthermore, pairwise functional co-dependency profiling was conducted against *CDK7* (transcriptional/cell-cycle kinase positive control) and *CDH1* (cell-adhesion tumor suppressor negative control).  
**Results:** Pan-cancer *PLK1* dependency scores clustered tightly around the universal essential benchmark (median Chronos = -0.91 to -0.94), with zero lineages exhibiting non-essentiality (Score $\ge 0.0$). Lineage-level differential sensitivity shifts remained restricted within a narrow band ($\Delta = -0.07$ to $+0.13$), yielding no statistically significant outliers following multiple testing correction ($q > 0.05$). Co-dependency analysis revealed a strong, statistically significant positive correlation with *CDK7* ($r \approx 0.42$, $p < 10^{-15}$), whereas *CDH1* exhibited complete phenotypic independence ($r \approx -0.02$, $p > 0.05$).  
**Conclusion:** *PLK1* represents an invariant pan-cancer common essential target across human malignancies. Rapidly cycling neoplastic models reliant on mitotic fidelity synchronously depend on *CDK7*-regulated transcriptional cycling, emphasizing that therapeutic intervention requires selective mitotic delivery rather than histology-guided stratification.

---

## 1. Introduction
Targeted oncological drug discovery relies heavily on distinguishing between lineage-selective dependencies (e.g., *SOX10* in melanoma, *AR* in prostate carcinoma) and ubiquitous housekeeping/mitotic machinery. *PLK1* regulates G2/M progression by phosphorylating substrates involved in spindle assembly checkpoint inactivation, chromosome segregation, and cytokinesis. 

Despite persistent interest in small-molecule PLK1 inhibition (e.g., Volasertib, Onvansertib), the clinical efficacy has frequently been hampered by systemic myelosuppression and narrow therapeutic windows. Here, we systematically interrogate the pan-cancer dependency profile of *PLK1* across 1,000+ cancer cell lines to establish its essentiality benchmark, lineage divergence, and molecular co-dependencies.

---

## 2. Materials and Computational Methods

### 2.1 Data Retrieval and Preprocessing
* **Screening Framework**: CRISPR knock-out dependency scores were retrieved from the Broad Institute Cancer Dependency Map (DepMap Project Achilles).
* **Algorithm Scale**: Dependency effects were calibrated via the **Chronos** algorithm:
  * **$0.0$**: Neutrality (knockout yields growth rates comparable to negative-control non-targeting sgRNAs).
  * **$-1.0$**: Universal Pan-Essential Benchmark (knockout induces lethality comparable to universal ribosomal and proteasomal housekeeping genes).
* **Cohort Harmonization**: DepMap model identifiers (`ModelID`) were mapped to primary disease annotations, lineage designations, and subtype metadata, discarding ambiguous/unannotated samples.

### 2.2 Statistical Inference & Hypothesis Testing
1. **Global Variance**: One-way Analysis of Variance (ANOVA) was implemented to assess between-group variance across primary cancer lineages against within-lineage residual variance.
2. **Lineage Differential Testing**: For each individual lineage cohort ($L$), differential dependency was evaluated against the background pan-cancer remainder cohort ($C \setminus L$) using Welch’s two-tailed unpaired $t$-test (robust against unequal sample sizes and heteroscedasticity):
   $$\Delta \text{ Chronos} = \mu_{\text{Lineage}} - \mu_{\text{Pan-Cancer}}$$
3. **Multiple Testing Correction**: P-values across lineages were adjusted using the Benjamini-Hochberg procedure to enforce False Discovery Rate control ($\alpha = 0.05$):
   $$q_i = \min \left( \frac{m}{i} p_{(i)}, q_{i+1} \right)$$

### 2.3 Genomic Biomarker & Functional Co-Dependency
To evaluate functional alignment, pairwise linear and monotonic co-dependency was calculated across matched cellular models:
* **Positive Co-Dependency Control (*CDK7*)**: Core catalytic component of TFIIH and CDK-activating kinase (CAK) complex.
* **Negative/Orthogonal Control (*CDH1*)**: Epithelial calcium-dependent adhesion molecule (E-cadherin).
* **Metrics**: Parametric Pearson correlation coefficient ($r$) and non-parametric Spearman rank correlation ($\rho$).

---

## 3. Results & Scientific Observations

### 3.1 Invariant Pan-Cancer Essentiality of *PLK1*
Distributional evaluation of *PLK1* Chronos scores across all examined cancer lineages demonstrates universal indispensability. 

* Median Chronos scores across solid and hematologic malignancies concentrated tightly between **$-0.75$ and $-1.05$**, straddling the canonical essential benchmark ($-1.0$).
* No cell line lineage demonstrated an interquartile range extending into the neutral or growth-tolerant zone ($\ge 0.0$), validating that mitotic progression in cultured tumor systems is universally *PLK1*-dependent.

![PLK1 Lineage Boxplot](results/figures/plk1_lineage_boxplot.png)

---

### 3.2 Uniform Lineage Sensitivity and Lack of Tissue-Specific Bias
To assess whether specific tumor types possess selective therapeutic vulnerabilities, differential effect sizes and significance metrics were mapped via volcano visualization.

* Effect size differentials ($\Delta \text{ Chronos}$) remained constrained to a narrow biological range ($-0.07 \le \Delta \le +0.13$).
* Following Benjamini-Hochberg adjustment, no lineage achieved statistically significant selective vulnerability ($q > 0.05$).
* *Observation*: Apparent uncorrected shifts in individual cohorts (e.g., prostate, fibroblast) reflect statistical sampling noise inherent to small sample sizes rather than genuine lineage-restricted hypersensitivity.

![PLK1 Volcano Plot](results/figures/plk1_volcano_plot.png)

---

### 3.3 Functional Co-Dependency: *CDK7* Coordination vs. *CDH1* Independence
Evaluating co-dependency profiles provides critical insights into the biological networks co-operating with *PLK1*.

* **Synchronous Mitotic & Transcriptional Vulnerability (*CDK7*)**:
  * *PLK1* and *CDK7* demonstrated a robust, statistically significant positive co-dependency ($r \approx 0.42$, $p < 10^{-15}$; Spearman $\rho \approx 0.41$).
  * *Biological Insight*: Neoplasms driven by rapid cell-cycle kinetics require high transcriptional output (driven by TFIIH/CDK7) in tandem with flawless centrosome dynamics (PLK1).
* **Orthogonal Phenotypic Independence (*CDH1*)**:
  * *CDH1* Chronos scores exhibited zero linear association with *PLK1* ($r \approx -0.02$, $p > 0.05$; Spearman $\rho \approx -0.02$).
  * *Biological Insight*: Epithelial-mesenchymal transition (EMT) dynamics and cellular adhesion status do not modulate baseline mitotic kinase lethality.

![PLK1 Co-Dependency Scatter](results/figures/plk1_cdk7_cdh1_correlation.png)

---

## 4. Quantitative Statistical Summary

| Contrast / Parameter | Gene / Cohort | Metric | Value | Scientific Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **Pan-Cancer Baseline** | *PLK1* | Median Chronos | $-0.91$ to $-0.94$ | Canonical common essential kinase |
| **Cross-Lineage Shift** | Pan-Cancer Lineages | Differential ($\Delta$) | $-0.07$ to $+0.13$ | Invariant cross-tissue dependency |
| **Multiple Testing** | All Primary Lineages | Benjamini-Hochberg | FDR $q > 0.05$ | Absence of lineage-restricted bias |
| **Functional Control 1** | *PLK1* vs. *CDK7* | Pearson $r$ / Spearman $\rho$ | $0.42$ / $0.41$ ($p < 10^{-15}$) | Significant cell-cycle co-dependency |
| **Functional Control 2** | *PLK1* vs. *CDH1* | Pearson $r$ / Spearman $\rho$ | $-0.02$ / $-0.02$ ($p = \text{n.s.}$) | Orthogonal phenotypic independence |
---

## 5. Conclusion & Translational Implications

1. **Universal Pan-Cancer Dependency**: Our systematic DepMap Chronos benchmarking conclusively demonstrates that *PLK1* operates as an invariant, pan-essential mitotic kinase across solid and hematologic tumor lineages. With median dependency scores tightly coupled around $-1.0$ and minimal differential variance ($\Delta \approx -0.07$ to $+0.13$, FDR $q > 0.05$), therapeutic sensitivity to *PLK1* knockout is fundamentally tissue-agnostic.
2. **Clinical Stratification Challenges**: Because *PLK1* lacks lineage-restricted hyper-vulnerabilities, conventional histology-guided patient selection is unlikely to yield preferential therapeutic windows. This molecular profile provides a mechanistic explanation for why historical small-molecule PLK1 inhibitors (e.g., Volasertib) encountered narrow therapeutic indices and dose-limiting myelosuppression.
3. **Co-Dependency & Combinatorial Rationale**: The strong functional synergy observed with *CDK7* ($r \approx 0.42$, $p < 10^{-15}$) indicates that tumors acutely dependent on mitotic spindle integrity simultaneously rely on transcriptional cycling machinery. Future therapeutic strategies should explore dual-pathway combinatorial regimens (e.g., PLK1 plus CDK7/transcription inhibition) or antibody-drug conjugates (ADCs) to bypass systemic toxicities while exploiting universal mitotic addiction.
---

## 6. Repository Structure
```text
├── data/
│   ├── raw/                # DepMap Achilles raw Chronos benchmarks & metadata
│   └── processed/          # Harmonized multi-gene dependency datasets
├── results/
│   ├── figures/            # High-resolution (300 DPI) publication visualizations
│   │   ├── plk1_lineage_boxplot.png
│   │   ├── plk1_volcano_plot.png
│   │   └── plk1_cdk7_cdh1_correlation.png
│   └── tables/             # Summary statistics, ANOVA, Welch t-test & correlation tables
├── src/
│   ├── download_depmap.py          # Data ingestion and Chronos matrix processing
│   ├── clean_lineage_data.py       # Metadata harmonization and lineage curation
│   ├── analyze_plk1.py             # Differential statistics & publication plotting
│   └── correlate_biomarkers.py     # Multi-gene co-dependency module
└── README.md
## References
1. Dempster, J. M., et al. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biology*, 22(1), 343.
2. Meyers, R. M., et al. (2017). Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. *Nature Genetics*, 49(12), 1779–1784.
3. Broad Institute Cancer Dependency Map (DepMap): [https://depmap.org/portal/](https://depmap.org/portal/)