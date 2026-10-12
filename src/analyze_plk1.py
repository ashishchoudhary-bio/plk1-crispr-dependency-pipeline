"""
src/analyze_plk1.py
Statistical testing (ANOVA, Welch's t-test, Benjamini-Hochberg FDR)
and publication-ready figures for PLK1 CRISPR dependency.
"""
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

PROCESSED_FILE = "data/processed/plk1_pan_cancer_cleaned.csv"
FIGURES_DIR = "results/figures"
TABLES_DIR = "results/tables"

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

def run_analysis():
    df = pd.read_csv(PROCESSED_FILE)
    print(f"[Loaded] {len(df)} cell lines across {df['Lineage'].nunique()} lineages.")

    # Minimum 8 cell lines per lineage for statistical power
    counts = df["Lineage"].value_counts()
    valid_lineages = counts[counts >= 8].index.tolist()
    sub_df = df[df["Lineage"].isin(valid_lineages)].copy()

    # 1. One-way ANOVA across all lineages
    grouped_data = [group["PLK1_Chronos_Score"].dropna().values for _, group in sub_df.groupby("Lineage")]
    f_stat, anova_p = stats.f_oneway(*grouped_data)
    print(f"\n[ANOVA] F-statistic = {f_stat:.4f}, p-value = {anova_p:.4e}")

    # 2. Lineage-wise Welch's t-test vs remainder of cohort
    stats_records = []
    for lineage, group in sub_df.groupby("Lineage"):
        target_scores = group["PLK1_Chronos_Score"].dropna()
        other_scores = sub_df[sub_df["Lineage"] != lineage]["PLK1_Chronos_Score"].dropna()
        
        t_stat, p_val = stats.ttest_ind(target_scores, other_scores, equal_var=False)
        stats_records.append({
            "Lineage": lineage,
            "N_CellLines": len(target_scores),
            "Mean_Chronos": np.round(target_scores.mean(), 4),
            "Std_Chronos": np.round(target_scores.std(), 4),
            "T_Statistic": np.round(t_stat, 4),
            "P_Value": p_val
        })

    stats_df = pd.DataFrame(stats_records)

    # Benjamini-Hochberg FDR correction
    _, q_vals, _, _ = multipletests(stats_df["P_Value"], alpha=0.05, method="fdr_bh")
    stats_df["FDR_Q_Value"] = np.round(q_vals, 6)

    # Sort by strongest dependency (most negative Chronos score)
    stats_df = stats_df.sort_values(by="Mean_Chronos", ascending=True)

    # Save summary table
    table_path = os.path.join(TABLES_DIR, "plk1_lineage_stats.csv")
    stats_df.to_csv(table_path, index=False)
    print(f"[Saved Table] Summary statistics written to {table_path}")

    # 3. Figure 1: Ranked Dependency Boxplot
    ranked_lineages = stats_df["Lineage"].tolist()
    plt.figure(figsize=(10, max(6, len(ranked_lineages) * 0.35)))
    sns.boxplot(
        data=sub_df,
        y="Lineage",
        x="PLK1_Chronos_Score",
        order=ranked_lineages,
        palette="viridis",
        fliersize=2,
        width=0.6
    )
    plt.axvline(x=0, color="gray", linestyle="--", linewidth=1, label="Non-essential (0)")
    plt.axvline(x=-1.0, color="crimson", linestyle="--", linewidth=1.2, label="Essential Benchmark (-1.0)")
    plt.title("PLK1 CRISPR Dependency Across Cancer Lineages (Chronos Score)", fontsize=12, fontweight="bold")
    plt.xlabel("PLK1 Chronos Dependency Score (More Negative = Higher Essentiality)", fontsize=10)
    plt.ylabel("Lineage", fontsize=10)
    plt.legend(loc="lower left", frameon=True)
    plt.tight_layout()
    
    boxplot_path = os.path.join(FIGURES_DIR, "plk1_lineage_boxplot.png")
    plt.savefig(boxplot_path, dpi=300)
    plt.close()
    print(f"[Saved Figure] Boxplot saved to {boxplot_path}")

    # 4. Figure 2: Lineage Effect Size vs Significance (Volcano Plot)
    plt.figure(figsize=(8, 6))
    pan_cancer_mean = sub_df["PLK1_Chronos_Score"].mean()
    stats_df["Effect_Diff"] = stats_df["Mean_Chronos"] - pan_cancer_mean
    stats_df["neg_log10_p"] = -np.log10(stats_df["P_Value"] + 1e-12)

    sns.scatterplot(
        data=stats_df,
        x="Effect_Diff",
        y="neg_log10_p",
        size="N_CellLines",
        sizes=(40, 200),
        hue="FDR_Q_Value",
        palette="coolwarm_r",
        legend=True
    )
    
    for _, row in stats_df.head(5).iterrows():
        plt.text(row["Effect_Diff"] + 0.005, row["neg_log10_p"] + 0.05, row["Lineage"], fontsize=9)

    plt.axhline(y=-np.log10(0.05), color="gray", linestyle=":", label="p = 0.05 cutoff")
    plt.title("Lineage Differential PLK1 Dependency vs Pan-Cancer Cohort", fontsize=12, fontweight="bold")
    plt.xlabel("Delta Chronos Score (Lineage Mean - Pan-Cancer Mean)", fontsize=10)
    plt.ylabel("-log10(p-value)", fontsize=10)
    plt.tight_layout()

    volcano_path = os.path.join(FIGURES_DIR, "plk1_volcano_plot.png")
    plt.savefig(volcano_path, dpi=300)
    plt.close()
    print(f"[Saved Figure] Volcano plot saved to {volcano_path}")

    print("\n--- Top Lineages Ranked by PLK1 Essentiality ---")
    print(stats_df[["Lineage", "N_CellLines", "Mean_Chronos", "P_Value", "FDR_Q_Value"]].head(8).to_string(index=False))

if __name__ == "__main__":
    run_analysis()