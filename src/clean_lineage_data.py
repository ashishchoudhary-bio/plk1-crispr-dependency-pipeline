"""
src/clean_lineage_data.py
Robust lineage retrieval and merging with PLK1 CRISPR dependency scores.
"""
import os
import urllib.request
import pandas as pd

RAW_DIR = "data/raw"
PROCESSED_DIR = "data/processed"
os.makedirs(PROCESSED_DIR, exist_ok=True)

METADATA_URL = "https://ndownloader.figshare.com/files/34008503"

def clean_and_merge():
    plk1_path = os.path.join(RAW_DIR, "plk1_depmap_raw.csv")
    meta_path = os.path.join(RAW_DIR, "Model_Metadata.csv")

    if not os.path.exists(plk1_path):
        raise FileNotFoundError(f"Missing required file: {plk1_path}")

    plk1_df = pd.read_csv(plk1_path)
    print(f"[Loaded] PLK1 scores: {len(plk1_df)} cell lines.")

    # Agar metadata file missing ya choti hai to direct CDN mirror se layein
    if not os.path.exists(meta_path) or os.path.getsize(meta_path) < 10000:
        print("[Downloading] Fetching authentic Broad DepMap sample info...")
        headers = {"User-Agent": "Mozilla/5.0"}
        req = urllib.request.Request(METADATA_URL, headers=headers)
        with urllib.request.urlopen(req) as resp, open(meta_path, "wb") as f:
            f.write(resp.read())
        print(f"[Saved] Metadata table saved to {meta_path}")

    meta_df = pd.read_csv(meta_path, low_memory=False, on_bad_lines="skip")
    print(f"[Metadata Columns]: {list(meta_df.columns[:10])}")

    # ID column detection
    id_candidates = ["ModelID", "DepMap_ID", "Broad_ID", "CCLE_Name", "stripped_cell_line_name"]
    id_col = next((c for c in id_candidates if c in meta_df.columns), meta_df.columns[0])

    # Lineage column detection
    lineage_candidates = [
        "OncotreeLineage", "lineage", "Lineage", "primary_disease", 
        "PrimaryDisease", "cancer_type", "sample_collection_site"
    ]
    lineage_col = next((c for c in lineage_candidates if c in meta_df.columns), None)

    if not lineage_col:
        for c in meta_df.columns:
            if any(term in c.lower() for term in ["lineage", "disease", "site", "tissue"]):
                lineage_col = c
                break

    print(f"[Selected] Mapping '{id_col}' -> 'ModelID' and '{lineage_col}' -> 'Lineage'")

    if lineage_col:
        sub_meta = meta_df[[id_col, lineage_col]].copy()
        sub_meta.columns = ["ModelID", "Lineage"]
        sub_meta["Lineage"] = sub_meta["Lineage"].fillna("Other / Unknown").str.strip()
        
        merged = pd.merge(plk1_df, sub_meta, on="ModelID", how="left")
        merged["Lineage"] = merged["Lineage"].fillna("Other / Unknown")
    else:
        merged = plk1_df.copy()
        merged["Lineage"] = "Pan-Cancer Cohort"

    merged = merged.dropna(subset=["PLK1_Chronos_Score"])
    out_path = os.path.join(PROCESSED_DIR, "plk1_pan_cancer_cleaned.csv")
    merged.to_csv(out_path, index=False)

    print(f"\n[Success] Processed dataset saved to: {out_path}")
    print(f"[Total Lines]: {len(merged)} cell lines across {merged['Lineage'].nunique()} lineages.")
    print("\n--- Lineage Summary ---")
    print(merged["Lineage"].value_counts().head(12))

if __name__ == "__main__":
    clean_and_merge()