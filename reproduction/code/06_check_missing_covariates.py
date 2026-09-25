import pyarrow.parquet as pq
import pandas as pd
import numpy as np

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df = pd.read_parquet(parquet_path)
df_rob = df[(df["country"] == "Poland") & (df["platform"].isin(["Versius", "da Vinci"]))].copy()

def get_pt_group(val):
    if not isinstance(val, str) or not val.strip():
        return np.nan
    val = val.strip()
    if "T2" in val:
        return "pT2"
    elif "T3a" in val:
        return "pT3a"
    elif "T3b" in val or "T4" in val:
        return "pT3b+"
    return np.nan

df_rob["pt_group"] = df_rob["pT_final"].apply(get_pt_group)
df_rob["isup_cat"] = df_rob["biopsy_isup"].apply(lambda x: 4 if x in [4, 5] else (x if pd.notna(x) else np.nan))
df_rob["log_psa"] = np.log(df_rob["psa_preop"])
df_rob["log_weight"] = np.log(df_rob["prostate_specimen_weight_g"])
df_rob["plnd"] = df_rob["pelvic_lymph_node_dissection_state"].map({"performed": 1, "not_performed": 0})
df_rob["ns"] = df_rob["nerve_sparing_performed"].astype(float) # 1.0, 0.0, nan

cols = ["age", "bmi", "log_psa", "isup_cat", "log_weight", "pt_group", "ns", "plnd", "operation_centre"]
for c in cols:
    n_miss = df_rob[c].isna().sum()
    print(f"Missing in {c}: {n_miss} / {len(df_rob)}")

print("\nMissing by platform:")
for p in ["Versius", "da Vinci"]:
    sub = df_rob[df_rob["platform"] == p]
    print(f"-- Platform {p} (N={len(sub)}) --")
    for c in cols:
        print(f"  {c}: {sub[c].isna().sum()} missing")
