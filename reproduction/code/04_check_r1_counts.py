import pyarrow.parquet as pq
import pandas as pd
import numpy as np

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df = pd.read_parquet(parquet_path)
df_pl = df[df["country"] == "Poland"].copy()

# Sort
df_pl = df_pl.sort_values(
    by=["surgery_date", "surgery_start_minutes_after_midnight", "operation_centre_within_day_order", "full_polish_series_case_sequence"],
    na_position="last"
).reset_index(drop=True)

df_pl["platform_n"] = df_pl.groupby("platform").cumcount() + 1
df_pl["year"] = pd.to_datetime(df_pl["surgery_date"]).dt.year

# pT group
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

df_pl["pt_group"] = df_pl["pT_final"].apply(get_pt_group)
df_pl["is_pt2"] = (df_pl["pt_group"] == "pT2")
df_pl["has_or_time"] = df_pl["or_time_min"].notna()

# Series definitions
series_dict = {
    "Versius": df_pl["platform"] == "Versius",
    "da Vinci": df_pl["platform"] == "da Vinci",
    "Restart da Vinci": (df_pl["platform"] == "da Vinci") & (df_pl["platform_n"] >= 20),
    "2024 Versius": (df_pl["platform"] == "Versius") & (df_pl["year"] == 2024),
    "2024 restart da Vinci": (df_pl["platform"] == "da Vinci") & (df_pl["platform_n"] >= 20) & (df_pl["year"] == 2024)
}

print("=== R1 Counts Summary ===")
for name, mask in series_dict.items():
    n_total = mask.sum()
    n_pt2 = (mask & df_pl["is_pt2"]).sum()
    n_or_time = (mask & df_pl["has_or_time"]).sum()
    print(f"{name}:")
    print(f"  Total cases: {n_total}")
    print(f"  pT2 cases: {n_pt2} (denominator: {n_total})")
    print(f"  Cases with or_time: {n_or_time} (denominator: {n_total})")
