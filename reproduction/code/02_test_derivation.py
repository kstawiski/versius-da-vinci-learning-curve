import pyarrow.parquet as pq
import pandas as pd
import numpy as np
import os

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df = pd.read_parquet(parquet_path)

# Filter Poland
df_pl = df[df["country"] == "Poland"].copy()
print(f"Total Polish rows: {len(df_pl)}")

# Platforms
print("Platform counts:")
print(df_pl["platform"].value_counts(dropna=False))

# Sorting columns null counts
sort_cols = ["surgery_date", "surgery_start_minutes_after_midnight", "operation_centre_within_day_order", "full_polish_series_case_sequence"]
for c in sort_cols:
    print(f"Nulls in {c}: {df_pl[c].isna().sum()}")

# Sort all 838 Polish operations:
# surgery_date, then surgery_start_minutes_after_midnight (missing last), then operation_centre_within_day_order, then full_polish_series_case_sequence
# In pandas:
df_pl = df_pl.sort_values(
    by=["surgery_date", "surgery_start_minutes_after_midnight", "operation_centre_within_day_order", "full_polish_series_case_sequence"],
    na_position="last"
).reset_index(drop=True)

# Check uniqueness of sequence
print("Is full_polish_series_case_sequence unique?", df_pl["full_polish_series_case_sequence"].is_unique)

# platform_n: running count within platform
df_pl["platform_n"] = df_pl.groupby("platform").cumcount() + 1

# dv_restart_n: da Vinci cases with platform_n >= 20
df_pl["dv_restart_n"] = np.where((df_pl["platform"] == "da Vinci") & (df_pl["platform_n"] >= 20), df_pl["platform_n"] - 19, np.nan)

# year = calendar year of surgery_date
df_pl["year"] = pd.to_datetime(df_pl["surgery_date"]).dt.year

# day_position = order within the surgeon's operating day (all platforms)
df_pl["day_position"] = df_pl.groupby("surgery_date").cumcount() + 1
df_pl["first_case"] = (df_pl["day_position"] == 1).astype(int)

# months_since_start = days since the first Polish operation / 30.4375
first_date = pd.to_datetime(df_pl["surgery_date"]).min()
df_pl["days_since_start"] = (pd.to_datetime(df_pl["surgery_date"]) - first_date).dt.days
df_pl["months_since_start"] = df_pl["days_since_start"] / 30.4375

# Data close = latest (surgery_date + fu_months_last_contact * 30.4375 days)
# Let's compute close date:
valid_fu = df_pl["fu_months_last_contact"].notna()
last_contact_dates = pd.to_datetime(df_pl["surgery_date"]) + pd.to_timedelta(df_pl["fu_months_last_contact"] * 30.4375, unit="D")
data_close = last_contact_dates.max()
print("Computed data close date exists:", pd.notna(data_close))

# Window-closed eligibility:
# 3-month continence: at least 120 days from surgery to data close
# 12-month continence: at least 425 days
# PSA persistence: at least 56 days
df_pl["days_to_data_close"] = (data_close - pd.to_datetime(df_pl["surgery_date"])).dt.days
df_pl["eligible_continence_3m"] = df_pl["days_to_data_close"] >= 120
df_pl["eligible_continence_12m"] = df_pl["days_to_data_close"] >= 425
df_pl["eligible_psa_persistence"] = df_pl["days_to_data_close"] >= 56

print("Window-closed eligibility counts:")
print("  3m continence eligible:", df_pl["eligible_continence_3m"].sum())
print("  12m continence eligible:", df_pl["eligible_continence_12m"].sum())
print("  psa persistence eligible:", df_pl["eligible_psa_persistence"].sum())

# Check R1 counts:
v_mask = df_pl["platform"] == "Versius"
dv_mask = df_pl["platform"] == "da Vinci"
dv_restart_mask = dv_mask & (df_pl["platform_n"] >= 20)
v_2024_mask = v_mask & (df_pl["year"] == 2024)
dv_restart_2024_mask = dv_restart_mask & (df_pl["year"] == 2024)

print("\nR1 Target Counts:")
print(f"Versius total: {v_mask.sum()} (target: 337)")
print(f"da Vinci total: {dv_mask.sum()} (target: 401)")
print(f"Restart series total: {dv_restart_mask.sum()} (target: 382)")
print(f"2024 Versius: {v_2024_mask.sum()} (target: 99)")
print(f"2024 restart da Vinci: {dv_restart_2024_mask.sum()} (target: 105)")
