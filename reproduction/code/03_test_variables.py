import pyarrow.parquet as pq
import pandas as pd
import numpy as np

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df = pd.read_parquet(parquet_path)
df_pl = df[df["country"] == "Poland"].copy()

print("--- margin_status ---")
print(df_pl["margin_status"].value_counts(dropna=False))

print("\n--- pT_final ---")
print(df_pl["pT_final"].value_counts(dropna=False))

print("\n--- postoperative_length_of_stay_days ---")
print(df_pl["postoperative_length_of_stay_days"].value_counts(dropna=False).head(10))

print("\n--- blood_loss_ml_implausible_repeated_value ---")
print(df_pl["blood_loss_ml_implausible_repeated_value"].value_counts(dropna=False))

print("\n--- readmission_30d ---")
print(df_pl["readmission_30d"].value_counts(dropna=False))

print("\n--- psa_persistence_eau ---")
print(df_pl["psa_persistence_eau"].value_counts(dropna=False))

print("\n--- pad_free_3m ---")
print(df_pl["pad_free_3m"].value_counts(dropna=False))

print("\n--- pad_free_12m ---")
print(df_pl["pad_free_12m"].value_counts(dropna=False))

print("\n--- clavien_dindo_highest_documented ---")
print(df_pl["clavien_dindo_highest_documented"].value_counts(dropna=False))

print("\n--- biopsy_isup ---")
print(df_pl["biopsy_isup"].value_counts(dropna=False))

print("\n--- nerve_sparing_performed ---")
print(df_pl["nerve_sparing_performed"].value_counts(dropna=False))

print("\n--- pelvic_lymph_node_dissection_state ---")
print(df_pl["pelvic_lymph_node_dissection_state"].value_counts(dropna=False))

print("\n--- operation_centre ---")
print(df_pl["operation_centre"].value_counts(dropna=False))

print("\n--- or_time_min nulls ---")
print("Null or_time_min:", df_pl["or_time_min"].isna().sum())
