import pyarrow.parquet as pq
import pandas as pd

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df = pd.read_parquet(parquet_path)
df_rob = df[(df["country"] == "Poland") & (df["platform"].isin(["Versius", "da Vinci"]))].copy()

print("Robotic psm counts:")
print(df_rob["margin_status"].value_counts(dropna=False))
