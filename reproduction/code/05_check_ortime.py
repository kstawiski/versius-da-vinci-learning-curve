import pyarrow.parquet as pq
import pandas as pd
import numpy as np

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df = pd.read_parquet(parquet_path)
df_rob = df[(df["country"] == "Poland") & (df["platform"].isin(["Versius", "da Vinci"]))].copy()

print("Total robotic cases:", len(df_rob))
print("or_time_min nulls:", df_rob["or_time_min"].isna().sum())
print("or_time_min summary statistics:")
print("Min:", df_rob["or_time_min"].min())
print("Max:", df_rob["or_time_min"].max())
print("Mean:", df_rob["or_time_min"].mean())
print("Median:", df_rob["or_time_min"].median())
print("Count by platform:")
print(df_rob.groupby("platform")["or_time_min"].count())
