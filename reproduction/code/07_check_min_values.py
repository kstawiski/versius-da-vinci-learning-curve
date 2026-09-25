import pyarrow.parquet as pq
import pandas as pd
import numpy as np

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df = pd.read_parquet(parquet_path)
df_rob = df[(df["country"] == "Poland") & (df["platform"].isin(["Versius", "da Vinci"]))].copy()

print("psa_preop min:", df_rob["psa_preop"].min(), "<=0 count:", (df_rob["psa_preop"] <= 0).sum())
print("weight min:", df_rob["prostate_specimen_weight_g"].min(), "<=0 count:", (df_rob["prostate_specimen_weight_g"] <= 0).sum())
