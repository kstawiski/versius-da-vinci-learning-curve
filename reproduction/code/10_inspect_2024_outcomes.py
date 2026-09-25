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
df_pl["is_restart_dv"] = (df_pl["platform"] == "da Vinci") & (df_pl["platform_n"] >= 20)

# 2024 cohort
c2024 = df_pl[df_pl["year"] == 2024].copy()
v2024 = c2024[c2024["platform"] == "Versius"]
dv2024 = c2024[c2024["is_restart_dv"]]

print(f"2024 Versius: {len(v2024)}")
print(f"2024 restart da Vinci: {len(dv2024)}")

# Data close
last_contact_dates = pd.to_datetime(df_pl["surgery_date"]) + pd.to_timedelta(df_pl["fu_months_last_contact"] * 30.4375, unit="D")
data_close = last_contact_dates.max()
c2024["days_to_close"] = (data_close - pd.to_datetime(c2024["surgery_date"])).dt.days
c2024["elig_3m"] = c2024["days_to_close"] >= 120
c2024["elig_12m"] = c2024["days_to_close"] >= 425
c2024["elig_psa"] = c2024["days_to_close"] >= 56

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

c2024["pt_group"] = c2024["pT_final"].apply(get_pt_group)
c2024["psm"] = c2024["margin_status"].map({"positive": 1, "negative": 0})
c2024["los"] = pd.to_numeric(c2024["postoperative_length_of_stay_days"], errors="coerce")
c2024["readmission_30d_num"] = c2024["readmission_30d"].map({True: 1, False: 0})
c2024["psa_persistence_eau_num"] = c2024["psa_persistence_eau"].map({True: 1, False: 0})

def parse_pad_free(val):
    if pd.isna(val): return np.nan
    s = str(val).strip().lower()
    if s in ["true", "1", "1.0"]: return 1
    if s in ["false", "0", "0.0"]: return 0
    return np.nan

c2024["pad_free_3m_num"] = c2024["pad_free_3m"].apply(parse_pad_free)
c2024["pad_free_12m_num"] = c2024["pad_free_12m"].apply(parse_pad_free)

for grp_name, grp_df in [("Versius 2024", c2024[c2024["platform"] == "Versius"]),
                         ("da Vinci restart 2024", c2024[c2024["is_restart_dv"]])]:
    print(f"\n=== {grp_name} (N={len(grp_df)}) ===")
    print("or_time_min: mean =", grp_df["or_time_min"].mean(), "n_obs =", grp_df["or_time_min"].notna().sum())
    
    pt2_df = grp_df[grp_df["pt_group"] == "pT2"]
    print("pT2 psm: pos =", (pt2_df["psm"] == 1).sum(), "neg =", (pt2_df["psm"] == 0).sum(), "missing =", pt2_df["psm"].isna().sum(), "total pT2 =", len(pt2_df))
    
    print("All-stage psm: pos =", (grp_df["psm"] == 1).sum(), "neg =", (grp_df["psm"] == 0).sum(), "missing =", grp_df["psm"].isna().sum())
    print("LOS: mean =", grp_df["los"].mean(), "n_obs =", grp_df["los"].notna().sum())
    print("30d readmission: events =", (grp_df["readmission_30d_num"] == 1).sum(), "non-events =", (grp_df["readmission_30d_num"] == 0).sum(), "missing =", grp_df["readmission_30d_num"].isna().sum())
    
    psa_sub = grp_df[grp_df["elig_psa"]]
    print("PSA persistence (window-closed): events =", (psa_sub["psa_persistence_eau_num"] == 1).sum(), "non-events =", (psa_sub["psa_persistence_eau_num"] == 0).sum(), "missing =", psa_sub["psa_persistence_eau_num"].isna().sum(), "eligible =", len(psa_sub))
    
    pf3_sub = grp_df[grp_df["elig_3m"]]
    print("Pad-free 3m (window-closed): events =", (pf3_sub["pad_free_3m_num"] == 1).sum(), "non-events =", (pf3_sub["pad_free_3m_num"] == 0).sum(), "missing =", pf3_sub["pad_free_3m_num"].isna().sum(), "eligible =", len(pf3_sub))
    
    pf12_sub = grp_df[grp_df["elig_12m"]]
    print("Pad-free 12m (window-closed): events =", (pf12_sub["pad_free_12m_num"] == 1).sum(), "non-events =", (pf12_sub["pad_free_12m_num"] == 0).sum(), "missing =", pf12_sub["pad_free_12m_num"].isna().sum(), "eligible =", len(pf12_sub))
