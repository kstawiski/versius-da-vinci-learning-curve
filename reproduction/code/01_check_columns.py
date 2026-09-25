import pyarrow.parquet as pq

parquet_path = "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
schema = pq.read_schema(parquet_path)
names = set(schema.names)

needed = [
    "country", "platform", "surgery_date", "surgery_start_minutes_after_midnight",
    "operation_centre_within_day_order", "full_polish_series_case_sequence",
    "fu_months_last_contact", "or_time_min", "margin_status", "pT_final",
    "postoperative_length_of_stay_days", "blood_loss_ml", "blood_loss_ml_implausible_repeated_value",
    "readmission_30d", "psa_persistence_eau", "pad_free_3m", "pad_free_12m",
    "clavien_dindo_highest_documented", "age", "bmi", "psa_preop", "biopsy_isup",
    "prostate_specimen_weight_g", "nerve_sparing_performed", "pelvic_lymph_node_dissection_state",
    "operation_centre"
]

print("Total columns in parquet:", len(schema.names))
missing = [c for c in needed if c not in names]
print("Missing columns:", missing)
for c in needed:
    if c in names:
        idx = schema.get_field_index(c)
        print(f"  {c}: {schema.field(idx).type}")
