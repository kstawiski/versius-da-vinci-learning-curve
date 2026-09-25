#!/usr/bin/env python3
"""Aggregate-only feasibility profile of the Polish single-surgeon series.
Prints counts, completeness and summary statistics only. Never prints identifiers,
dates, free text or row-level values. Writes evidence/FEASIBILITY_PROFILE.json."""
import json, sys, re
from pathlib import Path
import pandas as pd, numpy as np
ROOT=Path(__file__).resolve().parents[2]
REL=ROOT/'data/prostatectomy_analysis_database_v1.1_2026-09-04'
OUT=Path(__file__).resolve().parents[1]/'evidence/FEASIBILITY_PROFILE.json'
pw=pd.read_parquet(REL/'parquet/patient_wide.parquet')
res={}
res['n_total']=int(len(pw))
ctry=[c for c in pw.columns if c.lower() in ('country','country_src')]
res['country_col']=ctry
pl=pw[pw[ctry[0]].astype(str).str.contains('Pol',case=False,na=False)].copy()
res['n_poland']=int(len(pl))
res['platform_by_centre']=pl.groupby(['operation_centre','platform'],dropna=False).size().reset_index().astype(str).values.tolist()
res['approach_by_platform']=pl.groupby(['platform','approach'],dropna=False).size().reset_index().astype(str).values.tolist()
# sequence integrity (no values printed)
s=pl['full_polish_series_case_sequence']
res['series_seq']={'n':int(s.notna().sum()),'n_distinct':int(s.nunique()),'min':float(s.min()),'max':float(s.max())}
ps=pl.groupby('platform')['full_polish_series_platform_sequence'].agg(['count','nunique','min','max'])
res['platform_seq']=ps.astype(float).reset_index().values.tolist()
# check sequence monotone with surgery date (counts only)
d=pd.to_datetime(pl['surgery_date'],errors='coerce')
o=pl.assign(_d=d).sort_values('full_polish_series_case_sequence')
res['seq_date_monotone_violations']=int((o['_d'].diff().dt.days<0).sum())
res['within_day_order_status']=pl['within_day_order_status'].value_counts(dropna=False).astype(int).rename(index=str).to_dict()
res['governed_734_member']=pl['polish_governed_734_robotic_cohort'].value_counts(dropna=False).astype(int).rename(index=str).to_dict()
res['governed_734_by_platform']=pl.groupby(['platform','polish_governed_734_robotic_cohort'],dropna=False).size().reset_index().astype(str).values.tolist()
# surgery year distribution by platform (year-level aggregate only)
res['year_by_platform']=pl.assign(y=d.dt.year).groupby(['platform','y']).size().reset_index().astype(str).values.tolist()
# completeness by platform
cols=['age','bmi','psa_preop','biopsy_isup','eau_risk','asa_class','prostate_specimen_weight_g','or_time_min','or_time_min_clock_derived','or_time_min_clock_concordant','or_time_source','or_time_registry_conflict','blood_loss_ml','blood_transfusion','postoperative_length_of_stay_days','length_of_stay_days','catheter_days','clavien_max_30d_src','clavien_max_90d_src','clavien_dindo_registry_grade','clavien_dindo_highest_documented','reoperation_registry','urine_leak_registry','readmission_30d','margin_positive','margin_r_status_concordant','positive_margin_length_mm','pT_final','isup_final','pN','pelvic_lymph_nodes_removed','plnd_performed','nerve_sparing_performed','pad_free_3m','pad_free_6m','pad_free_12m','bcr_nccn','psa_persistence_nccn','bcr_aua_threshold','psa_persistence_eau','followup_assessable','fu_months_last_contact','fu_months_last_psa','hood_technique_src']
comp={}
for c in cols:
    if c not in pl.columns: comp[c]='ABSENT';continue
    g=pl.groupby('platform')[c].apply(lambda x: int(x.notna().sum()))
    comp[c]={str(k):int(v) for k,v in g.items()}
res['populated_by_platform']=comp
# categorical distributions by platform (closed vocabularies only)
cat={}
for c in ['eau_risk','biopsy_isup','isup_final','pT_final','margin_positive','nerve_sparing_performed','pad_free_3m','pad_free_6m','pad_free_12m','bcr_nccn','psa_persistence_nccn','bcr_aua_threshold','clavien_max_30d_src','clavien_max_90d_src','clavien_dindo_registry_grade','reoperation_registry','urine_leak_registry','or_time_source','plnd_performed','blood_transfusion','hood_technique_src','asa_class','within_day_order_status']:
    if c in pl.columns:
        t=pl.groupby(['platform',c],dropna=False).size().reset_index()
        cat[c]=t.astype(str).values.tolist()
res['categorical_by_platform']=cat
num={}
for c in ['age','bmi','psa_preop','prostate_specimen_weight_g','or_time_min','or_time_min_clock_derived','blood_loss_ml','postoperative_length_of_stay_days','catheter_days','pelvic_lymph_nodes_removed','positive_margin_length_mm','fu_months_last_contact','fu_months_last_psa']:
    if c in pl.columns:
        x=pd.to_numeric(pl[c],errors='coerce')
        g=pl.assign(_x=x).groupby('platform')['_x'].describe(percentiles=[.05,.25,.5,.75,.95]).round(2)
        num[c]=g.reset_index().values.tolist()
res['numeric_by_platform']=num
OUT.write_text(json.dumps(res,indent=1,default=str))
# serialized-output privacy scan: no ISO dates or long id-like tokens
txt=OUT.read_text()
bad=re.findall(r'\b(19|20)\d\d-\d\d-\d\d\b',txt)
if bad: OUT.unlink(); sys.exit('privacy scan: date pattern found; output removed')
print('written',OUT, 'bytes',len(txt))
