#!/usr/bin/env python3
"""Stage 1 derivation for manuscript6 (single-surgeon Versius and da Vinci learning curves).

Reads the locked release (read-only), rebuilds the surgeon's case order from verified
surgery dates and operating-room start times, derives exposures, outcomes and covariates,
and writes:
  analysis/restricted/ms6_analysis_dataset.parquet   (patient-level, mode 0600)
  evidence/DERIVATION_SUMMARY.json                   (aggregate only)
Prints aggregate counts only. Never prints identifiers, dates or free text.
"""
import json, os, re, sys, hashlib
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
ROOT = WS.parent
REL = ROOT / 'data/prostatectomy_analysis_database_v1.1_2026-09-04'
RESTRICTED = WS / 'analysis/restricted'
OUT_DS = RESTRICTED / 'ms6_analysis_dataset.parquet'
OUT_SUM = WS / 'evidence/DERIVATION_SUMMARY.json'

src = REL / 'parquet/patient_wide.parquet'
src_sha = hashlib.sha256(src.read_bytes()).hexdigest()
pw = pd.read_parquet(src)
pl = pw[pw['country'].astype(str) == 'Poland'].copy()
assert len(pl) == 838, len(pl)

def tbool(s):
    """Tolerant boolean parse: True/1/1.0 -> 1, False/0/0.0 -> 0, else NaN. Blank stays missing."""
    m = {'true': 1.0, '1': 1.0, '1.0': 1.0, 'false': 0.0, '0': 0.0, '0.0': 0.0}
    return s.astype('string').str.strip().str.lower().map(m).astype(float)

def num(s):
    return pd.to_numeric(s, errors='coerce')

d = pd.DataFrame(index=pl.index)
d['pid'] = pl['master_patient_id'].astype(str)
d['platform'] = pl['platform'].astype(str)
d['centre'] = pl['operation_centre'].map({'Salve Medica': 'SM', 'Bełchatów': 'BE'})
assert d['centre'].notna().all()
d['surgery_date'] = pd.to_datetime(pl['surgery_date'], errors='coerce')
assert d['surgery_date'].notna().all()
d['start_min'] = num(pl['surgery_start_minutes_after_midnight'])
d['within_day_centre_order'] = num(pl['operation_centre_within_day_order'])
d['released_series_seq'] = num(pl['full_polish_series_case_sequence'])
d['date_corrected'] = (pl['surgery_date_decision'].astype(str) == 'corrected_against_primary_registry').astype(int)
d['within_day_precision'] = pl['within_day_order_precision'].astype(str)

# ---- Rebuilt surgeon case order ----
o = d.sort_values(['surgery_date', 'start_min', 'within_day_centre_order', 'released_series_seq'],
                  na_position='last', kind='mergesort').copy()
o['series_n'] = np.arange(1, len(o) + 1)
o['robotic'] = o['platform'].isin(['Versius', 'da Vinci']).astype(int)
o['robotic_n'] = np.where(o['robotic'] == 1, o['robotic'].cumsum(), np.nan)
for p, tag in [('Versius', 'v'), ('da Vinci', 'dv'), ('Laparoscopic', 'lap')]:
    is_p = (o['platform'] == p).astype(int)
    o[f'prior_{tag}'] = is_p.cumsum() - is_p
o['platform_n'] = o.groupby('platform').cumcount() + 1
o['team_platform_n'] = o.groupby(['centre', 'platform']).cumcount() + 1
o['prior_robotic'] = o['prior_v'] + o['prior_dv']
# Same-day workload and position (surgeon level, all platforms)
o['same_day_cases'] = o.groupby('surgery_date')['series_n'].transform('size')
o['day_position'] = o.groupby('surgery_date').cumcount() + 1
# Unresolved ties: operations on the same date whose relative order is not fixed by distinct start times.
# A day with any missing start time among two or more operations is one block; otherwise operations
# sharing an identical start time form a block.
day_n = o.groupby('surgery_date')['series_n'].transform('size')
day_any_missing = o.groupby('surgery_date')['start_min'].transform(lambda x: x.isna().any())
dup_start = o.duplicated(['surgery_date', 'start_min'], keep=False) & o['start_min'].notna()
tie = (day_n > 1) & (day_any_missing | dup_start)
o['order_tie_unresolved'] = tie.astype(int)
o['tie_block'] = np.where(tie, o.groupby(['surgery_date']).ngroup() * 10000 + np.where(day_any_missing, 0, o['start_min'].fillna(0)), -1).astype(int)
# Calendar
t0 = o['surgery_date'].min()
o['months_since_start'] = (o['surgery_date'] - t0).dt.days / 30.4375
o['year'] = o['surgery_date'].dt.year
# da Vinci phases relative to the Versius programme
first_v = o.loc[o['platform'] == 'Versius', 'series_n'].min()
last_v = o.loc[o['platform'] == 'Versius', 'series_n'].max()
dvp = np.select([o['series_n'] < first_v, o['series_n'] < last_v], ['pre_versius', 'during_versius'], 'after_versius')
o['dv_phase'] = np.where(o['platform'] == 'da Vinci', dvp, None)
# Overlap window: from first da Vinci case after Versius case 300 ... pragmatic definition = calendar 2024
o['overlap_2024'] = (o['year'] == 2024).astype(int)

d = o

# ---- Covariates ----
d['age'] = num(pl['age'])
d['bmi'] = num(pl['bmi'])
d['psa'] = num(pl['psa_preop'])
d['biopsy_isup'] = num(pl['biopsy_isup'])
d['ct'] = pl['clinical_t_norm'].astype('string').replace('', pd.NA)
d['eau_risk'] = pl['eau_risk'].astype('string').replace('', pd.NA)
d['capra_s'] = num(pl['capra_s'])
w = num(pl['prostate_specimen_weight_g'])
d['specimen_weight_g'] = w
pt = pl['pT_final'].astype('string').str.strip().replace('', pd.NA)
d['pT'] = pt
d['pT_group'] = pd.Series(np.select([pt.str.contains(r'T2', na=False), pt.str.contains(r'T3a', na=False),
                           pt.str.contains(r'T3b|T4', na=False)], ['pT2', 'pT3a', 'pT3b+'], None), index=pt.index)
d['pT3plus'] = np.where(d['pT_group'].isna(), np.nan, (d['pT_group'] != 'pT2').astype(float))
d['isup_final'] = num(pl['isup_final'])
d['node_positive'] = tbool(pl['node_positive'])
d['nerve_sparing'] = tbool(pl['nerve_sparing_performed'])
d['plnd'] = pl['pelvic_lymph_node_dissection_state'].map({'performed': 1.0, 'not_performed': 0.0})
d['prior_turp'] = tbool(pl['prior_turp'])
d['neoadjuvant_hormonal'] = pl['neoadjuvant_hormonal_state'].map({'documented_present': 1.0, 'documented_absent': 0.0})

# ---- Outcomes ----
d['or_time'] = num(pl['or_time_min'])
d['or_time_source'] = pl['or_time_source'].astype(str)
d['or_time_clock_concordant'] = tbool(pl['or_time_min_clock_concordant'])
d['or_time_registry_conflict'] = tbool(pl['or_time_registry_conflict'])
ebl = num(pl['blood_loss_ml'])
impl = tbool(pl['blood_loss_ml_implausible_repeated_value']).fillna(0)
d['ebl_implausible_flag'] = impl
d['ebl'] = ebl.where(impl != 1)
d['los'] = num(pl['postoperative_length_of_stay_days'])
d['catheter_days'] = num(pl['catheter_days'])

cd_rank = {'none': 0, 'I': 1, 'II': 2, 'III': 3, 'IIIa': 3, 'IIIb': 3, 'IV': 4, 'IVa': 4, 'IVb': 4, 'V': 5}
def cd(col):
    s = pl[col].astype('string').str.strip()
    return s.map(cd_rank).astype(float)
d['cd30'] = cd('clavien_max_30d_src')
d['cd90'] = cd('clavien_max_90d_src')
d['cd_any'] = cd('clavien_dindo_highest_documented')
d['major_cd30'] = np.where(d['cd30'].isna(), np.nan, (d['cd30'] >= 3).astype(float))
d['major_cd90'] = np.where(d['cd90'].isna(), np.nan, (d['cd90'] >= 3).astype(float))
d['major_cd_any'] = np.where(d['cd_any'].isna(), np.nan, (d['cd_any'] >= 3).astype(float))
d['any_cd30'] = np.where(d['cd30'].isna(), np.nan, (d['cd30'] >= 1).astype(float))
d['cd30_ge2'] = np.where(d['cd30'].isna(), np.nan, (d['cd30'] >= 2).astype(float))
d['reoperation'] = pl['reoperation_registry'].map({'performed': 1.0, 'not_performed': 0.0})
d['urine_leak'] = pl['urine_leak_registry'].map({'documented_present': 1.0, 'documented_absent': 0.0})
d['readmission_30d'] = tbool(pl['readmission_30d'])
d['readmission_90d'] = tbool(pl['readmission_90d'])
d['transfusion_30d'] = tbool(pl['transfusion_30d'])

d['psm'] = pl['margin_status'].map({'positive': 1.0, 'negative': 0.0})
d['margin_r_discordant'] = (tbool(pl['margin_r_status_concordant']) == 0).astype(int)
d['psm_pt2'] = np.where(d['pT_group'] == 'pT2', d['psm'], np.nan)
d['psm_pt3'] = np.where(d['pT_group'].isin(['pT3a', 'pT3b+']), d['psm'], np.nan)
d['psm_length_mm'] = num(pl['positive_margin_length_mm'])
d['psm_length_relation'] = pl['margin_length_relation'].astype(str)
for site, col in [('apex', 'margin_at_apex'), ('base', 'margin_at_base'), ('posterior', 'margin_at_posterior'),
                  ('bladder_neck', 'margin_at_bladder_neck'), ('nvb', 'margin_at_neurovascular_bundle_any')]:
    d[f'psm_site_{site}'] = tbool(pl[col])
d['ln_removed'] = num(pl['pelvic_lymph_nodes_removed'])
d['ln_yield'] = np.where(d['plnd'] == 1, d['ln_removed'], np.nan)

for m in ['3m', '6m', '12m']:
    d[f'pad_free_{m}'] = tbool(pl[f'pad_free_{m}'])
    d[f'pads_{m}'] = num(pl[f'pads_{m}'])
    dayc = f'pads_{m}_days_postoperative'
    d[f'pads_{m}_day'] = num(pl[dayc]) if dayc in pl.columns else np.nan
d['psa_persistence_eau'] = tbool(pl['psa_persistence_eau'])
d['psa_persistence_nccn'] = tbool(pl['psa_persistence_nccn'])
d['bcr_nccn'] = tbool(pl['bcr_nccn'])
d['followup_assessable'] = tbool(pl['followup_assessable'])
d['fu_months_contact'] = num(pl['fu_months_last_contact'])
d['fu_months_psa'] = num(pl['fu_months_last_psa'])
d['fu_last_contact_date'] = pd.to_datetime(pl['fu_last_contact_date'], errors='coerce')

# ---- Source reconciliation assertions (index alignment and counts) ----
chk = pl.join(d[['pid', 'pT_group', 'psm', 'or_time', 'platform', 'centre']].rename(columns=lambda c: 'd_' + c))
assert (chk['d_pid'] == chk['master_patient_id'].astype(str)).all()
assert (chk['d_platform'] == chk['platform'].astype(str)).all()
src_pt3 = chk['pT_final'].astype(str).str.contains('T3|T4').sum()
assert int((chk['d_pT_group'].isin(['pT3a', 'pT3b+'])).sum()) == int(src_pt3), 'pT group misaligned'
assert int((chk['d_psm'] == 1).sum()) == int((chk['margin_status'].astype(str) == 'positive').sum())
assert np.allclose(chk['d_or_time'].fillna(-1), pd.to_numeric(chk['or_time_min'], errors='coerce').fillna(-1))
assert d.index.equals(o.index)

# ---- Write restricted dataset ----
RESTRICTED.mkdir(parents=True, exist_ok=True)
os.chmod(RESTRICTED, 0o700)
d = d.sort_values('series_n').reset_index(drop=True)
d.to_parquet(OUT_DS, index=False)
os.chmod(OUT_DS, 0o600)

# ---- Aggregate summary ----
rob = d[d['robotic'] == 1]
s = {'source_patient_wide_sha256': src_sha, 'n_poland': int(len(d)),
     'platform_counts': d['platform'].value_counts().to_dict(),
     'robotic_n': int(len(rob)),
     'order_rebuild': {
         'released_vs_rebuilt_positions_differ': int((d['released_series_seq'] != d['series_n']).sum()),
         'max_abs_position_shift': int((d['released_series_seq'] - d['series_n']).abs().max()),
         'date_corrected_cases': int(d['date_corrected'].sum()),
         'unresolved_same_day_ties': int(d['order_tie_unresolved'].sum()),
         'unresolved_ties_robotic': int(d.loc[d.robotic == 1, 'order_tie_unresolved'].sum()),
         'tie_blocks': int(d.loc[d.tie_block >= 0, 'tie_block'].nunique()),
         'within_day_precision': d['within_day_precision'].value_counts().to_dict()},
     'exposure_structure': {
         'lap_before_first_versius': int(d.loc[d['platform'] == 'Versius', 'prior_lap'].min()),
         'dv_before_first_versius': int(d.loc[d['platform'] == 'Versius', 'prior_dv'].min()),
         'dv_phase_counts': rob.loc[rob['platform'] == 'da Vinci', 'dv_phase'].value_counts().to_dict(),
         'versius_cases_by_prior_dv_band': pd.cut(rob.loc[rob.platform == 'Versius', 'prior_dv'], [-1, 0, 13, 46, 100, 200, 401]).astype(str).value_counts().sort_index().to_dict(),
         'dv_cases_by_prior_versius_band': pd.cut(rob.loc[rob.platform == 'da Vinci', 'prior_v'], [-1, 0, 100, 200, 300, 336, 337]).astype(str).value_counts().sort_index().to_dict(),
         'overlap_2024_by_platform': rob[rob['overlap_2024'] == 1]['platform'].value_counts().to_dict(),
         'team_platform_max': {f'{a}|{b}': int(v) for (a, b), v in d.groupby(['centre', 'platform'])['team_platform_n'].max().items()}},
     'populated_by_platform': {c: rob.groupby('platform')[c].apply(lambda x: int(x.notna().sum())).to_dict()
                               for c in ['or_time', 'ebl', 'los', 'cd30', 'cd90', 'cd_any', 'reoperation', 'urine_leak',
                                         'readmission_30d', 'readmission_90d', 'psm', 'psm_pt2', 'psm_length_mm', 'ln_yield',
                                         'pad_free_3m', 'pad_free_6m', 'pad_free_12m', 'psa_persistence_eau',
                                         'psa_persistence_nccn', 'bcr_nccn', 'bmi', 'psa', 'biopsy_isup', 'ct',
                                         'specimen_weight_g', 'pT_group', 'isup_final', 'nerve_sparing', 'plnd']},
     'flags': {'ebl_implausible_2400_set_missing': int(d['ebl_implausible_flag'].sum()),
               'or_time_not_clock_source': d.loc[d.robotic == 1, 'or_time_source'].value_counts().to_dict(),
               'or_time_clock_discordant': int((rob['or_time_clock_concordant'] == 0).sum()),
               'margin_r_discordant': int(rob['margin_r_discordant'].sum())}}
OUT_SUM.write_text(json.dumps(s, indent=1, default=str))
txt = OUT_SUM.read_text()
if re.search(r'\b(19|20)\d\d-\d\d-\d\d\b', txt) or re.search(r'[0-9a-f]{8}-[0-9a-f]{4}-', txt):
    OUT_SUM.unlink(); sys.exit('privacy scan failed; summary removed')
print(json.dumps({k: s[k] for k in ['n_poland', 'robotic_n', 'order_rebuild', 'exposure_structure', 'flags']}, indent=1, default=str))
