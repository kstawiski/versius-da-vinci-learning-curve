#!/usr/bin/env python3
"""POST HOC endpoint provenance and cross-source agreement for the revision (aggregate output only).
Reads the locked release read-only and the restricted analysis dataset, and writes
  analysis/results/R_endpoint_verification.json
  tables/tableS_verification.md   (primary-endpoint sources and agreement by hospital and series)
  tables/tableS_missingness.md    (observed values by variable and series)
Covers the operative-time clock fields, pT stage, margin status against the pathology R status and the
structured margin-length field, the clock start time of the first operation of each day by hospital and
year, and the number of EAU PSA persistence classifications that depend on the rule that every value in
days 28-56 must reach 0.1 ng/mL. Never prints identifiers, dates of individual operations or free text."""
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
REL = WS.parent / 'data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet'
d = pd.read_parquet(WS / 'analysis/restricted/ms6_analysis_dataset.parquet')
r = d[d.robotic == 1].copy()
r['series'] = np.where(r.platform == 'Versius', 'Versius', np.where(r.platform_n < 20, 'da Vinci 1-19', 'da Vinci restart'))
cols = ['master_patient_id', 'surgery_date_decision', 'robot_model', 'platform_conflict', 'or_time_source', 'or_time_min_clock_concordant',
        'or_time_registry_conflict', 'or_start_time', 'pt_stage_registry_conflict', 'R_status', 'margin_r_status_concordant',
        'margin_length_relation', 'positive_margin_length_mm', 'margin_length_registry_conflict']
pw = pd.read_parquet(REL / 'patient_wide.parquet', columns=cols).rename(columns={'master_patient_id': 'pid'})
r = r.drop(columns=[c for c in pw.columns if c != 'pid' and c in r.columns]).merge(pw, on='pid', how='left', validate='one_to_one')
assert len(r) == 738

def tb(s):
    return s.astype('string').str.strip().str.lower().map({'true': 1.0, '1': 1.0, '1.0': 1.0, 'false': 0.0, '0': 0.0, '0.0': 0.0}).astype(float)
for c in ['or_time_min_clock_concordant', 'or_time_registry_conflict', 'pt_stage_registry_conflict', 'margin_r_status_concordant',
          'margin_length_registry_conflict', 'platform_conflict']:
    r[c] = tb(r[c])
rel = r.margin_length_relation.astype('string').fillna('').str.strip()
r['mlen_field'] = (rel != '').astype(float)

def frac(num, den):
    return f"{int(num)}/{int(den)} ({100 * num / den:.0f}%)" if den else '–'

strata = [('SalveMedica', 'SM', 'Versius'), ('SalveMedica', 'SM', 'da Vinci restart'), ('Bełchatów', 'BE', 'Versius'), ('Bełchatów', 'BE', 'da Vinci restart')]
res = {'label': 'POST HOC endpoint provenance and cross-source agreement', 'strata': {}}
rows = {}
for hn, hc, sr in strata:
    g = r[(r.centre == hc) & (r.series == sr)]
    pos, neg = g[g.psm == 1], g[g.psm == 0]
    clock = g[g.or_time_source == 'registry_start_end_clock']
    s = {
        'operations': len(g),
        'date_corrected_against_registry': int((g.surgery_date_decision.astype(str) == 'corrected_against_primary_registry').sum()),
        'robot_model_generation_named': int(g.robot_model.astype(str).isin(['da Vinci X', 'da Vinci Xi', 'da Vinci Si', 'da Vinci SP']).sum()),
        'robot_model_generation_values': g.robot_model.astype(str).value_counts().to_dict(),
        'platform_conflict_flag': int((g.platform_conflict == 1).sum()),
        'or_time_recorded': int(g.or_time.notna().sum()),
        'or_time_from_start_end_clock': len(clock),
        'or_time_duration_only': int((g.or_time_source == 'registry_duration_only').sum()),
        'clock_concordant': int((g.or_time_min_clock_concordant == 1).sum()), 'clock_concordance_assessed': int(g.or_time_min_clock_concordant.notna().sum()),
        'or_time_registry_conflict': int((g.or_time_registry_conflict == 1).sum()),
        'pT_recorded': int(g.pT_group.notna().sum()), 'pT_registry_conflict': int((g.pt_stage_registry_conflict == 1).sum()),
        'margin_recorded': int(g.psm.notna().sum()), 'margin_positive': len(pos), 'margin_negative': len(neg),
        'R_concordant_positive': int((pos.margin_r_status_concordant == 1).sum()), 'R_assessed_positive': int(pos.margin_r_status_concordant.notna().sum()),
        'R_concordant_negative': int((neg.margin_r_status_concordant == 1).sum()), 'R_assessed_negative': int(neg.margin_r_status_concordant.notna().sum()),
        'positive_with_length_field': int((pos.mlen_field == 1).sum()),
        'positive_with_length_value': int(pos.positive_margin_length_mm.notna().sum()),
        'negative_with_length_field': int((neg.mlen_field == 1).sum()),
        'length_field_present': int((g.mlen_field == 1).sum()),
        'margin_length_registry_conflict': int((g.margin_length_registry_conflict == 1).sum()),
    }
    res['strata'][f'{hn} {sr}'] = s
    rows[f'{hn}, {sr}'] = {
        'Operations': str(s['operations']),
        'Surgery date corrected against the operative registry': str(s['date_corrected_against_registry']),
        'Operative time recorded': frac(s['or_time_recorded'], s['operations']),
        'Operative time from registry start and end clock times': frac(s['or_time_from_start_end_clock'], s['operations']),
        'Recorded duration within 15 min of the clock difference': frac(s['clock_concordant'], s['clock_concordance_assessed']),
        'Operative time differing between registry files': str(s['or_time_registry_conflict']),
        'pT stage recorded': frac(s['pT_recorded'], s['operations']),
        'pT stage differing between registry files': str(s['pT_registry_conflict']),
        'Margin status recorded': frac(s['margin_recorded'], s['operations']),
        'Positive margins with R1 in the pathology R field': frac(s['R_concordant_positive'], s['R_assessed_positive']),
        'Negative margins with R0 in the pathology R field': frac(s['R_concordant_negative'], s['R_assessed_negative']),
        'Positive margins with a margin-length entry': frac(s['positive_with_length_field'], s['margin_positive']),
        'Negative margins with a margin-length entry': frac(s['negative_with_length_field'], s['margin_negative']),
    }
T = pd.DataFrame(rows)
md = '| Check | ' + ' | '.join(T.columns) + ' |\n|' + '|'.join('---' for _ in range(len(T.columns) + 1)) + '|\n'
md += '\n'.join('| ' + idx + ' | ' + ' | '.join(T.loc[idx]) + ' |' for idx in T.index)
(WS / 'tables/tableS_verification.md').write_text(md + '\n')

# registry margin-length entry by hospital, robot, period and margin status
r['period3'] = np.select([r.year <= 2023, r.year == 2024], ['2021-2023', '2024'], '2025-2026')
mb = {}
for (hc, sr, per), g in r[r.series != 'da Vinci 1-19'].groupby(['centre', 'series', 'period3']):
    pos, neg = g[g.psm == 1], g[g.psm == 0]
    mb[f'{hc} {sr} {per}'] = {'positive_with_entry': [int((pos.mlen_field == 1).sum()), len(pos)], 'negative_with_entry': [int((neg.mlen_field == 1).sum()), len(neg)]}
res['margin_length_entry_by_period'] = mb

# first operation of each day: clock start time by hospital and year (stability of the recorded interval)
st = r.or_start_time.astype('string').str.extract(r'^(\d{1,2}):(\d{2})').astype(float)
r['start_min_clock'] = st[0] * 60 + st[1]
fc = r[(r.day_position == 1) & r.start_min_clock.notna()]
clk = {}
for (hc, yr), g in fc.groupby(['centre', 'year']):
    m = g.start_min_clock.median(); q1, q3 = g.start_min_clock.quantile([.25, .75])
    fmt = lambda v: f"{int(v // 60):02d}:{int(round(v % 60)):02d}"
    clk[f'{hc} {yr}'] = {'n_days': len(g), 'median': fmt(m), 'q1': fmt(q1), 'q3': fmt(q3)}
res['first_case_start_clock_by_hospital_year'] = clk
conc = {}
for (hc, yr), g in r.groupby(['centre', 'year']):
    a = g.or_time_min_clock_concordant
    conc[f'{hc} {yr}'] = {'concordant': int((a == 1).sum()), 'assessed': int(a.notna().sum()), 'median_or_time': float(g.or_time.median())}
res['clock_concordance_by_hospital_year'] = conc

# EAU PSA persistence: dependence on the every-value rule
p = pd.read_parquet(REL / 'psa_long.parquet', columns=['master_patient_id', 'days_from_surgery', 'comparator', 'endpoint_use_decision',
                                                     'psa_point_value_ng_ml', 'psa_upper_bound_ng_ml'])
p = p[p.master_patient_id.isin(set(r.pid))]
psa_long_ids = set(p.master_patient_id)
p = p[p.days_from_surgery.between(28, 56)]
psa_window_ids = set(p.master_patient_id)
use = p.endpoint_use_decision.astype(str)
p = p[use.str.contains('eligible_for') & ~use.str.startswith(('held_', 'excluded_'))]
cmpv = p.comparator.astype(str).str.strip()
point = pd.to_numeric(p.psa_point_value_ng_ml, errors='coerce').astype(float).to_numpy(); ub = pd.to_numeric(p.psa_upper_bound_ng_ml, errors='coerce').astype(float).to_numpy()
cmpv = cmpv.to_numpy()
p['ge01'] = np.where(cmpv == '<', np.where(ub <= 0.1, 0.0, np.nan), np.where(~np.isnan(point), (point >= 0.1).astype(float), np.nan))
agg = p.groupby('master_patient_id').agg(n=('ge01', 'size'), n_known=('ge01', lambda x: x.notna().sum()), any_ge=('ge01', lambda x: (x == 1).any()),
                                        all_ge=('ge01', lambda x: x.notna().all() and (x == 1).all()))
agg = agg.join(r.set_index('pid')[['psa_persistence_eau', 'series']])
agg = agg[agg.psa_persistence_eau.notna()]
res['psa_rule'] = {
    'classified_patients': int(len(agg)),
    'recomputed_all_rule_agrees_with_release': int(((agg.all_ge.astype(float)) == agg.psa_persistence_eau).sum()),
    'patients_with_two_or_more_values': {s: int(((agg.series == s) & (agg.n >= 2)).sum()) for s in ['Versius', 'da Vinci restart']},
    'not_persistent_under_all_rule_but_one_value_ge_0_1': {s: int(((agg.series == s) & agg.any_ge & ~agg.all_ge.astype(bool) & (agg.psa_persistence_eau == 0)).sum())
                                                           for s in ['Versius', 'da Vinci restart']},
    'persistent_by_series': {s: [int(agg[(agg.series == s)].psa_persistence_eau.sum()), int((agg.series == s).sum())] for s in ['Versius', 'da Vinci restart']}}

# Apply the PSA outcome table's follow-up rule to the release classifications.
close = (d.surgery_date + pd.to_timedelta(d.fu_months_contact * 30.4375, unit='D')).max()
release_classified = r[((close - r.surgery_date).dt.days >= 56) & r.psa_persistence_eau.notna()]
audited_ids = set(agg.index)
psa_rule = res['psa_rule']
psa_rule['persistent_in_release_by_series'] = {}
psa_rule['classified_with_values_in_audit_by_series'] = {}
psa_rule['classified_without_audited_values'] = {}
psa_rule['classified_in_audit_outside_table_followup_by_series'] = {}
for s, (_, audit_n) in psa_rule['persistent_by_series'].items():
    g = release_classified[release_classified.series == s]
    missing = g.loc[~g.pid.isin(audited_ids), 'pid']
    has_long_row = missing.isin(psa_long_ids)
    in_window = missing.isin(psa_window_ids)
    psa_rule['persistent_in_release_by_series'][s] = [int(g.psa_persistence_eau.sum()), len(g)]
    psa_rule['classified_with_values_in_audit_by_series'][s] = audit_n
    psa_rule['classified_without_audited_values'][s] = {
        'total': len(missing),
        'no_psa_long_row': int((~has_long_row).sum()),
        'psa_long_rows_but_none_on_days_28_to_56': int((has_long_row & ~in_window).sum()),
        'rows_on_days_28_to_56_but_none_eligible': int(in_window.sum()),
    }
    psa_rule['classified_in_audit_outside_table_followup_by_series'][s] = int(
        ((agg.series == s) & ~agg.index.isin(g.pid)).sum())

# missingness by variable and series (observed values)
mv = [('Age', 'age'), ('Body mass index', 'bmi'), ('Preoperative PSA', 'psa'), ('Biopsy ISUP grade group', 'biopsy_isup'),
      ('Specimen weight', 'specimen_weight_g'), ('pT stage', 'pT_group'), ('Final ISUP grade group', 'isup_final'), ('Nerve sparing', 'nerve_sparing'),
      ('Lymph-node dissection', 'plnd'), ('Operative time', 'or_time'), ('Margin status', 'psm'), ('Estimated blood loss', 'ebl'),
      ('Postoperative stay', 'los'), ('Clavien-Dindo grade', 'cd_any'), ('Reoperation', 'reoperation'), ('Transfusion within 30 days', 'transfusion_30d'),
      ('Readmission within 30 days', 'readmission_30d')]
mrows = []
for lab, v in mv:
    row = {'Variable': lab}
    for sr in ['Versius', 'da Vinci restart']:
        g = r[r.series == sr]; row[f'{sr}, observed'] = frac(g[v].notna().sum(), len(g))
    mrows.append(row)
M = pd.DataFrame(mrows)
mmd = '| ' + ' | '.join(M.columns) + ' |\n|' + '|'.join('---' for _ in M.columns) + '|\n' + '\n'.join('| ' + ' | '.join(str(x) for x in row) + ' |' for row in M.values)
(WS / 'tables/tableS_missingness.md').write_text(mmd + '\n')
res['missingness'] = M.to_dict(orient='records')
res['data_close'] = str(pd.to_datetime(d.fu_last_contact_date).max().date())
ym = lambda m: r.loc[m, 'surgery_date'].min().strftime('%Y-%m')
res['robot_start_months'] = {'Versius SM': ym((r.platform == 'Versius') & (r.centre == 'SM')), 'Versius BE': ym((r.platform == 'Versius') & (r.centre == 'BE')),
                             'restart': ym(r.series == 'da Vinci restart'), 'da Vinci BE': ym((r.platform == 'da Vinci') & (r.centre == 'BE'))}
res['robot_model_all_da_vinci'] = r.loc[r.platform == 'da Vinci', 'robot_model'].astype(str).value_counts().to_dict()

for f in [WS / 'tables/tableS_verification.md', WS / 'tables/tableS_missingness.md']:
    if re.search(r'\b(19|20)\d\d-\d\d-\d\d\b|mpid_', f.read_text()): sys.exit(f'privacy scan failed {f}')
(WS / 'analysis/results/R_endpoint_verification.json').write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str))
print(md); print(); print(json.dumps({k: res[k] for k in ['first_case_start_clock_by_hospital_year', 'clock_concordance_by_hospital_year', 'psa_rule', 'data_close']}, indent=1, ensure_ascii=False))
print(); print(mmd)
