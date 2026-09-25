#!/usr/bin/env python3
"""Manuscript tables (aggregate only): Table 1 case mix, Table 2 outcomes by series and phase,
Supplementary tables for margins by period and hospital reporting. Writes tables/*.md and *.csv."""
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
d = pd.read_parquet(WS / 'analysis/restricted/ms6_analysis_dataset.parquet')
close = (d.surgery_date + pd.to_timedelta(d.fu_months_contact * 30.4375, unit='D')).max()
d['days_open'] = (close - d.surgery_date).dt.days
d['series'] = np.select([d.platform == 'Versius', (d.platform == 'da Vinci') & (d.platform_n < 20), d.platform == 'da Vinci'],
                        ['Versius', 'da Vinci, cases 1-19', 'da Vinci restart'], 'Laparoscopic')
OUT = WS / 'tables'; OUT.mkdir(exist_ok=True)

def mi(x, dec=0):
    x = pd.to_numeric(x, errors='coerce').dropna()
    if not len(x): return '-'
    f = f"{{:.{dec}f}}"
    return f"{f.format(x.median())} ({f.format(x.quantile(.25))}-{f.format(x.quantile(.75))})"
def npct(mask, denom_mask=None):
    m = mask if denom_mask is None else mask[denom_mask]
    m = m.dropna().astype(float)
    return f"{int(m.sum())}/{len(m)} ({100 * m.mean():.0f}%)" if len(m) else '-'

def t1(g):
    r = {}
    r['Operations, n'] = str(len(g))
    r['Age, years'] = mi(g.age)
    r['Body mass index, kg/m²'] = mi(g.bmi, 1)
    r['Preoperative PSA, ng/mL'] = mi(g.psa, 1)
    r['Specimen weight, g'] = mi(g.specimen_weight_g)
    r['Biopsy ISUP grade group ≥3'] = npct(g.biopsy_isup.where(g.biopsy_isup.isna(), g.biopsy_isup >= 3))
    er = g.eau_risk
    for lab, keys in [('EAU risk low', ['low']), ('EAU risk intermediate', ['favourable_intermediate', 'unfavourable_intermediate']),
                      ('EAU risk high or locally advanced', ['high', 'locally_advanced'])]:
        r[lab] = npct(er.where(er.isna(), er.isin(keys)))
    r['pT3 or higher'] = npct(g.pT3plus)
    r['Final ISUP grade group ≥3'] = npct(g.isup_final.where(g.isup_final.isna(), g.isup_final >= 3))
    r['pN1 (among dissections)'] = npct(g.node_positive, g.plnd == 1)
    r['Nerve sparing'] = npct(g.nerve_sparing)
    r['Pelvic lymph-node dissection'] = npct(g.plnd)
    r['Operated at Bełchatów'] = npct((g.centre == 'BE').astype(float))
    r['Follow-up to last contact, months'] = mi(g.fu_months_contact, 1)
    return pd.Series(r)
order = ['Versius', 'da Vinci restart', 'da Vinci, cases 1-19']
T1 = pd.DataFrame({s: t1(d[d.series == s]) for s in order})
T1['2024, Versius'] = t1(d[(d.year == 2024) & (d.platform == 'Versius')])
T1['2024, da Vinci restart'] = t1(d[(d.year == 2024) & (d.series == 'da Vinci restart')])
T1.to_csv(OUT / 'table1_case_mix.csv')

rob = d[d.robotic == 1].copy()
rob['phase'] = np.where(rob.platform == 'Versius',
                        'Versius ' + pd.cut(rob.platform_n, [0, 50, 100, 200, 337], labels=['1-50', '51-100', '101-200', '201-337']).astype(str),
                        np.where(rob.platform_n < 20, 'da Vinci 1-19',
                                 'Restart ' + pd.cut(rob.platform_n - 19, [0, 50, 100, 200, 382], labels=['1-50', '51-100', '101-200', '201-382']).astype(str)))
def t2(g):
    r = {}
    r['Operations, n'] = str(len(g))
    r['Operative time, min'] = mi(g.or_time)
    r['Estimated blood loss, mL'] = mi(g.ebl)
    r['Postoperative stay, days'] = mi(g.los)
    r['Positive margin, all stages'] = npct(g.psm)
    r['Positive margin, pT2'] = npct(g.psm_pt2)
    r['Positive margin, pT3 or higher'] = npct(g.psm_pt3)
    r['Nodes removed (dissections)'] = mi(g.ln_yield)
    r['Clavien-Dindo grade ≥III'] = npct(g.major_cd_any)
    r['Reoperation'] = npct(g.reoperation)
    r['Transfusion within 30 days'] = npct(g.transfusion_30d)
    r['Readmission within 30 days'] = npct(g.readmission_30d)
    r['PSA persistence (EAU)'] = npct(g.psa_persistence_eau, g.days_open >= 56)
    r['At most one pad per day, 3 months'] = npct(g.pads_3m.where(g.pads_3m.isna(), g.pads_3m <= 1), g.days_open >= 120)
    r['At most one pad per day, 12 months'] = npct(g.pads_12m.where(g.pads_12m.isna(), g.pads_12m <= 1), g.days_open >= 425)
    r['Pad-free, 12 months'] = npct(g.pad_free_12m, g.days_open >= 425)
    return pd.Series(r)
cols = {'Versius': rob[rob.platform == 'Versius'], 'da Vinci restart': rob[rob.series == 'da Vinci restart']}
for ph in ['Versius 1-50', 'Versius 51-100', 'Versius 101-200', 'Versius 201-337', 'Restart 1-50', 'Restart 51-100', 'Restart 101-200', 'Restart 201-382']:
    cols[ph] = rob[rob.phase == ph]
T2 = pd.DataFrame({k: t2(v) for k, v in cols.items()})
T2.to_csv(OUT / 'table2_outcomes.csv')

# Supplementary: margins by calendar period and hospital, and structured reporting
rob['period'] = np.select([rob.year <= 2023, rob.year == 2024], ['2021-2023', '2024'], '2025-2026')
rob['structured'] = (rob.psm_length_relation.fillna('').replace('None', '') != '').astype(float)
rows = []
for (per, c, pf), g in rob.groupby(['period', 'centre', 'platform']):
    rows.append({'Period': per, 'Hospital': 'SalveMedica' if c == 'SM' else 'Bełchatów', 'Robot': pf, 'Operations': len(g),
                 'Structured margin report': f"{100 * g.structured.mean():.0f}%", 'Positive margin, pT2': npct(g.psm_pt2),
                 'Positive margin, all stages': npct(g.psm), 'Nerve sparing': npct(g.nerve_sparing), 'pT3 or higher': npct(g.pT3plus)})
S1 = pd.DataFrame(rows)
S1.to_csv(OUT / 'tableS_margins_by_period.csv', index=False)

def to_md(df, index=True):
    df = df.reset_index() if index else df
    head = '| ' + ' | '.join(str(c) for c in df.columns) + ' |'
    sep = '|' + '|'.join('---' for _ in df.columns) + '|'
    body = ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in df.values]
    return '\n'.join([head, sep] + body)
(OUT / 'table1_case_mix.md').write_text(to_md(T1.rename_axis('Characteristic')))
(OUT / 'table2_outcomes.md').write_text(to_md(T2.rename_axis('Outcome')))
(OUT / 'tableS_margins_by_period.md').write_text(to_md(S1, index=False))
for p in OUT.glob('*'):
    if p.is_file() and re.search(r'\b(19|20)\d\d-\d\d-\d\d\b', p.read_text(errors='ignore')):
        p.unlink(); sys.exit(f'privacy scan failed {p}')
print(T1.to_string()); print(); print(T2.to_string())
