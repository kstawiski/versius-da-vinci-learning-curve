#!/usr/bin/env python3
"""Stage 1 exploratory data analysis (aggregate only).

Case mix by platform, hospital and platform-case block; outcome completeness by
platform, hospital and year; outliers. Outcome VALUES by case order are not
computed here (prespecification of the learning-curve analysis). Writes CSVs to
tables/eda/ and a JSON summary to evidence/EDA_SUMMARY.json.
"""
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd

WS = Path(__file__).resolve().parents[1]
d = pd.read_parquet(WS / 'analysis/restricted/ms6_analysis_dataset.parquet')
OUT = WS / 'tables/eda'; OUT.mkdir(parents=True, exist_ok=True)
rob = d[d.robotic == 1].copy()
rob['block50'] = ((rob['platform_n'] - 1) // 50 * 50 + 1).astype(int).astype(str) + '-' + (((rob['platform_n'] - 1) // 50 + 1) * 50).astype(int).astype(str)

def med_iqr(x):
    x = x.dropna()
    if len(x) == 0: return ''
    return f"{x.median():.1f} ({x.quantile(.25):.1f}-{x.quantile(.75):.1f}) n={len(x)}"
def pct(x, val=1):
    x = x.dropna()
    return f"{int((x == val).sum())}/{len(x)} ({100 * (x == val).mean():.1f}%)" if len(x) else ''

cont = ['age', 'bmi', 'psa', 'specimen_weight_g']
binv = ['pT3plus', 'nerve_sparing', 'plnd', 'node_positive', 'prior_turp']
def casemix(g):
    r = {'n': len(g)}
    for c in cont: r[c] = med_iqr(g[c])
    r['biopsy_isup>=3'] = pct((g['biopsy_isup'] >= 3).where(g['biopsy_isup'].notna()).astype(float))
    r['isup_final>=3'] = pct((g['isup_final'] >= 3).where(g['isup_final'].notna()).astype(float))
    for c in binv: r[c] = pct(g[c])
    r['ct_T3+'] = pct(g['ct'].str.contains('T3|T4', na=False).where(g['ct'].notna()).astype(float))
    r['ct_T1'] = pct(g['ct'].str.contains('T1', na=False).where(g['ct'].notna()).astype(float))
    for k in ['low', 'favourable_intermediate', 'unfavourable_intermediate', 'high', 'locally_advanced']:
        r[f'eau_{k}'] = pct((g['eau_risk'] == k).where(g['eau_risk'].notna()).astype(float))
    return pd.Series(r)

t1 = rob.groupby('platform').apply(casemix).T
t1['all_robotic'] = casemix(rob)
t1.to_csv(OUT / 'table1_casemix_by_platform.csv')
t1c = rob.groupby(['platform', 'centre']).apply(casemix).T
t1c.columns = [f'{a}|{b}' for a, b in t1c.columns]
t1c.to_csv(OUT / 'casemix_by_platform_centre.csv')
tb = rob.groupby(['platform', 'block50'], sort=False).apply(casemix).T
tb.columns = [f'{a}|{b}' for a, b in tb.columns]
tb.to_csv(OUT / 'casemix_by_platform_block50.csv')
t24 = rob[rob.overlap_2024 == 1].groupby('platform').apply(casemix).T
t24.to_csv(OUT / 'casemix_overlap2024.csv')

# Outcome completeness (share populated) by platform x centre x year
outc = ['or_time', 'ebl', 'los', 'cd30', 'cd90', 'cd_any', 'reoperation', 'urine_leak', 'readmission_30d',
        'psm', 'ln_yield', 'pad_free_3m', 'pad_free_12m', 'psa_persistence_eau', 'psa_persistence_nccn', 'specimen_weight_g', 'bmi', 'ct']
comp = rob.groupby(['platform', 'centre', 'year'])[outc].agg(lambda x: x.notna().mean()).round(3)
comp.insert(0, 'n', rob.groupby(['platform', 'centre', 'year']).size())
comp.to_csv(OUT / 'completeness_by_platform_centre_year.csv')
# Clavien-Dindo documentation pattern by platform x centre x year (grade distribution incl. blank)
cdd = rob.assign(cd30c=rob['cd30'].map({0: 'none', 1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V'}).fillna('blank'))
cdtab = pd.crosstab([cdd.platform, cdd.centre, cdd.year], cdd.cd30c)
cdtab.to_csv(OUT / 'clavien30_documentation_by_platform_centre_year.csv')
# margin-length reporting format by platform x centre
ml = pd.crosstab([rob.platform, rob.centre], rob.psm_length_relation)
ml.to_csv(OUT / 'margin_length_format_by_platform_centre.csv')
# Outliers (counts only)
outl = {}
for c, lo, hi in [('or_time', 60, 360), ('los', 1, 10), ('ebl', 0, 1500), ('psa', 0, 50), ('specimen_weight_g', 15, 150), ('bmi', 17, 42), ('ln_yield', 0, 40)]:
    x = rob[c]
    outl[c] = {'below': int((x < lo).sum()), 'above': int((x > hi).sum()), 'min': float(x.min()), 'max': float(x.max()), 'bounds': [lo, hi]}
# follow-up by platform
fu = rob.groupby('platform')['fu_months_contact'].describe().round(1)
fu.to_csv(OUT / 'followup_by_platform.csv')
summ = {'outliers': outl,
        'or_time_skew_by_platform': rob.groupby('platform')['or_time'].skew().round(2).to_dict(),
        'same_day_cases_dist': rob['same_day_cases'].value_counts().sort_index().to_dict(),
        'files': sorted(p.name for p in OUT.glob('*.csv'))}
(WS / 'evidence/EDA_SUMMARY.json').write_text(json.dumps(summ, indent=1, default=str))
for p in list(OUT.glob('*.csv')) + [WS / 'evidence/EDA_SUMMARY.json']:
    if re.search(r'\b(19|20)\d\d-\d\d-\d\d\b', p.read_text()):
        p.unlink(); sys.exit(f'privacy scan failed: {p.name}')
print(json.dumps(summ, indent=1, default=str))
