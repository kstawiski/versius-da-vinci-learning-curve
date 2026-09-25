#!/usr/bin/env python3
"""Print aggregate primary results compactly (no patient-level content exists in these files)."""
import json, sys
from pathlib import Path
R = Path(__file__).resolve().parents[1] / 'analysis/results'
pl = json.load(open(R / 'A2_plateaus.json'))
print('== A1/A2 plateaus ==')
for k, v in pl.items():
    b = v.get('boot_primary', {})
    b10 = v.get('boot_sens10', {})
    f = (lambda x: f"{x:.3f}") if 'psm' in k else (lambda x: f"{x:.1f}")
    print(f"{k:22s} n={v['n_model']} n*={v['nstar']} (10min {v.get('nstar_sens10')}) boot med {b.get('median')} CI {b.get('ci95')} never {b.get('n_never')}/{b.get('n_ok')} case1share {b.get('share_case1')} | "
          f"c1 {f(v['value_case1'])} c25 {f(v['value_case25'])} c50 {f(v['value_case50'])} c100 {f(v['value_case100'])} last {f(v['value_last'])} asym {f(v['asymptote'])} "
          f"slope100/10 {v['slope_last100_per10']} edf {v['edf_median']:.2f} p {v['smooth_p_median']:.2g}/{v['smooth_p_max']:.2g}"
          + (f" | boot10 med {b10.get('median')} CI {b10.get('ci95')}" if b10 else ''))
c = json.load(open(R / 'A3_A4_B1_contrasts.json'))
def show(x, scale=1):
    if x.get('not_estimable'): return f"{x['label']}: NOT ESTIMABLE (n {x['n_ref']}/{x['n_cmp']})"
    s = scale
    extra = f" events {x.get('events_ref')}/{x['n_ref']} vs {x.get('events_cmp')}/{x['n_cmp']}" if 'events_ref' in x else f" n {x['n_ref']} vs {x['n_cmp']} crude {x['crude_ref']:.2f} vs {x['crude_cmp']:.2f} med {x.get('median_ref')} vs {x.get('median_cmp')}"
    eq = f" EQUIV={x['equivalent']} (±{x['eq_margin']})" if 'equivalent' in x else ''
    return f"{x['label']}: {x['estimate']*s:.2f} (95% {x['ci95'][0]*s:.2f} to {x['ci95'][1]*s:.2f}; 90% {x['ci90'][0]*s:.2f} to {x['ci90'][1]*s:.2f}) p={x['p']:.3g}{eq}{extra}"
print('\n== A3 ==')
for k in ['a_or', 'a_psm_pt2', 'a_psm_all']:
    print(show(c['A3'][k], 1 if k == 'a_or' else 100))
for t in c['A3']['slopes_first50']:
    print(f"  {t['term']}: {t['est']:.2f} ({t['ci95'][0]:.2f} to {t['ci95'][1]:.2f}) p={t['p']:.3g}")
print('  A3c', c['A3']['c'])
for blk in ['A4', 'B1']:
    print(f'\n== {blk} ==')
    for k, x in c[blk].items():
        if k in ('rare', 'observation', 'calendar_adjusted', 'definitions'): continue
        print(show(x, 1 if k in ('or_time', 'los', 'ebl') else 100))
    if 'calendar_adjusted' in c[blk]:
        for k, x in c[blk]['calendar_adjusted'].items(): print('  cal:', show(x, 1 if k == 'or_time' else 100))
    for k, x in c[blk]['rare'].items(): print(f"  rare {k}: {x['events']}/{x['n']} (missing {x['n_missing']})")
    print('  observation', c[blk]['observation'])
    if 'definitions' in c[blk]: print('  defs', c[blk]['definitions'])
print('\ndata close gap days', c['data_close_days_after_last_surgery'])
