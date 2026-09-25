#!/usr/bin/env python3
"""Table 3 (post hoc, labelled): positive margins by hospital, robot and calendar period, with the share of
structured margin reports. Aggregate only."""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
WS = Path(__file__).resolve().parents[1]
d = pd.read_parquet(WS / 'analysis/restricted/ms6_analysis_dataset.parquet'); r = d[d.robotic == 1].copy()
h1 = r.surgery_date.dt.month <= 6
r['period'] = np.select([r.year <= 2023, r.year == 2024, (r.year == 2025) & h1], ['2021–2023', '2024', 'Jan–Jun 2025'], 'Jul 2025–Feb 2026')
r['robot'] = np.where(r.platform == 'Versius', 'Versius', np.where(r.platform_n < 20, 'da Vinci, cases 1–19', 'da Vinci restart'))
r['structured'] = (r.psm_length_relation.fillna('').replace('None', '') != '').astype(float)
def np_(x):
    x = x.dropna(); return f"{int(x.sum())}/{len(x)} ({100 * x.mean():.0f}%)" if len(x) else '–'
rows = []
order = ['2021–2023', '2024', 'Jan–Jun 2025', 'Jul 2025–Feb 2026']
for c, cn in [('SM', 'SalveMedica'), ('BE', 'Bełchatów')]:
    for per in order:
        for rb in ['Versius', 'da Vinci restart', 'da Vinci, cases 1–19']:
            g = r[(r.centre == c) & (r.period == per) & (r.robot == rb)]
            if len(g) == 0: continue
            rows.append({'Hospital': cn, 'Period': per, 'Robot': rb, 'Operations': len(g),
                         'Structured margin report': f"{100 * g.structured.mean():.0f}%",
                         'Positive margin, pT2': np_(g.psm_pt2), 'Positive margin, all stages': np_(g.psm),
                         'pT3 or higher': np_(g.pT3plus), 'Nerve sparing': np_(g.nerve_sparing)})
T = pd.DataFrame(rows)
T.to_csv(WS / 'tables/table3_margins_hospital_period.csv', index=False)
md = '| ' + ' | '.join(T.columns) + ' |\n|' + '|'.join('---' for _ in T.columns) + '|\n' + '\n'.join('| ' + ' | '.join(str(v) for v in row) + ' |' for row in T.values)
(WS / 'tables/table3_margins_hospital_period.md').write_text(md)
if re.search(r'\b(19|20)\d\d-\d\d-\d\d\b', md): sys.exit('privacy scan failed')
from math import sqrt
def wilson(x, n, z=1.959964):
    p = x / n; den = 1 + z * z / n; c = (p + z * z / (2 * n)) / den; h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den; return c - h, c + h
def newcombe(x1, n1, x2, n2):
    p1, p2 = x1 / n1, x2 / n2; l1, u1 = wilson(x1, n1); l2, u2 = wilson(x2, n2); dd = p1 - p2
    return dd, dd - sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), dd + sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)
def cnt(mask):
    x = r.loc[mask, 'psm_pt2'].dropna(); return int(x.sum()), len(x)
SM, BE = r.centre == 'SM', r.centre == 'BE'
vers, rest = r.robot == 'Versius', r.robot == 'da Vinci restart'
early_rest = (r.year == 2024) | ((r.year == 2025) & h1)
contr = [('SalveMedica, restart Jan 2024–Jun 2025 vs Versius 2021–2023', cnt(SM & rest & early_rest), cnt(SM & vers & (r.year <= 2023))),
         ('SalveMedica, restart Jul 2025–Feb 2026 vs Versius 2021–2023', cnt(SM & rest & ~early_rest), cnt(SM & vers & (r.year <= 2023))),
         ('Bełchatów, Versius 2024 vs Versius 2021–2023', cnt(BE & vers & (r.year == 2024)), cnt(BE & vers & (r.year <= 2023))),
         ('Bełchatów, restart (all) vs Versius 2024', cnt(BE & rest), cnt(BE & vers & (r.year == 2024)))]
lines = ['| Post hoc contrast, pT2 positive margin | Events | Difference, points (95% CI) |', '|---|---|---|']
out = {}
for lab, a, b in contr:
    dd, lo, hi = newcombe(a[0], a[1], b[0], b[1])
    lines.append(f"| {lab} | {a[0]}/{a[1]} vs {b[0]}/{b[1]} | {100 * dd:.1f} ({100 * lo:.1f} to {100 * hi:.1f}) |")
    out[lab] = {'events_a': a, 'events_b': b, 'diff': dd, 'ci95': [lo, hi]}
(WS / 'tables/table3_posthoc_contrasts.md').write_text('\n'.join(lines))
import json; (WS / 'analysis/results/table3_posthoc_contrasts.json').write_text(json.dumps(out, indent=1))
print(md); print('\n'.join(lines))
