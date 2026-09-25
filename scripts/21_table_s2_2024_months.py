#!/usr/bin/env python3
"""Table S2: operations in 2024 by month, hospital and robot, from the aggregate 2024 overlap output."""
import json
from pathlib import Path
WS = Path(__file__).resolve().parents[1]
o = json.load(open(WS / 'analysis/results/S_overlap_2024.json'))
cnt = {}
for r in o['monthly_counts']:
    cnt[(r['hospital'], int(r['month']), r['robot'])] = int(r['Freq'])
months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
cols = [('SM', 'Versius'), ('SM', 'restart'), ('BE', 'Versius'), ('BE', 'restart')]
hdr = '| Month, 2024 | SalveMedica, Versius | SalveMedica, da Vinci restart | Bełchatów, Versius | Bełchatów, da Vinci restart |'
rows = [hdr, '|---|---|---|---|---|']
tot = [0] * 4
for i, m in enumerate(months, 1):
    vals = [cnt.get((h, i, r), 0) for h, r in cols]
    tot = [a + b for a, b in zip(tot, vals)]
    rows.append(f'| {m} | ' + ' | '.join(str(v) if v else '–' for v in vals) + ' |')
rows.append('| Total | ' + ' | '.join(str(v) for v in tot) + ' |')
(WS / 'tables/tableS2_2024_by_month.md').write_text('\n'.join(rows) + '\n')
print('\n'.join(rows))
