#!/usr/bin/env python3
"""Single authoritative table of every between-group contrast (A3a, A4, A4 calendar-adjusted, B1).
Binary outcomes use the full adjusted model unless the prespecified sparse rule applies (either group with
fewer than 10 events, or all events), in which case the Firth estimate is used; zero-event outcomes are crude
only. Every figure, table and summary reads this file. Aggregate only."""
import csv, json
from pathlib import Path
R = Path(__file__).resolve().parents[1] / 'analysis/results'
c = json.load(open(R / 'A3_A4_B1_contrasts.json'))
sp = json.load(open(R / 'A4_B1_sparse_binary.json'))
CONT = {'or_time', 'los', 'ebl'}
rows = []
def add(block, outcome, x, scale):
    if x is None: return
    if x.get('not_estimable'):
        rows.append(dict(block=block, outcome=outcome, source='not estimable', n_ref=x.get('n_ref'), n_cmp=x.get('n_cmp'))); return
    rows.append(dict(block=block, outcome=outcome, source='full model', est=x['estimate'] * scale, lo95=x['ci95'][0] * scale, hi95=x['ci95'][1] * scale,
                     lo90=x['ci90'][0] * scale, hi90=x['ci90'][1] * scale, p=x['p'], n_ref=x['n_ref'], n_cmp=x['n_cmp'],
                     events_ref=x.get('events_ref'), events_cmp=x.get('events_cmp'), crude_ref=x.get('crude_ref'), crude_cmp=x.get('crude_cmp'),
                     equivalent=x.get('equivalent'), eq_margin=(x.get('eq_margin') or 0) * (scale if outcome not in CONT else 1) if x.get('eq_margin') else None))
for k, x in c['A3'].items():
    if k.startswith('a_'): add('A3a initial phase', {'a_or': 'or_time', 'a_psm_pt2': 'psm_pt2', 'a_psm_all': 'psm_all'}[k], x, 1 if k == 'a_or' else 100)
for blk, lab in [('A4', 'A4 concurrent 2024'), ('B1', 'B1 mature phases')]:
    for k, x in c[blk].items():
        if k in ('rare', 'observation', 'calendar_adjusted', 'definitions', 'ebl'): continue  # blood loss is descriptive only
        skey = f"{blk} {k if k != 'psm_all' else 'psm'}"
        e = sp.get(skey)
        if e is not None and k not in CONT and k not in ('psm_pt2',):
            if 'rd' in e:
                rows.append(dict(block=lab, outcome=k, source='Firth (sparse rule)', est=e['rd'] * 100, lo95=e['ci95'][0] * 100, hi95=e['ci95'][1] * 100,
                                 lo90=e['ci90'][0] * 100, hi90=e['ci90'][1] * 100, p=e['p'], n_ref=e['n_ref'], n_cmp=e['n_cmp'],
                                 events_ref=e['events_ref'], events_cmp=e['events_cmp'])); continue
            if 'no events' in e.get('note', ''):
                rows.append(dict(block=lab, outcome=k, source='crude only (no events)', n_ref=e['n_ref'], n_cmp=e['n_cmp'], events_ref=0, events_cmp=0)); continue
            if 'too small' in e.get('note', ''):
                rows.append(dict(block=lab, outcome=k, source='not estimable', n_ref=e.get('n_ref'), n_cmp=e.get('n_cmp'))); continue
        add(lab, k, x, 1 if k in CONT else 100)
    if blk == 'A4':
        for k, x in c['A4']['calendar_adjusted'].items():
            add('A4 concurrent 2024, calendar-adjusted', k, x, 1 if k in CONT else 100)
b1 = json.load(open(R / 'B1_mature_bootstrap.json'))['posthoc_fixed_cutoff']
for k, x in b1.items():
    if k == 'label': continue
    add('B1f post hoc, after the Versius learning phase', k, x, 1 if k in CONT else 100)
keys = ['block', 'outcome', 'source', 'est', 'lo95', 'hi95', 'lo90', 'hi90', 'p', 'n_ref', 'n_cmp', 'events_ref', 'events_cmp', 'crude_ref', 'crude_cmp', 'equivalent', 'eq_margin']
with open(R / 'authoritative_contrasts.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); [w.writerow({k: r.get(k) for k in keys}) for r in rows]
for r in rows:
    if 'est' in r:
        print(f"{r['block']:40s} {r['outcome']:22s} {r['source']:20s} {r['est']:8.2f} ({r['lo95']:.2f} to {r['hi95']:.2f}) 90% ({r['lo90']:.2f} to {r['hi90']:.2f}) n {r['n_ref']}/{r['n_cmp']}" + (f" ev {r['events_ref']}/{r['events_cmp']}" if r.get('events_ref') is not None else '') + (f" EQ={r['equivalent']}" if r.get('equivalent') is not None else ''))
    else:
        print(f"{r['block']:40s} {r['outcome']:22s} {r['source']} n {r.get('n_ref')}/{r.get('n_cmp')}")
