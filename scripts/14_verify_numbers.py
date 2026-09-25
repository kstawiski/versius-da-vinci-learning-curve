#!/usr/bin/env python3
"""Deterministic number binding: every key manuscript number is recomputed from the aggregate result files
and the exact expected string must appear in the manuscript. Usage: 14_verify_numbers.py MANUSCRIPT.md"""
import json, sys, re
from pathlib import Path
import pandas as pd
WS = Path(__file__).resolve().parents[1]; R = WS / 'analysis/results'
M = Path(sys.argv[1]).read_text()
cv = pd.read_csv(R / 'A1_standardized_curves.csv')
def cur(series, oc, n, col='est'): return float(cv[(cv.series == series) & (cv.outcome == oc) & (cv.n == n)][col].iloc[0])
ac = pd.read_csv(R / 'authoritative_contrasts.csv')
def con(block, oc): return ac[(ac.block == block) & (ac.outcome == oc)].iloc[0]
s = json.load(open(R / 'secondary_sensitivity.json')); ph = json.load(open(R / 'posthoc_learning_metrics.json'))
pl = json.load(open(R / 'A2_plateaus.json')); c3 = json.load(open(R / 'A3_A4_B1_contrasts.json'))
tc = pd.read_csv(R / 'B2_team_curves.csv'); t3 = json.load(open(R / 'table3_posthoc_contrasts.json'))
b1 = json.load(open(R / 'B1_mature_bootstrap.json'))
m = lambda v: f"{v:.1f}".replace('-', '−')
r0 = lambda v: f"{round(v):.0f}".replace('-', '−')
checks = []
def need(label, text): checks.append((label, text, text in M))
# Versius curve
need('V c1', f"{r0(cur('Versius','or',1))} min at case 1")
need('V c1 CI', f"(95% CI {r0(cur('Versius','or',1,'lo'))}–{r0(cur('Versius','or',1,'hi'))})")
need('V c20', f"{r0(cur('Versius','or',20))} min at case 20")
need('V c50', f"{r0(cur('Versius','or',50))} min at case 50")
need('V c50 CI', f"(95% CI {r0(cur('Versius','or',50,'lo'))}–{r0(cur('Versius','or',50,'hi'))})")
vv = cv[(cv.series == 'Versius') & (cv.outcome == 'or') & (cv.n.between(50, 300))].est
need('V 50-300 range', f"between {r0(vv.min())} and {r0(vv.max())} min up to case 300")
sg = s['S_seg']['Versius']
need('seg abstract', f"breakpoint case {r0(sg['breakpoint'])}, 95% CI {r0(sg['boot_ci95'][0])}–{r0(sg['boot_ci95'][1])}")
need('seg results', f"case {r0(sg['breakpoint'])} (bootstrap 95% CI {r0(sg['boot_ci95'][0])}–{r0(sg['boot_ci95'][1])})")
need('seg slopes', f"by {abs(round(sg['slope_before_pct_per10'])):.0f}% per 10 cases before this point and by {abs(sg['slope_after_pct_per10']):.1f}% per 10 cases after it")
pv = ph['Versius']
need('ph improvement', f"total improvement was {r0(pv['point']['improvement'])} min (95% CI {r0(pv['ci95']['improvement'][0])}–{r0(pv['ci95']['improvement'][1])})")
need('ph frac80', f"case {r0(pv['point']['frac80'])} (95% CI {r0(pv['ci95']['frac80'][0])}–{r0(pv['ci95']['frac80'][1])})")
need('plateau', f"plateau rule gave case {pl['Versius or']['nstar']}")
# restart
need('R c1', f"{r0(cur('dVrestart','or',1))} min at restart case 1 (95% CI {r0(cur('dVrestart','or',1,'lo'))}–{r0(cur('dVrestart','or',1,'hi'))})")
need('R asym', f"long-run level of {r0(pl['dVrestart or']['asymptote'])} min")
pd_ = ph['dVrestart']
need('R improvement', f"{pd_['point']['improvement']:.1f} min (95% CI {r0(pd_['ci95']['improvement'][0])} to {r0(pd_['ci95']['improvement'][1])})")
a = con('A3a initial phase', 'or_time')
need('A3a', f"{-a.est:.1f} min shorter than the first 100 Versius operations (95% CI {-a.hi95:.1f}–{-a.lo95:.1f}")
sl = {t['term']: t for t in c3['A3']['slopes_first50']}['case1_diff_min']
rr = cv[(cv.series == 'dVrestart') & (cv.outcome == 'or')]; pk = rr.loc[rr.est.idxmax()]
need('R transient rise', f"rose transiently to {r0(pk.est)} min around restart case {int(pk.n)}, in early 2025, and was {r0(cur('dVrestart','or',300))} min by case 300")
f = con('B1f post hoc, after the Versius learning phase', 'or_time')
need('B1f or', f"{m(f.est)} min, 90% CI {m(f.lo90)} to {m(f.hi90)}")
fp = con('B1f post hoc, after the Versius learning phase', 'psm_pt2')
# 2024 (v3 wording: direction stated explicitly)
x = con('A4 concurrent 2024', 'or_time'); need('A4 or', f"took {m(x.est)} min longer than Versius operations (95% CI {m(x.lo95)} to {m(x.hi95)})")
x = con('A4 concurrent 2024', 'psm_pt2'); need('A4 pT2', f"{m(x.est)} percentage points higher in pT2 disease (95% CI {m(x.lo95)} to {m(x.hi95)})")
x = con('A4 concurrent 2024', 'psm_all'); need('A4 all', f"{m(x.est)} points higher at all stages (95% CI {m(x.lo95)} to {m(x.hi95)})")
x = con('A4 concurrent 2024', 'los'); need('A4 los', f"{x.est:.1f} days longer on da Vinci (95% CI {m(x.lo95).replace('−0.0','−0.0')}")
x = con('A4 concurrent 2024', 'psa_persistence_eau'); need('A4 psa', f"PSA persistence was {m(-x.est)} points lower (difference {m(x.est)}, 95% CI {m(x.lo95)} to {m(x.hi95)}")
a2 = con('A3a initial phase', 'psm_pt2'); need('abstract A3a pT2', f"pT2) disease were {a2.est:.1f} percentage points higher than on Versius (95% CI {a2.lo95:.1f}–{a2.hi95:.1f})")
need('abstract A3a or', f"{-a.est:.1f} min shorter than the first 100 on Versius (95% CI {-a.hi95:.1f}–{-a.lo95:.1f})")
x = con('A4 concurrent 2024', 'or_time'); need('A4 or 90', f"({m(x.lo90)} to {m(x.hi90)} min)")
# margins
a = con('A3a initial phase', 'psm_pt2'); need('A3a pT2', f"{a.est:.1f} percentage points more pT2 margins than the first 100 Versius operations (95% CI {a.lo95:.1f}–{a.hi95:.1f})")
need('B1f pT2 body', f"remained {fp.est:.1f} points higher on the restart (95% CI {fp.lo95:.1f}–{fp.hi95:.1f}")
need('B1f pT2 abstract', f"pT2 margins remained {fp.est:.1f} points higher (95% CI {fp.lo95:.1f}–{fp.hi95:.1f})")
need('R psm c100/c200', f"rose to {r0(100 * cur('dVrestart','psm_pt2',100))}% by case 100 and fell to {r0(100 * cur('dVrestart','psm_pt2',200))}% by case 200")
need('R psm p', f"(P = {pl['dVrestart psm_pt2']['smooth_p_median']:.2f})")
k = list(t3.keys())
need('SM contrasts', f"+{100 * t3[k[0]]['diff']:.1f} points (95% CI {100 * t3[k[0]]['ci95'][0]:.1f}–{100 * t3[k[0]]['ci95'][1]:.1f}) and {m(100 * t3[k[1]]['diff'])} points (95% CI {m(100 * t3[k[1]]['ci95'][0])} to {m(100 * t3[k[1]]['ci95'][1])})")
need('BE contrast', f"difference of {100 * t3[k[2]]['diff']:.1f} points (95% CI {100 * t3[k[2]]['ci95'][0]:.1f}–{100 * t3[k[2]]['ci95'][1]:.1f})")
need('BE restart', f"({m(100 * t3[k[3]]['diff'])} points, 95% CI {m(100 * t3[k[3]]['ci95'][0])} to {m(100 * t3[k[3]]['ci95'][1])})")
cu = s['S_cusum']['Versius_racusum']; need('racusum', f"{cu['observed']} observed against {r0(cu['expected'])} expected")
# team curves
def team(pf, h, x): return float(tc[(tc.platform == pf) & (tc.hospital == h) & (tc.x == x)].est.iloc[0])
need('SM V', f"fell from {r0(team('Versius','SM',1))} min to {r0(team('Versius','SM',50))} min by team case 50 and {r0(team('Versius','SM',100))} min by team case 100")
need('BE V first', f"first operations took {r0(team('Versius','BE',1))} min")
bev = tc[(tc.platform == 'Versius') & (tc.hospital == 'BE') & (tc.x <= 100)].est
need('BE V range', f"between {r0(bev.min())} and {r0(bev.max())} min up to team case 100")
smd = tc[(tc.platform == 'da Vinci') & (tc.hospital == 'SM')]
need('SM dV', f"began at {r0(smd.est.iloc[0])} min (95% CI {r0(smd.lo.iloc[0])}–{r0(smd.hi.iloc[0])}), against a long-run level of {r0(smd.sort_values('x').est.tail(50).mean())} min for that team, and dipped to {r0(smd.est.min())} min by team case {int(smd.x.iloc[smd.est.values.argmin()])}")
bed = tc[(tc.platform == 'da Vinci') & (tc.hospital == 'BE')]
need('BE dV', f"The Bełchatów da Vinci team's first operations took {r0(bed.est.iloc[0])} min, the same as that team's first Versius operations, and its time fell gradually to {r0(bed.est.iloc[-1])} min by team case {int(bed.x.max())}")
checks.append(('BE dV first equals BE V first', f"{r0(bed.est.iloc[0])} == {r0(team('Versius','BE',1))}", r0(bed.est.iloc[0]) == r0(team('Versius','BE',1))))
ts = json.load(open(R / 'B2_team_curves_summary.json'))
need('BE dV start', f"began after {ts['da Vinci|BE']['surgeon_platform_case_at_team_start'] - 20} restart operations")
gv = pd.merge(tc[(tc.platform == 'Versius') & (tc.hospital == 'SM')], tc[(tc.platform == 'Versius') & (tc.hospital == 'BE')], on='x')
gd = pd.merge(smd, bed, on='x')
need('gaps', f"up to {r0((gv.est_y - gv.est_x).max())} min longer on Versius and up to {r0((gd.est_y - gd.est_x).max())} min longer on da Vinci")
need('abstract gap', f"differed by up to {r0(max((gv.est_y - gv.est_x).max(), (gd.est_y - gd.est_x).max()))} min between hospitals")
need('discussion gap', f"differed by up to {r0(max((gv.est_y - gv.est_x).max(), (gd.est_y - gd.est_x).max()))} min between the two hospitals at equal team experience")
need('decline 1-50', f"time changed by {r0(ph['dVrestart']['point']['decline_1_50'])} min on the restart (95% CI {r0(ph['dVrestart']['ci95']['decline_1_50'][0])} to {r0(ph['dVrestart']['ci95']['decline_1_50'][1])}) and fell by {r0(ph['Versius']['point']['decline_1_50'])} min on Versius (95% CI {r0(ph['Versius']['ci95']['decline_1_50'][0])}–{r0(ph['Versius']['ci95']['decline_1_50'][1])})")
cs = json.load(open(R / 'B1f_cutoff_sensitivity.json'))
checks.append(('cutoff 46 rerun equals B1f', f"{cs['after_46']['or_time']['estimate']:.6f} == {f.est:.6f}", abs(cs['after_46']['or_time']['estimate'] - f.est) < 1e-6))
need('cutoff sens', f"differed by {m(cs['after_30']['or_time']['estimate'])} and {m(cs['after_43']['or_time']['estimate'])} min, and pT2 margins were {100 * cs['after_30']['psm_pt2']['estimate']:.1f} and {100 * cs['after_43']['psm_pt2']['estimate']:.1f} points higher")
# figure 1 volume peak
q = pd.read_csv(WS / 'figures/source_data/fig1A_quarterly_counts.csv'); qq = q.groupby('q').n.sum()
peak = qq.idxmax()
checks.append(('volume peak quarter (info)', f"{peak} with {qq.max()} operations", peak.startswith('2025')))
# bootstrap statement
need('boot 98%', f"Both plateau cutoffs were reached in {round(100 * b1['share_reached']):.0f}% of bootstrap replicates, and among these the operative-time difference stayed within the equivalence margin in {round(100 * b1['or_time']['share_within_15']):.0f}%")
ok = all(c[2] for c in checks)
for c in checks: print(('OK   ' if c[2] else 'FAIL ') + c[0] + ' :: ' + c[1])
print('ALL_BOUND' if ok else 'BINDING_FAILURES')
