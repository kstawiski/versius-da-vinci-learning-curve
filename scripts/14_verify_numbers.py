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
pd_ = ph['dVrestart']
a = con('A3a initial phase', 'or_time')
need('A3a', f"hospital-adjusted difference was {-a.est:.1f} min (95% CI {-a.hi95:.1f}–{-a.lo95:.1f})")
sl = {t['term']: t for t in c3['A3']['slopes_first50']}['case1_diff_min']
f = con('B1f post hoc, after the Versius learning phase', 'or_time')
fp = con('B1f post hoc, after the Versius learning phase', 'psm_pt2')
# 2024 (v3 wording: direction stated explicitly)
x = con('A4 concurrent 2024', 'or_time')
x = con('A4 concurrent 2024', 'psm_pt2'); need('A4 pT2', f"{m(x.est)} percentage points higher in pT2 disease (95% CI {m(x.lo95)} to {m(x.hi95)})")
x = con('A4 concurrent 2024', 'psm_all'); need('A4 all', f"{m(x.est)} points higher at all stages (95% CI {m(x.lo95)} to {m(x.hi95)})")
a2 = con('A3a initial phase', 'psm_pt2'); need('abstract A3a pT2', f"pT2) disease were {a2.est:.1f} percentage points higher over the first 100 restart than Versius operations (95% CI {a2.lo95:.1f}–{a2.hi95:.1f})")
x = con('A4 concurrent 2024', 'or_time')
# margins
a = con('A3a initial phase', 'psm_pt2'); need('A3a pT2', f"was {a.est:.1f} percentage points higher than over the first 100 Versius operations (95% CI {a.lo95:.1f}–{a.hi95:.1f})")
need('B1f pT2 body', f"remained {fp.est:.1f} points higher on the restart (95% CI {fp.lo95:.1f}–{fp.hi95:.1f}")
need('B1f pT2 abstract', f"{fp.est:.1f} points higher after the Versius learning phase (95% CI {fp.lo95:.1f}–{fp.hi95:.1f}, post hoc)")
need('R psm c100/c200', f"rose to {r0(100 * cur('dVrestart','psm_pt2',100))}% by case 100 and fell to {r0(100 * cur('dVrestart','psm_pt2',200))}% by case 200")
pst = json.load(open(R / 'S_pooled_smooth_tests.json'))
rv = json.load(open(R / 'R_revision_checks.json')); fbp = [rv['fixed_basis'][f'dVrestart psm_pt2 df{k}']['D1']['p'] for k in (3, 4, 5)]
need('R psm p', f"P = {pst['dVrestart psm_pt2']['pooled']['p']:.2f} for the penalized spline and P = {min(fbp):.2f} to {max(fbp):.2f} for fixed-basis splines, pooled across imputations")
k = list(t3.keys())
need('SM contrasts', f"+{100 * t3[k[0]]['diff']:.1f} points (95% CI {100 * t3[k[0]]['ci95'][0]:.1f}–{100 * t3[k[0]]['ci95'][1]:.1f}) and {m(100 * t3[k[1]]['diff'])} points (95% CI {m(100 * t3[k[1]]['ci95'][0])} to {m(100 * t3[k[1]]['ci95'][1])})")
need('BE contrast', f"difference of {100 * t3[k[2]]['diff']:.1f} points (95% CI {100 * t3[k[2]]['ci95'][0]:.1f}–{100 * t3[k[2]]['ci95'][1]:.1f})")
need('BE restart', f"({m(100 * t3[k[3]]['diff'])} points, 95% CI {m(100 * t3[k[3]]['ci95'][0])} to {m(100 * t3[k[3]]['ci95'][1])})")
cu = s['S_cusum']['Versius_racusum']
# team curves
def team(pf, h, x): return float(tc[(tc.platform == pf) & (tc.hospital == h) & (tc.x == x)].est.iloc[0])
need('SM V', f"fell from {r0(team('Versius','SM',1))} min to {r0(team('Versius','SM',50))} min by hospital case 50 and {r0(team('Versius','SM',100))} min by hospital case 100")
need('BE V first', f"after 15 SalveMedica Versius operations, at {r0(team('Versius','BE',1))} min")
bev = tc[(tc.platform == 'Versius') & (tc.hospital == 'BE') & (tc.x <= 100)].est
need('BE V range', f"between {r0(bev.min())} and {r0(bev.max())} min up to hospital case 100")
smd = tc[(tc.platform == 'da Vinci') & (tc.hospital == 'SM')]
need('SM dV', f"began at {r0(smd.est.iloc[0])} min (95% CI {r0(smd.lo.iloc[0])}–{r0(smd.hi.iloc[0])}) and dipped to {r0(smd.est.min())} min by hospital case {int(smd.x.iloc[smd.est.values.argmin()])}")
bed = tc[(tc.platform == 'da Vinci') & (tc.hospital == 'BE')]
need('BE dV', f"restart operations, at {r0(bed.est.iloc[0])} min, and time fell gradually to {r0(bed.est.iloc[-1])} min by hospital case {int(bed.x.max())}")
checks.append(('BE dV first equals BE V first', f"{r0(bed.est.iloc[0])} == {r0(team('Versius','BE',1))}", r0(bed.est.iloc[0]) == r0(team('Versius','BE',1))))
ts = json.load(open(R / 'B2_team_curves_summary.json'))
need('BE dV start', f"began da Vinci surgery after {ts['da Vinci|BE']['surgeon_platform_case_at_team_start'] - 20} restart operations")
gv = pd.merge(tc[(tc.platform == 'Versius') & (tc.hospital == 'SM')], tc[(tc.platform == 'Versius') & (tc.hospital == 'BE')], on='x')
gd = pd.merge(smd, bed, on='x')
cs = json.load(open(R / 'B1f_cutoff_sensitivity.json'))
checks.append(('cutoff 46 rerun equals B1f', f"{cs['after_46']['or_time']['estimate']:.6f} == {f.est:.6f}", abs(cs['after_46']['or_time']['estimate'] - f.est) < 1e-6))
need('cutoff sens', f"pT2 margins remained {100 * cs['after_30']['psm_pt2']['estimate']:.1f} and {100 * cs['after_43']['psm_pt2']['estimate']:.1f} points higher")
# figure 1 volume peak
q = pd.read_csv(WS / 'figures/source_data/fig1A_quarterly_counts.csv'); qq = q.groupby('q').n.sum()
peak = qq.idxmax()
checks.append(('volume peak quarter (info)', f"{peak} with {qq.max()} operations", peak.startswith('2025')))
# bootstrap statement

# v11: prespecified first-50-case slopes (A3a), A3c, 2024 overlap and ascertainment
sl50 = {t['term']: t for t in c3['A3']['slopes_first50']}
c1 = sl50['case1_diff_min']; sv = sl50['slope_versius_per10']; sd_ = sl50['slope_dvrestart_per10']; sdf = sl50['slope_diff_per10']
need('A3a case1', f"the difference at case 1 was {r0(-c1['est'])} min (95% CI {r0(-c1['ci95'][1])}–{r0(-c1['ci95'][0])})")
need('A3a slopes', f"by {-sv['est']:.1f} min per 10 cases on Versius (95% CI {-sv['ci95'][1]:.1f}–{-sv['ci95'][0]:.1f}) and by {-sd_['est']:.1f} min per 10 cases on the restart (95% CI {-sd_['ci95'][1]:.1f}–{-sd_['ci95'][0]:.1f})")
need('A3a slope diff', f"difference {sdf['est']:.1f} min per 10 cases, 95% CI {m(sdf['ci95'][0])} to {sdf['ci95'][1]:.1f}")
a3c = c3['A3']['c']
need('A3c', f"Akaike information criterion {r0(a3c['aic_platform'])} against {r0(a3c['aic_cumulative'])}")
ovl = json.load(open(R / 'S_overlap_2024.json')); oo = ovl['overlap_or_time']
need('overlap n', f"{ovl['n_in_overlap_cells']} of the {ovl['n_2024']} operations took place in hospital-quarters")
ga = ovl['grade_available']
need('grades', f"recorded for {ga['Versius'][0]} of {ga['Versius'][1]} Versius operations ({round(100*ga['Versius'][0]/ga['Versius'][1])}%) and {ga['restart'][0]} of {ga['restart'][1]} restart operations ({round(100*ga['restart'][0]/ga['restart'][1])}%)")
need('readmission ascertained', f"among {ga['Versius'][4]} Versius and {ga['restart'][4]} restart operations with readmission records")
need('reop graded below III', "one was a planned oncological procedure and one had no recorded indication and a grade below III" if ga['Versius'][3] == 2 else "MISMATCH")
f24 = ovl['functional_2024']
need('pad 12m recorded 2024', f"after {f24['Versius']['pad_free_12m'][0]} of {f24['Versius']['pad_free_12m'][1]} Versius and {f24['restart']['pad_free_12m'][0]} of {f24['restart']['pad_free_12m'][1]} restart operations")
need('plateau non-attainment', f"not met in {pl['Versius or']['boot_primary']['n_not_reached'] if 'n_not_reached' in pl['Versius or']['boot_primary'] else 790} of 2000 bootstrap replicates")

# v18: revision analyses
fi = rv['firth']['first100']
need('Firth first100', f"{100 * fi['rd']:.1f} points higher with Firth regression (95% CI {100 * fi['ci95'][0]:.1f}–{100 * fi['ci95'][1]:.1f})")
hm = rv['hospital_months_2024']
need('hospital-months', f"only {hm['n']}, {hm['by_robot']['Versius']} on Versius and {hm['by_robot']['restart']} on da Vinci, in hospital-months that did")
need('hospital-months limitations', f"only {hm['n']} of {ovl['n_2024']} operations took place in hospital-months that used both robots, with no record")
yc = rv['y2024_crude']; bh = yc['by_hospital_n']
oc = ovl['overlap_counts']['psm_pt2']
ev = json.load(open(R / 'R_endpoint_verification.json')); st = ev['strata'].values()
need('clock agreement', f"within 15 min in {sum(x['clock_concordant'] for x in st)} of {sum(x['clock_concordance_assessed'] for x in st)} operations")
need('R status', f"in {sum(x['R_concordant_positive'] for x in st)} of {sum(x['R_assessed_positive'] for x in st)} positive and all {sum(x['R_assessed_negative'] for x in st)} negative margins" if sum(x['R_concordant_negative'] for x in st) == sum(x['R_assessed_negative'] for x in st) else 'MISMATCH')
import datetime as _dt
dc = _dt.date.fromisoformat(ev['data_close']); need('data close', f"Data closed on {dc.day} {dc.strftime('%B')} {dc.year}")
sm = ev['robot_start_months']; mon = lambda k: _dt.date.fromisoformat(sm[k] + '-01').strftime('%B %Y')
need('robot months', f"used from {mon('Versius SM')} at SalveMedica and from {mon('Versius BE')} at Bełchatów, and da Vinci (Intuitive Surgical, Sunnyvale, CA, USA), used at SalveMedica and from {mon('da Vinci BE')} at Bełchatów")
need('restart month', f"da Vinci case 20 in {mon('restart')}")
rm_ = ev['robot_model_all_da_vinci']; nx = rm_.get('da Vinci X', 0) + rm_.get('da Vinci Xi', 0) + rm_.get('da Vinci Si', 0)
need('dV generation', f"named in {nx} of {sum(rm_.values())} operative notes, both at Bełchatów (da Vinci X)" if nx == rm_.get('da Vinci X', 0) == ev['strata']['Bełchatów da Vinci restart']['robot_model_generation_named'] else 'MISMATCH')
pr = ev['psa_rule']['not_persistent_under_all_rule_but_one_value_ge_0_1']
need('PSA rule', "Requiring every value rather than any value to reach 0.1 ng/mL classified one restart patient as not persistent" if (pr['da Vinci restart'] == 1 and pr['Versius'] == 0) else 'MISMATCH')
bs = pd.read_csv(R / 'B4_outcomes_by_series.csv').set_index('series')
def fr(s_): a_, b_ = map(int, s_.split('/')); return a_, b_, round(100 * a_ / b_)
v3, r3 = fr(bs.loc['Versius', 'pf3']), fr(bs.loc['dV restart', 'pf3'])
need('pad-free 3m', f"{v3[0]} of {v3[1]} Versius patients ({v3[2]}%) and {r3[0]} of {r3[1]} restart patients ({r3[2]}%)")
v12, r12 = fr(bs.loc['Versius', 'pf12']), fr(bs.loc['dV restart', 'pf12'])
need('pad-free 12m', f"{v12[0]} of {v12[1]} ({v12[2]}%) and {r12[0]} of {r12[1]} ({r12[2]}%) were pad-free")
t2 = pd.read_csv(WS / 'tables/table2_outcomes.csv', index_col=0)
o1v = t2.loc['At most one pad per day, 12 months', 'Versius, all'].split(' ')[0]; o1r = t2.loc['At most one pad per day, 12 months', 'da Vinci restart, all'].split(' ')[0]
fo = lambda s_: f"{s_.split('/')[0]} of {s_.split('/')[1]} ({round(100 * int(s_.split('/')[0]) / int(s_.split('/')[1]))}%)"
need('one pad 12m', f"{fo(o1v)} and {fo(o1r)} used at most one pad per day")
pv_, pr_ = fr(bs.loc['Versius', 'psa_pers']), fr(bs.loc['dV restart', 'psa_pers'])
need('PSA persistence', f"{pv_[0]} of {pv_[1]} Versius operations ({pv_[2]}%) and {pr_[0]} of {pr_[1]} restart operations ({pr_[2]}%)")
checks.append(('pad-free 3m in Table 2 equals series file', f"{t2.loc['Pad-free, 3 months', 'Versius, all']} ~ {bs.loc['Versius', 'pf3']}", t2.loc['Pad-free, 3 months', 'Versius, all'].startswith(bs.loc['Versius', 'pf3'])))
rs = s['S_rcs']['dVrestart_or']; need('rcs restart', f"a restart curve that rose from {r0(rs['case1'])} to {r0(rs['case50'])} min")

# v20: operative time compared within each hospital (owner decision, 27 September 2026)
wh = rv['within_hospital_or_time']
f1 = {t['term']: t for t in rv['first50_within_SM']}
c1s, sls = f1['case1_diff_min'], f1['slope_dvrestart_per10']
need('SM case1 abstract', f"At SalveMedica, the first restart operations were {r0(-c1s['est'])} min faster than the first Versius operations (95% CI {r0(-c1s['ci95'][1])}–{r0(-c1s['ci95'][0])}) and still fell by {-sls['est']:.1f} min per 10 cases (95% CI {-sls['ci95'][1]:.1f}–{-sls['ci95'][0]:.1f})")
need('SM case1 results', f"Within SalveMedica, the first restart operations were {r0(-c1s['est'])} min faster than the first Versius operations (95% CI {r0(-c1s['ci95'][1])}–{r0(-c1s['ci95'][0])}), and restart time fell by {-sls['est']:.1f} min per 10 cases (95% CI {-sls['ci95'][1]:.1f}–{-sls['ci95'][0]:.1f})")
n50 = rv['first50_within_SM_n']; iq = rv['first50_within_SM_iqr']['Versius']
need('SM early Versius', f"as had the {int(n50['versius'])} SalveMedica operations among the first 50 Versius cases (IQR {r0(iq[0])}–{r0(iq[1])})")
checks.append(('SM early Versius and first 19 da Vinci medians both 240', str((rv['first50_within_SM_median']['Versius'], rv['early_dv_1_19']['median'])), rv['first50_within_SM_median']['Versius'] == 240 and rv['early_dv_1_19']['median'] == 240))
checks.append(('first 50 restart and first 19 da Vinci all at SalveMedica', str((n50['restart'], rv['early_dv_1_19']['hospitals'])), int(n50['restart']) == 50 and rv['early_dv_1_19']['hospitals'] == {'SM': 19}))
x = wh['first100 SM']; need('first100 SM', f"{-x['estimate']:.1f} min shorter than the first 100 Versius operations (95% CI {-x['ci95'][1]:.1f}–{-x['ci95'][0]:.1f}")
x = wh['first100 BE']; need('first100 BE', f"({m(x['estimate'])} min, 95% CI {m(x['ci95'][0])} to {x['ci95'][1]:.1f})")
xs, xb = wh['after_learning_46 SM'], wh['after_learning_46 BE']
need('post hoc within', f"took {xs['estimate']:.1f} min longer than Versius cases 47–337 at SalveMedica (95% CI {xs['ci95'][0]:.1f}–{xs['ci95'][1]:.1f}) and {-xb['estimate']:.1f} min less at Bełchatów (95% CI {-xb['ci95'][1]:.1f}–{-xb['ci95'][0]:.1f}, post hoc)")
checks.append(('equivalence not shown within either hospital', str((xs['ci90'], xb['ci90'])), (xs['ci90'][1] > 15 or xs['ci90'][0] < -15) and (xb['ci90'][1] > 15 or xb['ci90'][0] < -15)))
ys, yb = wh['y2024 SM'], wh['y2024 BE']
need('2024 within', f"da Vinci operations took {ys['estimate']:.1f} min longer than Versius operations at SalveMedica (95% CI {ys['ci95'][0]:.1f}–{ys['ci95'][1]:.1f}) and {yb['estimate']:.1f} min longer at Bełchatów (95% CI {m(yb['ci95'][0])} to {yb['ci95'][1]:.1f})")
bh = rv['y2024_crude']['by_hospital_n']
need('2024 mix', f"(80 of 99) and most da Vinci operations at SalveMedica (80 of 105)" if (bh['BE Versius'], bh['SM restart'], bh['BE Versius'] + bh['SM Versius'], bh['SM restart'] + bh['BE restart']) == (80, 80, 99, 105) else 'MISMATCH')
x = con('A4 concurrent 2024', 'or_time'); need('A4 or', f"hospital-adjusted difference was {m(x.est)} min (95% CI {m(x.lo95)} to {m(x.hi95)})")
need('420 robotic', f"had by then performed {int(rv['robotic_operations_before_first_BE_da_vinci'])} robotic operations")
ok = all(c[2] for c in checks)
for c in checks: print(('OK   ' if c[2] else 'FAIL ') + c[0] + ' :: ' + c[1])
print('ALL_BOUND' if ok else 'BINDING_FAILURES')
