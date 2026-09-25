#!/usr/bin/env python3
"""Supplementary Table S1: sensitivity and robustness analyses (aggregate). Reads result JSON files only."""
import json
from pathlib import Path
WS = Path(__file__).resolve().parents[1]; R = WS / 'analysis/results'
pl = json.load(open(R / 'A2_plateaus.json')); s = json.load(open(R / 'secondary_sensitivity.json'))
ph = json.load(open(R / 'posthoc_learning_metrics.json')); b1 = json.load(open(R / 'B1_mature_bootstrap.json'))
c = json.load(open(R / 'A3_A4_B1_contrasts.json'))
f = lambda v, d=0: '–' if v is None else ('not reached' if isinstance(v, float) and v == float('inf') else f"{v:.{d}f}")
def upper(b):
    u = b['unconditional']['p97.5']
    return 'not reached' if u is None or u == 'Inf' or (isinstance(u, (int, float)) and u > 10 ** 6) else f"{u:.0f}"
def nstar_text(key):
    v = pl[key]; b = v['boot_primary']
    return (f"case {v['nstar']} (95% CI {f(b['unconditional']['p2.5'])} to {upper(b)}; rule not met in "
            f"{b['n_not_reached']} of {b['n_ok']} replicates)") if v['nstar'] is not None else f"not reached (rule not met in {b['n_not_reached']} of {b['n_ok']} replicates)"
rows = []
rows.append(('Prespecified plateau rule, operative time within 15 min of the final-50-case mean', nstar_text('Versius or'), nstar_text('dVrestart or') + ', marking the end of a transient rise'))
rows.append(('Plateau rule, 10 min', 'not reached' if pl['Versius or']['nstar_sens10'] is None else f"case {pl['Versius or']['nstar_sens10']}",
             'not reached' if pl['dVrestart or']['nstar_sens10'] is None else f"case {pl['dVrestart or']['nstar_sens10']}"))
sg = s['S_seg']; dav = json.load(open(R / 'S_seg_davies_adjusted.json'))
def segtext(k):
    x = sg[k]
    return (f"case {f(x['breakpoint'], 1)} (bootstrap 95% CI {f(x['boot_ci95'][0])}–{f(x['boot_ci95'][1])}, {x['boot_n_failed']} of "
            f"{x['boot_n_ok'] + x['boot_n_failed']} fits failed). {f(x['slope_before_pct_per10'], 1)}% per 10 cases before (95% CI {f(x['slope_before_ci95'][0], 1)} to {f(x['slope_before_ci95'][1], 1)}), "
            f"{f(x['slope_after_pct_per10'], 2)}% after (95% CI {f(x['slope_after_ci95'][0], 2)} to {f(x['slope_after_ci95'][1], 2)}). Davies test P {'< 0.001' if dav[k]['davies_p_max'] < 0.001 else '= ' + str(round(dav[k]['davies_p_max'], 3))} in every imputation")
rows.append(('Segmented regression, adjusted, pooled over 20 imputations', segtext('Versius'), segtext('dVrestart')))
rows.append(('Segmented regression, unadjusted', f"case {f(sg['Versius_unadjusted']['breakpoint'], 1)}", f"case {f(sg['dVrestart_unadjusted']['breakpoint'], 1)}"))
pv, pd_ = ph['Versius'], ph['dVrestart']
def frac(p):
    v = p['point'].get('frac80'); lo, hi = p['ci95'].get('frac80', [None, None])
    return 'not applicable, no initial improvement' if v is None else f"case {v:.0f} (95% CI {f(lo)}–{f(hi)})"
rows.append(('Post hoc, case by which 80% of the initial improvement was reached', frac(pv) + f", total improvement {f(pv['point']['improvement'])} min (95% CI {f(pv['ci95']['improvement'][0])}–{f(pv['ci95']['improvement'][1])})",
             f"no initial learning phase, total improvement {f(pd_['point']['improvement'], 1)} min (95% CI {f(pd_['ci95']['improvement'][0])} to {f(pd_['ci95']['improvement'][1])})"))
r = s['S_rcs']
rows.append(('Restricted cubic splines instead of penalized spline', f"{f(r['Versius_or']['case1'])} min at case 1, {f(r['Versius_or']['case50'])} at case 50, {f(r['Versius_or']['case100'])} at case 100",
             f"{f(r['dVrestart_or']['case1'])} min at case 1, {f(r['dVrestart_or']['case50'])} at case 50"))
n = s['S_noNSLND']
rows.append(('Without nerve-sparing and lymph-node dissection terms', f"{f(n['Versius_or']['case1'])} min at case 1, {f(n['Versius_or']['case50'])} at case 50",
             f"{f(n['dVrestart_or']['case1'])} min at case 1, {f(n['dVrestart_or']['case50'])} at case 50"))
cc = s['S_cc']
rows.append((f"Complete-case analysis ({cc['n_complete']} operations)",
             f"{f(cc['Versius_or']['case1'])} min at case 1, {f(cc['Versius_or']['case50'])} at case 50",
             f"{f(cc['dVrestart_or']['case1'])} min at case 1, {f(cc['dVrestart_or']['case50'])} at case 50. First 100 restart minus first 100 Versius {f(cc['A3a_or']['estimate'], 1)} min (95% CI {f(cc['A3a_or']['ci95'][0], 1)} to {f(cc['A3a_or']['ci95'][1], 1)}). 2024 difference {f(cc['A4_or']['estimate'], 1)} min (95% CI {f(cc['A4_or']['ci95'][0], 1)} to {f(cc['A4_or']['ci95'][1], 1)})"))
dv = s['S_dvall']
rows.append(('da Vinci curve over all 401 da Vinci cases', '–', f"{f(dv['or']['case1'])} min at case 1 (2021), {f(dv['or']['case50'])} at case 50"))
t = s['S_ties']
rows.append((f"Random order within {t['n_blocks']} unresolved tie blocks ({t['n_operations_in_blocks']} operations, {t['n_robotic_in_blocks']} robotic), {t['n_perm']} permutations",
             f"plateau rule case {f(t['Versius_or_nstar_range'][0])}–{f(t['Versius_or_nstar_range'][1])}; first-100 contrast unchanged",
             f"plateau rule case {f(t['dVrestart_or_nstar_range'][0])} in every permutation"))
ca = c['A4']['calendar_adjusted']
rows.append(('2024 comparison with a calendar-month spline', '–', f"operative time {f(ca['or_time']['estimate'], 1)} min (95% CI {f(ca['or_time']['ci95'][0], 1)} to {f(ca['or_time']['ci95'][1], 1)}), pT2 margin {f(ca['psm_pt2']['estimate'] * 100, 1)} points (95% CI {f(ca['psm_pt2']['ci95'][0] * 100, 1)} to {f(ca['psm_pt2']['ci95'][1] * 100, 1)})"))
rows.append(('Mature phases with plateau cutoffs re-selected in each bootstrap replicate (1000 replicates)', '–',
             f"plateaus reached on both robots in {b1['share_reached'] * 100:.0f}% of replicates. Among them, operative time 90% interval {f(b1['or_time']['p5'], 1)} to {f(b1['or_time']['p95'], 1)} min, pT2 margin 90% interval {f(b1['psm_pt2']['p5'] * 100, 1)} to {f(b1['psm_pt2']['p95'] * 100, 1)} points"))
fx = b1['posthoc_fixed_cutoff']
rows.append(('Post hoc, Versius cases 47–337 against all restart cases', '–',
             f"operative time {f(fx['or_time']['estimate'], 1)} min (90% CI {f(fx['or_time']['ci90'][0], 1)} to {f(fx['or_time']['ci90'][1], 1)}), pT2 margin {f(fx['psm_pt2']['estimate'] * 100, 1)} points (90% CI {f(fx['psm_pt2']['ci90'][0] * 100, 1)} to {f(fx['psm_pt2']['ci90'][1] * 100, 1)})"))
cs = json.load(open(R / 'B1f_cutoff_sensitivity.json'))
for k, first in (('after_30', 31), ('after_43', 44)):
    x = cs[k]
    rows.append((f'Post hoc, Versius cases {first}–337 against all restart cases', '–',
                 f"operative time {f(x['or_time']['estimate'], 1)} min (90% CI {f(x['or_time']['ci90'][0], 1)} to {f(x['or_time']['ci90'][1], 1)}), pT2 margin {f(x['psm_pt2']['estimate'] * 100, 1)} points (90% CI {f(x['psm_pt2']['ci90'][0] * 100, 1)} to {f(x['psm_pt2']['ci90'][1] * 100, 1)})"))
b3 = s['B3']
rows.append(('Previous robotic case on the other robot within 7 days, 2024', '–', f"{int(b3['n_switch'])} of {int(b3['n_total'])} operations, {f(b3['est_min'], 1)} min (95% CI {f(b3['ci95'][0], 1)} to {f(b3['ci95'][1], 1)})"))
cu = s['S_cusum']
rows.append(('Risk-adjusted CUSUM (odds ratio 2), positive margin at any stage', f"observed {cu['Versius_racusum']['observed']}, expected {f(cu['Versius_racusum']['expected'])}, signals {cu['Versius_racusum']['signals']}",
             f"observed {cu['dVrestart_racusum']['observed']}, expected {f(cu['dVrestart_racusum']['expected'])}, signals {cu['dVrestart_racusum']['signals']}" + (f", first at case {cu['dVrestart_racusum']['first_signal_case']}" if cu['dVrestart_racusum']['first_signal_case'] else '')))
lv, ld = cu['Versius_lccusum_pT2'], cu['dVrestart_lccusum_pT2']
rows.append(('LC-CUSUM, pT2 margin, 10% acceptable and 20% unacceptable', f"{'competence signal' if lv['competence_signal'] else 'no competence signal'} ({lv['events']}/{lv['n_pT2']}); achieved alpha {lv['achieved_alpha']:.3f}, beta {lv['achieved_beta']:.3f}",
             f"{'competence signal' if ld['competence_signal'] else 'no competence signal'} ({ld['events']}/{ld['n_pT2']}); achieved alpha {ld['achieved_alpha']:.3f}, beta {ld['achieved_beta']:.3f}"))
md = '| Analysis | Versius | da Vinci restart |\n|---|---|---|\n' + '\n'.join(f"| {a} | {b} | {c_} |" for a, b, c_ in rows)
rm = s['S_racusum_model']
md += (f"\n\nRisk model for the risk-adjusted CUSUM. {rm['n']} operations, {rm['events']} positive margins, apparent c-statistic {rm['auc_apparent']:.2f}, "
       f"optimism-corrected {rm['auc_optimism_corrected']:.2f}.")
md = md.replace(' -', ' −').replace('(-', '(−').replace('to -', 'to −')
(WS / 'tables/tableS1_sensitivity.md').write_text(md); print(md)
