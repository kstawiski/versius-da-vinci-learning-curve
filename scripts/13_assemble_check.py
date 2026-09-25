#!/usr/bin/env python3
"""Assemble the manuscript (references, tables) and run deterministic checks:
prose mechanics (no colon, semicolon or em dash in prose; sentences <= 40 words), citation order,
word counts against BJUI limits. Usage: 13_assemble_check.py SOURCE.md OUTPUT.md"""
import re, sys, json
from pathlib import Path
WS = Path(__file__).resolve().parents[1]
src, out = Path(sys.argv[1]), Path(sys.argv[2])
t = src.read_text()
refs = (WS / 'references/REFERENCES_v2.md').read_text().strip()
def dash(s): return re.sub(r'(?<=\d)-(?=\d)', '–', s)
T1 = dash((WS / 'tables/table1_case_mix.md').read_text().strip())
T2 = dash((WS / 'tables/table2_outcomes.md').read_text().strip())
T3 = dash((WS / 'tables/table3_margins_hospital_period.md').read_text().strip())
pc = json.load(open(WS / 'analysis/results/table3_posthoc_contrasts.json'))
def ci(k):
    x = pc[k]; f = lambda v: f"{100 * v:.1f}".replace('-', '−')
    return f"{f(x['diff'])} points (95% CI {f(x['ci95'][0])} to {f(x['ci95'][1])})"
k = list(pc.keys())
tables = f"""**Table 1.** Patient and tumour characteristics by robotic series and, for 2024, by robot.

{T1}

Values are median (IQR) or n/N (%). Denominators exclude missing values. The da Vinci restart comprises da Vinci cases 20–401, operated from 2024. Da Vinci cases 1–19 were operated in 2021–2023. EAU, European Association of Urology. ISUP, International Society of Urological Pathology.

**Table 2.** Perioperative, pathological and functional outcomes by series and learning phase.

{T2}

Values are median (IQR) or n/N (%). Phases are the surgeon's case numbers on Versius and on the da Vinci restart. Denominators exclude missing values. Clavien-Dindo grades below III are not shown because the registry stopped recording them routinely in 2024. Grades are as recorded in the registry, and two of the six Versius reoperations carried a grade below III. Blood loss was mostly recorded as rounded estimates and is shown descriptively. PSA persistence includes men operated at least 56 days before data close. Continence at 3 and 12 months includes men operated at least 120 and 425 days before data close, respectively. A dash marks a phase in which no patient had reached the 12-month window.

**Table 3.** Positive surgical margins by hospital, robot and calendar period.

{T3}

Values are n/N (%). The column headed structured margin report gives the share of pathology reports that recorded the margin length in a structured field. Post hoc comparisons in pT2 disease used 95% CIs by the Newcombe method, and a positive difference means a higher rate in the group named first in each comparison. At SalveMedica, the da Vinci restart from January 2024 to June 2025 against Versius in 2021–2023 differed by {ci(k[0])}, and the restart from July 2025 onward by {ci(k[1])}. At Bełchatów, Versius in 2024 against Versius in 2021–2023 differed by {ci(k[2])}, and the whole restart against Versius in 2024 by {ci(k[3])}."""
t = t.replace('REFERENCES_PLACEHOLDER', refs).replace('TABLES_PLACEHOLDER', tables)
out.write_text(t)

# ---------- checks ----------
lines = t.split('\n')
prose, section = [], None
skip_sections = {'References', 'Tables'}
for ln in lines:
    if ln.startswith('## '): section = ln[3:].strip(); continue
    if ln.startswith('#') or not ln.strip() or ln.startswith('|'): continue
    if section is None or section in skip_sections: continue  # title-page apparatus and excluded surfaces
    if re.match(r'^(Running title|Article type|Authors|Affiliations|Corresponding author|Keywords|Author contributions)\.', ln): continue  # title-page apparatus and the CRediT role list
    prose.append((section, ln))
issues = []
for sec, ln in prose:
    body = re.sub(r'\[OWNER INPUT[^\]]*\]', '', ln)
    body = re.sub(r'https?://\S+|doi:\S+', '', body)
    if re.search(r'(?<!\d):(?!\d)', body): issues.append(('colon', sec, ln[:90]))
    if ';' in body: issues.append(('semicolon', sec, ln[:90]))
    if '—' in body: issues.append(('em dash', sec, ln[:90]))
    for s in re.split(r'(?<=[.!?])\s+(?=[A-Z(\[])', body.replace('**', '')):
        w = len(s.split())
        if w > 40: issues.append((f'sentence {w} words', sec, s[:90]))
# citation order
body_text = t.split('## References')[0]
cites = []
for m in re.finditer(r'\[(\d+(?:[–,\s]+\d+)*)\]', body_text):
    for part in re.split(r',\s*', m.group(1)):
        if '–' in part: a, b = map(int, part.split('–')); cites += list(range(a, b + 1))
        else: cites.append(int(part))
first = []
for c in cites:
    if c not in first: first.append(c)
order_ok = first == sorted(first) and first == list(range(1, len(first) + 1))
nref = len(refs.split('\n'))
def words(sec_start, sec_end):
    seg = t.split(sec_start)[1].split(sec_end)[0]
    seg = re.sub(r'^#+.*$', '', seg, flags=re.M)
    return len(re.sub(r'\[\d+(?:[–,\s]+\d+)*\]', '', seg).split())
abstract_words = words('## Abstract', 'Keywords.')
main_words = words('## Introduction', '## Declarations')
title = lines[0].lstrip('# ').strip()
print(json.dumps({'title_chars': len(title), 'abstract_words': abstract_words, 'main_text_words': main_words,
                  'references': nref, 'citation_order_ok': order_ok, 'first_citation_order': first,
                  'max_cited': max(cites) if cites else 0, 'mechanics_issues': len(issues)}, indent=1))
for i in issues: print('ISSUE', i)
