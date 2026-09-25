#!/usr/bin/env python3
"""Build the numbered Vancouver reference list (BJUI style) from Crossref metadata for the cited keys,
in order of first citation. Writes references/REFERENCES_v2.md and references/references_v2.json."""
import csv, json, sys, time, urllib.request, urllib.parse, html, re
from pathlib import Path
WS = Path(__file__).resolve().parents[1]
ORDER = ['FernandezConejo2026BJU', 'Almeida2025', 'Chen2025', 'Reitano2025', 'Chierigo2026', 'Ficarra2024', 'Hammadeh2026',
         'Antonelli2025Minerva', 'Yu2026', 'Gavi2025', 'Liu2025', 'Roberts2025', 'vonElm2007', 'McCulloch2009', 'Dindo2004',
         'Wood2011', 'Duan1983', 'vanBuuren2011', 'Muggeo2003', 'Firth1993', 'Heinze2002', 'Steiner2000', 'Biau2008',
         'Newcombe1998', 'Bravi2024PSM', 'Vickers2007', 'Meneghetti2026', 'Mangano2021']
doi = {}
for f in ['evidence/literature/references.csv', 'evidence/literature/references_methods.csv']:
    for r in csv.DictReader(open(WS / f)): doi[r['key']] = r['doi']
out, lines = [], []
for i, k in enumerate(ORDER, 1):
    d = doi[k]
    req = urllib.request.Request('https://api.crossref.org/works/' + urllib.parse.quote(d, safe=''), headers={'User-Agent': 'ms6-refs (mailto:none@example.org)'})
    m = json.load(urllib.request.urlopen(req, timeout=40))['message']
    au = m.get('author', [])
    names = []
    for a in au:
        fam = a.get('family') or a.get('name', '')
        ini = ''.join(p[0] for p in re.split(r'[\s\-\.]+', a.get('given', '')) if p)
        names.append(f"{fam} {ini}".strip())
    auth = ', '.join(names) if len(names) <= 6 else ', '.join(names[:3]) + ', et al'
    title = html.unescape(re.sub(r'<[^>]+>', '', (m.get('title') or [''])[0])).strip().rstrip('.')
    jn = (m.get('short-container-title') or m.get('container-title') or [''])[0]
    yr = (m.get('published-print') or m.get('published-online') or m.get('issued'))['date-parts'][0][0]
    vol, iss, pg = m.get('volume'), m.get('issue'), m.get('page') or m.get('article-number')
    loc = f"{yr}"
    if vol: loc += f";{vol}"
    if iss: loc += f"({iss})"
    if pg: loc += f":{pg}"
    ref = f"{auth}. {title}. {jn} {loc}. doi:{d}"
    out.append({'n': i, 'key': k, 'doi': d, 'reference': ref})
    lines.append(f"{i}. {ref}")
    time.sleep(0.2)
JMAP = {'BJU International': 'BJU Int', 'Int. braz j urol.': 'Int Braz J Urol', 'J Robotic Surg': 'J Robot Surg', 'European Urology Open Science': 'Eur Urol Open Sci',
        'The Lancet': 'Lancet', 'Annals of Surgery': 'Ann Surg', 'Journal of the Royal Statistical Society Series B: Statistical Methodology': 'J R Stat Soc Series B Stat Methodol',
        'Journal of the American Statistical Association': 'J Am Stat Assoc', 'J. Stat. Soft.': 'J Stat Softw', 'Statistics in Medicine': 'Stat Med', 'Statist. Med.': 'Stat Med',
        'British Journal of Surgery': 'Br J Surg', 'European Urology Oncology': 'Eur Urol Oncol', 'Robotics Computer Surgery': 'Int J Med Robot',
        'JNCI Journal of the National Cancer Institute': 'J Natl Cancer Inst'}
AFIX = {'ANTONELLI A, VECCIA A, MALANDRA S': 'Antonelli A, Veccia A, Malandra S', 'Buuren Sv, Groothuis-Oudshoorn K': 'van Buuren S, Groothuis-Oudshoorn K',
        'FIRTH D': 'Firth D', 'Chen Sy, Ma Y, Yang Jw, Liu H, Wang L, Li Xr': 'Chen SY, Ma Y, Yang JW, Liu H, Wang L, Li XR', 'Fernández‐Conejo': 'Fernández-Conejo', 'Steiner SH. Monitoring surgical': 'Steiner SH, Cook RJ, Farewell VT, Treasure T. Monitoring surgical'}  # Crossref lists one author; PubMed 12933566 lists four
TFIX = {'inR': 'in R', 'break‐points': 'break-points', 'Versius‐assisted': 'Versius-assisted', '18‐Month of Follow‐Up': '18-Month of Follow-Up'}
lines = []
for o in out:
    r = o['reference']
    for a, b in list(AFIX.items()) + list(TFIX.items()): r = r.replace(a, b)
    for a, b in JMAP.items(): r = r.replace('. ' + a + ' ', '. ' + b + ' ')
    r = r.replace('BJU Int 2026:bju.70405', 'BJU Int 2026').replace('Minerva Urol Nephrol 2025;77(5).', 'Minerva Urol Nephrol 2025;77(5).')
    o['reference'] = r; lines.append(f"{o['n']}. {r}")
(WS / 'references/references_v2.json').write_text(json.dumps(out, indent=1, ensure_ascii=False))
(WS / 'references/REFERENCES_v2.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
