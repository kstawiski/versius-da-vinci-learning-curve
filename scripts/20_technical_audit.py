#!/usr/bin/env python3
"""Deterministic technical audit ledgers for the assembled manuscript and supplement:
link resolution (DOIs and URLs), abstract-to-body number fidelity, table percentage arithmetic, display-item
callout inventory with first-mention order, and citation coverage per Introduction and Discussion paragraph.
Usage: 20_technical_audit.py MANUSCRIPT_ASSEMBLED.md SUPPLEMENT.md OUT.json"""
import json, re, sys, urllib.request, datetime
from pathlib import Path
ms, sup, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
T = ms.read_text(); S = sup.read_text()
res = {'date': datetime.date.today().isoformat(), 'manuscript': str(ms), 'supplement': str(sup)}

# 1. Link resolution
refs = T.split('## References', 1)[1].split('## Tables', 1)[0]
links = []
for m in re.finditer(r'doi:(10\.\S+?)(?=\.?\s|\.?$)', refs, re.M):
    links.append(('doi', m.group(1).rstrip('.')))
for m in re.finditer(r'https?://[^\s)]+', T.split('## References', 1)[0] + T.split('## Tables', 1)[1]):
    links.append(('url', m.group(0).rstrip('.')))
ledger = []
for kind, v in links:
    url = f'https://doi.org/{v}' if kind == 'doi' else v
    try:
        req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'Mozilla/5.0 (link check)'})
        r = urllib.request.urlopen(req, timeout=30); status, final = r.status, r.geturl()
    except urllib.error.HTTPError as e:
        status, final = e.code, getattr(e, 'url', url)
    except Exception as e:
        status, final = type(e).__name__, ''
    ledger.append({'kind': kind, 'id': v, 'status': status, 'final_url': final[:150]})
res['links'] = ledger

# 2. Abstract-to-body fidelity: every number in the abstract must appear in the body, tables or legends.
abstract = T.split('## Abstract', 1)[1].split('Keywords.', 1)[0]
body = T.split('## Introduction', 1)[1]
nums = sorted(set(re.findall(r'−?\d+(?:\.\d+)?', abstract)), key=lambda x: float(x.replace('−', '-')))
res['abstract_numbers_missing_from_body'] = [n for n in nums if n not in body]

# 3. Table arithmetic: n/N (%) cells recomputed with rounding tolerance.
bad = []
for label, text in [('main', T.split('## Tables', 1)[1].split('## Figure legends', 1)[0]), ('supplement', S)]:
    for m in re.finditer(r'(\d+)/(\d+) \((\d+(?:\.\d+)?)%\)', text):
        n, d, p = int(m.group(1)), int(m.group(2)), float(m.group(3))
        if d == 0 or n > d or abs(100 * n / d - p) > 0.51:
            bad.append({'where': label, 'cell': m.group(0)})
for m in re.finditer(r'(\d+) of (\d+) [^()]{0,40}\((\d+)%\)', body):
    n, d, p = int(m.group(1)), int(m.group(2)), float(m.group(3))
    if n > d or abs(100 * n / d - p) > 0.51:
        bad.append({'where': 'text', 'cell': m.group(0)})
res['percentage_errors'] = bad

# 4. Callout inventory and first-mention order (main text before the reference list).
text_main = T.split('## Introduction', 1)[1].split('## References', 1)[0]
inv = {}
for kind, pat in [('Figure', r'Figures? (\d)(?:[A-D])?(?:\s*(?:and|–)\s*(\d))?'), ('Table', r'Tables? (\d)(?:\s*(?:and|–)\s*(\d))?'),
                  ('Table S', r'Table S(\d)'), ('Figure S', r'Figures? S(\d)'), ('Appendix S', r'Appendix S(\d)')]:
    order = []
    for m in re.finditer(pat, text_main):
        for g in m.groups():
            if g and g not in order: order.append(g)
    inv[kind] = order
res['callouts_first_mention_order'] = inv
res['supplied_items'] = {'Figure': ['1', '2', '3'], 'Table': ['1', '2', '3'], 'Table S': ['1'], 'Figure S': ['1', '2'], 'Appendix S': ['1', '2']}

# 5. Citation coverage per paragraph in Introduction and Discussion.
cov = []
for sec in ['Introduction', 'Discussion']:
    blk = T.split(f'## {sec}', 1)[1].split('\n## ', 1)[0]
    paras = [p.strip() for p in blk.split('\n\n') if p.strip() and not p.strip().startswith('#')]
    for i, p in enumerate(paras, start=1):
        c = len(re.findall(r'\[\d+(?:[–,\s]+\d+)*\]', p))
        cov.append({'section': sec, 'paragraph': i, 'of': len(paras), 'citations': c, 'opening': p[:90]})
res['citation_coverage'] = cov
out.write_text(json.dumps(res, ensure_ascii=False, indent=1))
print(json.dumps({k: (v if k not in ('links', 'citation_coverage') else None) for k, v in res.items()}, ensure_ascii=False, indent=1))
print('links not 200:', [l for l in ledger if l['status'] != 200])
print('zero-citation paragraphs:', [(c['section'], c['paragraph'], c['of'], c['opening']) for c in cov if c['citations'] == 0])
