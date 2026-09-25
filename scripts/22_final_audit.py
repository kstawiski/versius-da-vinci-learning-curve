#!/usr/bin/env python3
"""Final technical audit ledgers on the exact submission package and its sources (manuscript-technicalities).
Covers prose mechanics, abbreviation ledger, vocabulary screen, pronouns, reference list rendering,
citation ledger, callouts and panels, P-value and interval style, DOCX and PDF integrity, TIFF properties.
Usage: 22_final_audit.py MANUSCRIPT_ASSEMBLED.md SUPPLEMENT.md PACKAGE_DIR OUT.json"""
import json, re, subprocess, sys, zipfile
from pathlib import Path
ms, sup, pkg, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
WS = Path(__file__).resolve().parents[1]
T = ms.read_text(); S = sup.read_text(); C = (WS / 'submission/cover-letter.md').read_text()
R = {}

def prose_blocks(text, skip=('References', 'Tables')):
    sec, blocks = None, []
    for ln in text.split('\n'):
        if ln.startswith('## '): sec = ln[3:].strip(); continue
        if ln.startswith('#') or not ln.strip() or ln.startswith('|') or ln.startswith('!['): continue
        if sec in skip: continue
        if re.match(r'^(Running title|Article type|Keywords|Author contributions)\.', ln) or ln.startswith('^') or ln.startswith('**'): continue
        blocks.append((sec or 'front', ln))
    return blocks

def sentences(ln):
    ln = re.sub(r'\b(e\.g|i\.e|vs|et al|Fig|No)\.', lambda m: m.group(0).replace('.', '§'), ln)
    return [s.replace('§', '.') for s in re.split(r'(?<=[.!?])\s+(?=[A-Z(\[])', ln)]

# 1 prose mechanics
mech = {'long_35_50': [], 'long_over_50': [], 'punct': {}}
for name, text in [('manuscript', T), ('supplement', S), ('cover_letter', C)]:
    words = 0; counts = {'colon': 0, 'semicolon': 0, 'em_dash': 0}
    for sec, ln in prose_blocks(text):
        body = re.sub(r'https?://\S+|doi:\S+|\[PLACEHOLDER[^\]]*\]', '', ln)
        body = re.sub(r'\d:\d', '', body)
        words += len(body.split())
        counts['colon'] += len(re.findall(r':', body)); counts['semicolon'] += body.count(';'); counts['em_dash'] += body.count('—')
        for s in sentences(ln):
            n = len(s.split())
            if n > 50: mech['long_over_50'].append((name, sec, n, s[:100]))
            elif n > 35: mech['long_35_50'].append((name, sec, n, s[:100]))
    mech['punct'][name] = {**counts, 'words': words, 'per_1000': {k: round(1000 * v / max(words, 1), 2) for k, v in counts.items()}}
R['mechanics'] = mech

# 2 abbreviation ledger (abstract and body separately)
abstract = T.split('## Abstract', 1)[1].split('Keywords.', 1)[0]
body = T.split('## Introduction', 1)[1].split('## References', 1)[0]
def abbrev_ledger(text):
    led = {}
    for m in re.finditer(r'\b([A-Z][A-Z0-9-]{1,}[a-z]?)\b', text):
        a = m.group(1)
        if a in led or a in {'CI', 'P', 'USA', 'UK', 'CA', 'MA', 'II', 'III', 'IV', 'IIIa', 'IIIb', 'NA'} or re.fullmatch(r'[IVX]+', a): continue
        start = m.start(); ctx = text[max(0, start - 150):start + 60]
        defined = bool(re.search(r'\(' + re.escape(a) + r'[\),]', text[:start + len(a) + 2])) or bool(re.search(re.escape(a) + r'\s*\(', text[start:start + len(a) + 3]))
        led[a] = {'first_pos': start, 'defined_at_or_before_first_use': defined, 'context': ctx.replace('\n', ' ')[-120:]}
    return led
R['abbreviations_abstract'] = {k: v['defined_at_or_before_first_use'] for k, v in abbrev_ledger(abstract).items()}
R['abbreviations_body'] = {k: v['defined_at_or_before_first_use'] for k, v in abbrev_ledger(body).items()}

# 3 vocabulary screen and pronouns
VOC = r'\b(null|attenuated toward|bounded|canonical|ascertainment effect|variance expansion|quantile|cross-tabulat|Bernoulli|robust|compelling|noteworthy|pivotal|nuanced|underscore|sheds? light|testament|leverage|utili[sz]e|operationali[sz]e|landscape|paradigm|multifaceted|holistic|tapestry|delve|evidence ceiling|provenance|upstream|downstream|it is important to note|taken together|in the context of|could potentially)\b'
R['vocabulary_hits'] = [(n, sec, re.search(VOC, ln, re.I).group(0), ln[:110]) for n, text in [('manuscript', T), ('supplement', S), ('cover_letter', C)]
                        for sec, ln in prose_blocks(text) if re.search(VOC, ln, re.I)]
R['gendered_pronouns'] = [(n, ln[:110]) for n, text in [('manuscript', T), ('supplement', S), ('cover_letter', C)]
                          for sec, ln in prose_blocks(text) if re.search(r'\b(he|his|him|she|her|hers)\b', ln, re.I)]

# 4 reference list rendering
refs = [r for r in T.split('## References', 1)[1].split('## Tables', 1)[0].strip().split('\n') if r.strip()]
rl = []
for r in refs:
    n = int(r.split('.', 1)[0]); body_r = r.split('. ', 1)[1]
    authors = body_r.split('. ', 1)[0]
    n_auth = len([a for a in authors.split(', ') if a and a != 'et al'])
    rl.append({'n': n, 'n_authors_listed': n_auth, 'et_al': 'et al' in authors, 'rule_ok': (n_auth <= 6 and 'et al' not in authors) or (n_auth == 3 and 'et al' in authors),
               'has_doi': 'doi:' in r, 'has_year': bool(re.search(r'\b(19|20)\d\d\b', r)), 'has_volume': bool(re.search(r';\d+', r)), 'text': r[:140]})
R['references'] = {'count': len(refs), 'rule_violations': [x for x in rl if not x['rule_ok']], 'no_doi': [x['n'] for x in rl if not x['has_doi']],
                   'no_volume': [x['text'] for x in rl if not x['has_volume']], 'doi_duplicates': len(refs) - len({re.search(r'doi:(\S+)', r).group(1) for r in refs if 'doi:' in r})}
cites = []
for m in re.finditer(r'\[(\d+(?:[–,\s]+\d+)*)\]', T.split('## References')[0]):
    for part in re.split(r',\s*', m.group(1)):
        if '–' in part: a, b = map(int, part.split('–')); cites += list(range(a, b + 1))
        else: cites.append(int(part))
first = list(dict.fromkeys(cites))
R['citations'] = {'first_occurrence_order_ok': first == list(range(1, len(first) + 1)), 'uncited': sorted(set(range(1, len(refs) + 1)) - set(cites)),
                  'markers_without_entry': sorted(set(cites) - set(range(1, len(refs) + 1)))}

# 5 callouts with panels
main = T.split('## Introduction', 1)[1].split('## References', 1)[0]
R['panel_callouts'] = sorted(set(re.findall(r'Figure (\d[A-D])', main)))
legends = T.split('## Figure legends', 1)[1].split('## Supplementary material', 1)[0]
R['legend_panels'] = {f: sorted(set(re.findall(r'\(([A-D])\)', blk))) for f, blk in re.findall(r'\*\*Figure (\d)\.\*\*(.*?)(?=\*\*Figure|\Z)', legends, re.S)}

# 6 P values and intervals
R['p_value_forms'] = sorted(set(re.findall(r'\bP\s*[=<>]\s*[0-9.]+|\bp\s*[=<>]\s*[0-9.]+', T + S)))
R['exact_zero_p'] = re.findall(r'P\s*=\s*0(?:\.0+)?\b(?!\.\d*[1-9])', T + S)
R['hyphen_as_minus'] = [m.group(0) for m in re.finditer(r'(?<![\w-])-\d', T.split('## References')[0])][:10]

# 7 DOCX integrity
dx = {}
for f in sorted(pkg.glob('*.docx')):
    z = zipfile.ZipFile(f); names = z.namelist(); x = z.read('word/document.xml').decode(errors='ignore')
    core = z.read('docProps/core.xml').decode(errors='ignore') if 'docProps/core.xml' in names else ''
    dx[f.name] = {'tracked_changes': bool(re.search(r'<w:(ins|del) ', x)), 'comments': 'word/comments.xml' in names,
                  'line_numbers': 'lnNumType' in x, 'placeholder_style': x.count('PlaceholderText'),
                  'creator': (re.search(r'<dc:creator>(.*?)</dc:creator>', core) or [None, ''])[1],
                  'last_modified_by': (re.search(r'<cp:lastModifiedBy>(.*?)</cp:lastModifiedBy>', core) or [None, ''])[1],
                  'text_placeholders': len(re.findall(r'\[PLACEHOLDER', re.sub(r'<[^>]+>', '', x)))}
R['docx'] = dx

# 8 PDF integrity
pdf = {}
for f in sorted(pkg.glob('*.pdf')):
    fonts = subprocess.run(['pdffonts', str(f)], capture_output=True, text=True).stdout.split('\n')[2:]
    txt = subprocess.run(['pdftotext', str(f), '-'], capture_output=True, text=True).stdout
    info = subprocess.run(['pdfinfo', str(f)], capture_output=True, text=True).stdout
    pdf[f.name] = {'pages': int(re.search(r'Pages:\s+(\d+)', info).group(1)), 'fonts_not_embedded': [l for l in fonts if l.strip() and ' no ' in l[:80]],
                   'replacement_chars': txt.count('\ufffd'), 'placeholder_text': txt.count('[PLACEHOLDER')}
R['pdf'] = pdf

# 9 TIFF properties
from PIL import Image
tif = {}
for f in sorted(pkg.glob('*.tiff')):
    im = Image.open(f); dpi = im.info.get('dpi')
    tif[f.name] = {'pixels': im.size, 'dpi': tuple(round(v) for v in dpi) if dpi else None,
                   'width_mm': round(im.size[0] / dpi[0] * 25.4, 1) if dpi else None, 'compression': im.info.get('compression'), 'mode': im.mode}
R['tiff'] = tif
out.write_text(json.dumps(R, ensure_ascii=False, indent=1, default=str))
print(json.dumps({k: (v if k not in ('abbreviations_abstract', 'abbreviations_body') else {a: d for a, d in v.items() if not d})
                  for k, v in R.items()}, ensure_ascii=False, indent=1, default=str)[:9000])
