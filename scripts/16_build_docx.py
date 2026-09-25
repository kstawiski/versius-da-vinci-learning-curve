#!/usr/bin/env python3
"""Build DOCX files: manuscript (tables inline, figures submitted separately) and supplement (figures embedded).
Usage: 16_build_docx.py MANUSCRIPT_ASSEMBLED.md SUPPLEMENT.md OUTDIR"""
import subprocess, sys, re
from pathlib import Path
from docx import Document
from docx.shared import Pt
ms, sup, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
def build(src, dst, resource):
    subprocess.run(['pandoc', str(src), '-f', 'markdown-smart', '-t', 'docx', '-o', str(dst), '--resource-path', str(resource)], check=True)
    d = Document(dst)
    for st in ['Normal', 'Body Text', 'First Paragraph', 'Compact']:
        if st in [s.name for s in d.styles]:
            f = d.styles[st].font; f.name = 'Times New Roman'; f.size = Pt(12)
            pf = d.styles[st].paragraph_format; pf.line_spacing = 2.0 if st != 'Compact' else 1.0
    d.save(dst)
build(ms, out / 'Manuscript.docx', ms.parent)
build(sup, out / 'Supplementary_material.docx', sup.parent)
for f in ['fig1_series_structure', 'fig2_learning_curves', 'fig3_contrasts']:
    src = ms.parent.parent / 'figures' / f'{f}.tiff'
    n = {'fig1_series_structure': 'Figure1', 'fig2_learning_curves': 'Figure2', 'fig3_contrasts': 'Figure3'}[f]
    (out / f'{n}.tiff').write_bytes(src.read_bytes())
print('built', sorted(p.name for p in out.iterdir()))
