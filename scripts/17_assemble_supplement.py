#!/usr/bin/env python3
"""Assemble the supplement from its source: Table S1 from tables/tableS1_sensitivity.md and the reproduction
table from tables/tableS_repro.md. Usage: 17_assemble_supplement.py SOURCE.md OUTPUT.md"""
import sys
from pathlib import Path
WS = Path(__file__).resolve().parents[1]
src, out = Path(sys.argv[1]), Path(sys.argv[2])
t = src.read_text()
assert t.count('TABLE_S1_PLACEHOLDER') == 1 and t.count('REPRO_TABLE_PLACEHOLDER') == 1
if 'TABLE_S2_PLACEHOLDER' in t:
    t = t.replace('TABLE_S2_PLACEHOLDER', (WS / 'tables/tableS2_2024_by_month.md').read_text().strip())
t = t.replace('TABLE_S1_PLACEHOLDER', (WS / 'tables/tableS1_sensitivity.md').read_text().strip())
t = t.replace('REPRO_TABLE_PLACEHOLDER', (WS / 'tables/tableS_repro.md').read_text().strip())
out.write_text(t)
print('supplement assembled', out)
