import os
import re

target_dir = "<PROJECT_ROOT>/manuscript6/analysis/repro_results"
iso_date_pattern = re.compile(r'\b(19\d\d|20\d\d)-(0[1-9]|1[0-2])-(0[1-9]|[12]\d|3[01])\b')
id_pattern = re.compile(r'\b(PL-\d+|patient_\d+|case_\d{4,})\b', re.IGNORECASE)

print(f"Scanning directory for privacy compliance: {target_dir}")
violating_files = []

for root, dirs, files in os.walk(target_dir):
    for f in files:
        fpath = os.path.join(root, f)
        # We allow reference to raw database folder path or dates like 2026-09-04 / 2026-09-25 in folder names / seeds,
        # but check for any patient date or row-level identifier
        with open(fpath, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
        for idx, line in enumerate(lines, 1):
            # Exclude known release names/paths or script seeds
            clean_line = line.replace("2026-09-04", "").replace("2026-09-25", "").replace("20260925", "")
            m_date = iso_date_pattern.search(clean_line)
            m_id = id_pattern.search(clean_line)
            if m_date or m_id:
                violating_files.append((fpath, idx, m_date.group(0) if m_date else m_id.group(0)))

if not violating_files:
    print("PASS: Zero patient dates or patient identifiers detected outside restricted directory.")
else:
    print(f"FAIL: Detected {len(violating_files)} potential privacy issues:")
    for vf in violating_files:
        print(f"  {vf[0]}:{vf[1]} -> {vf[2]}")
