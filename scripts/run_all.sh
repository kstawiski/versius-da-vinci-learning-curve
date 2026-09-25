#!/bin/bash
# Full manuscript6 pipeline (aggregate outputs only). Run from the manuscript6 directory.
set -e
cd "$(dirname "$0")/.."
python3 scripts/00_feasibility_profile.py
python3 scripts/01_derive.py
python3 scripts/02_eda.py > /dev/null
MS6_NBOOT=${MS6_NBOOT:-2000} MS6_MIMP=20 MS6_NCORE=${MS6_NCORE:-16} Rscript scripts/03_primary.R
Rscript scripts/03b_sparse_binary.R
MS6_NBOOT_B1=${MS6_NBOOT_B1:-1000} MS6_NCORE=${MS6_NCORE:-16} Rscript scripts/03c_mature_boot.R
Rscript scripts/03d_cutoff_sensitivity.R
Rscript scripts/03e_pooled_smooth_tests.R
Rscript scripts/03f_overlap_2024.R
MS6_MIMP=20 MS6_NCORE=${MS6_NCORE:-16} MS6_NSEGBOOT=${MS6_NSEGBOOT:-500} MS6_NPERM=${MS6_NPERM:-200} Rscript scripts/04_secondary.R
Rscript scripts/04b_davies.R
MS6_NBOOT_PH=${MS6_NBOOT_PH:-1000} MS6_NCORE=${MS6_NCORE:-16} Rscript scripts/06_posthoc_learning_metrics.R
Rscript scripts/07_team_curves.R
python3 scripts/08_tables.py > /dev/null
python3 scripts/10_table3_margins.py > /dev/null
python3 scripts/11_authoritative_contrasts.py > /dev/null
python3 scripts/09_table_sensitivity.py > /dev/null || echo "TABLE_S1_NEEDS_UPDATE"
python3 scripts/21_table_s2_2024_months.py > /dev/null
Rscript scripts/05_fig1_timeline.R
Rscript scripts/05_fig2_curves.R
Rscript scripts/05_fig3_contrasts.R
Rscript scripts/05_figS_cusum.R
if [ -d manuscript ]; then  # manuscript assembly, only inside the manuscript workspace
  python3 scripts/13_assemble_check.py manuscript/MANUSCRIPT_v17.md manuscript/MANUSCRIPT_v17_assembled.md
  python3 scripts/14_verify_numbers.py manuscript/MANUSCRIPT_v17_assembled.md | tail -1
  python3 scripts/17_assemble_supplement.py manuscript/SUPPLEMENT_v7_source.md manuscript/SUPPLEMENT_v7.md
fi
echo PIPELINE_COMPLETE
