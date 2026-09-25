# Code and aggregate outputs

From Versius to da Vinci: learning curve and skill transfer in 738 single-surgeon prostatectomies.

The archive contains analysis code and aggregate outputs only. No patient-level rows, identifiers or source documents are included. The patient-level research database is restricted under institutional approval and cannot be shared.

## Contents

- `ANALYSIS_PLAN.md`: the statistical analysis plan, finalised before outcomes were examined by case order. Deviations are listed in Appendix S1 of the Supplementary material.
- `scripts/`: analysis scripts, run in numeric order by `scripts/run_all.sh` against the restricted database.
  - `01_derive.py`: analysis dataset, case order and learning-phase indices.
  - `lc_common.R`: shared model, imputation, standardization and contrast functions.
  - `03_primary.R`, `03b_sparse_binary.R`, `03c_mature_boot.R`, `03d_cutoff_sensitivity.R`: learning curves, plateau and segmented estimates, robot comparisons and cutoff sensitivity.
  - `04_secondary.R`, `04b_davies.R`, `06_posthoc_learning_metrics.R`, `07_team_curves.R`: sensitivity analyses, CUSUM charts, post hoc learning metrics and hospital-team curves.
  - `05_*.R`: figures. `08`–`11` and `09`: tables and the authoritative contrast table. `12`–`14`, `16`, `17`: references, manuscript assembly, number binding and document build.
- `analysis/results/`: aggregate results written by the scripts (JSON and CSV).
- `figures/source_data/`: aggregate data behind each figure.
- `tables/`: the rendered tables.

## Environment

R 4.6.1 with mgcv 1.9.4, mice 3.19.0, segmented 2.2.2, sandwich 3.1.2, logistf 1.26.1, ggplot2 4.0.3 and patchwork 1.3.2. Python 3.14.3 with pandas 2.3.3, NumPy 2.5.2 and statsmodels 0.14.6.

## Reproduction

With access to the restricted database, run `bash scripts/run_all.sh` from the archive root after pointing `01_derive.py` at the database release. Figures and tables can be rebuilt from `analysis/results/` and `figures/source_data/` alone. Multiple imputation and bootstrap estimates reproduce within Monte Carlo error because they depend on fixed random seeds and the installed package versions.

## Licence

MIT, see `LICENSE`.
