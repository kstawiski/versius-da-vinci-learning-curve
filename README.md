# Code and aggregate outputs

From Versius to da Vinci: operative time and surgical margins in 738 single-surgeon prostatectomies. Version 1.1.3.

The archive contains analysis code and aggregate outputs only. No patient-level rows, identifiers or source documents are included. The patient-level research database is restricted under institutional approval and cannot be shared.

## Contents

- `ANALYSIS_PLAN.md`: the statistical analysis plan (version 2), locked before outcomes were examined by case order and after aggregate case mix and completeness had been reviewed. `plan/SAP_CHANGELOG.md` records every change from version 1, the implementation notes and the post hoc additions, and governs where the two files differ. Deviations are also listed in Appendix S1 of the Supplementary material.
- `scripts/`: analysis scripts, run in numeric order by `scripts/run_all.sh` against the restricted database.
  - `01_derive.py`: analysis dataset, case order and learning-phase indices.
  - `lc_common.R`: shared model, imputation, standardization and contrast functions.
  - `03_primary.R`: learning curves, plateaus, first-100 comparisons, first-50-case slope models (A3a), cumulative versus platform-specific experience (A3c) and the 2024 comparison.
  - `03b`–`03f`: sparse-outcome models, the mature-phase bootstrap, cutoff sensitivity, smooth-term tests pooled across imputations (D2 rule) and the 2024 overlap and ascertainment analysis.
  - `04_secondary.R`, `04b_davies.R`, `06_posthoc_learning_metrics.R`, `07_team_curves.R`: sensitivity analyses, CUSUM charts, post hoc learning metrics and hospital-specific curves.
  - `05_*.R`: figures. `08`–`11` and `21`: tables. `12`–`14`, `16`, `17`, `19`, `20`: references, manuscript assembly, number binding, document build and technical audit.
- `analysis/results/`: aggregate results written by the scripts (JSON and CSV).
- `figures/source_data/`: aggregate data behind the figures.
- `tables/`: the rendered tables.
- `reproduction/`: the independent re-implementation (code and aggregate results). Absolute local paths in this code were replaced by `<PROJECT_ROOT>`; nothing else was changed.

## Environment

R 4.6.1 with mgcv 1.9.4, mice 3.19.0, segmented 2.2.2, sandwich 3.1.2, logistf 1.26.1, ggplot2 4.0.3 and patchwork 1.3.2. Python 3.14.3 with pandas 2.3.3, NumPy 2.5.2 and statsmodels 0.14.6.

## Reproduction

With access to the restricted database, run `bash scripts/run_all.sh` from the archive root after pointing `01_derive.py` at the database release. Figure 3 and the tables can be checked against `analysis/results/` and `figures/source_data/`. Figures 1, 2 and S1 are drawn from the restricted patient-level data and cannot be rebuilt from this archive alone. Multiple imputation and bootstrap estimates reproduce within Monte Carlo error because they depend on fixed random seeds and the installed package versions.

## Licence

MIT, see `LICENSE`.
