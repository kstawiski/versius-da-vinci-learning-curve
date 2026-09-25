# Independent Reproduction Report: Surgical Learning-Curve Analysis (Manuscript 6)

## Progress and Reproduction Status
- [x] Initial environment verification: R 4.6.1 and Python 3.12 verified with all statistical packages (`mgcv`, `mice`, `sandwich`, `arrow`, `dplyr`, `logistf`, `pandas`, `numpy`, `statsmodels`, `scipy`).
- [x] Restricted intermediate environment secured: `<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/` initialized with permissions `0700` and patient-level files set to mode `0600`. Custom `TMPDIR` isolated within the restricted directory.
- [x] Case order and cohort derivation: 838 Polish operations derived according to the locked hierarchy; exact reconciliation of all R1 series counts (Versius 337, da Vinci 401, da Vinci restart 382, 2024 Versius 99, 2024 restart da Vinci 105).
- [x] Covariate imputation: MICE with predictive mean matching (PMM), 20 imputations,Auxiliary outcomes included in the imputation model but excluded from outcome estimation.
- [x] R2 Versius learning curves: Penalized regression splines (GAM REML, $k=10$ Gaussian for $\log(\text{or\_time})$ with smearing retransformation; $k=5$ Binomial for pT2 positive margin). Evaluated across 20 imputations.
- [x] R3 da Vinci restart learning curves: Identical GAM REML specifications applied to `dv_restart_n` ($N=382$).
- [x] R4 A3a transfer analysis: Restart da Vinci cases 1–100 versus Versius cases 1–100 for operative time (HC3 robust linear regression) and pT2 positive margins (marginally standardized logistic risk difference).
- [x] R5 A4 concurrent 2024 platform comparison: Adjusted differences for operative time, pT2 margins, all-stage margins, length of stay, 30-day readmission, EAU PSA persistence, and pad-free continence at 3 and 12 months with 95% and 90% confidence intervals.
- [x] R6 Bootstrap uncertainty: 600 case-resampling replicates of the Versius cohort with stochastic imputation to estimate 95% and 90% percentile intervals for plateau case number $n^*$.
- [x] Deliverables generated: `REPRO_RESULTS.json` serialized, code pipeline structured under `code/`, and privacy audit confirmed clean.

---

## Executive Summary of Reproduced Numbers

### R1. Cohort and Data Completeness
| Series | Total Cases | pT2 Cases ($n$ / $N$, %) | Operative Time ($n$ / $N$, %) |
| :--- | :---: | :---: | :---: |
| **Versius (All)** | 337 | 190 / 337 (56.4%) | 335 / 337 (99.4%) |
| **da Vinci (All)** | 401 | 253 / 401 (63.1%) | 401 / 401 (100.0%) |
| **da Vinci Restart** | 382 | 239 / 382 (62.6%) | 382 / 382 (100.0%) |
| **Versius 2024** | 99 | 50 / 99 (50.5%) | 99 / 99 (100.0%) |
| **da Vinci Restart 2024** | 105 | 63 / 105 (60.0%) | 105 / 105 (100.0%) |

### R2. Versius Learning Curves (Primary Efficiency and Quality Endpoints)
- **Operative Time Curve** (Gaussian GAM on log minutes, $k=10$, REML, Duan smearing factor):
  - **Asymptote** (mean over final 50 cases, 288–337): **159.68 minutes**
  - **$n^*$ Plateau Threshold (15-minute band)**: **Case 111** (standardized time reaches 174.65 min and remains $\le 174.68$ min thereafter).
  - **$n^*$ Plateau Threshold (10-minute band)**: **Case 337** (the curve does not stay within 10 minutes of the asymptote until the final case due to continued late-series downward drift).
  - **Slope over final 100 cases**: **$-0.2336$ minutes/case** (showing persistent gradual shortening).
  - **Standardized curve values at milestones**:
    - Case 1: **289.42 min**
    - Case 25: **200.65 min**
    - Case 50: **168.63 min**
    - Case 100: **177.86 min**
    - Case 337 (last): **146.15 min**
- **pT2 Positive Surgical Margin Curve** (Binomial GAM, $k=5$, REML, standardized over pT2 cases):
  - **Asymptote** (mean over final 50 cases): **34.58%**
  - **$n^*$ Plateau Threshold (5 percentage points band, $\le 39.58\%$)**: **Case 286**
  - **Standardized curve values at milestones**:
    - Case 1: **21.27%**
    - Case 25: **18.89%**
    - Case 50: **16.71%**
    - Case 100: **13.84%**
    - Case 337 (last): **39.23%**

### R3. da Vinci Restart Series Learning Curves
- **Operative Time Curve** (Gaussian GAM on log minutes, $k=10$, REML, Duan smearing factor):
  - **Asymptote** (mean over final 50 cases, 333–382): **157.62 minutes**
  - **$n^*$ Plateau Threshold (15-minute band)**: **Case 182**
  - **$n^*$ Plateau Threshold (10-minute band)**: **Case 196**
  - **Slope over final 100 cases**: **$+0.0317$ minutes/case** (stable plateau).
  - **Standardized curve values at milestones**:
    - Case 1: **158.76 min**
    - Case 25: **146.71 min**
    - Case 50: **145.81 min**
    - Case 100: **173.92 min**
    - Case 382 (last): **154.43 min**
- **pT2 Positive Surgical Margin Curve** (Binomial GAM, $k=5$, REML, standardized over pT2 cases):
  - **Asymptote** (mean over final 50 cases): **33.79%**
  - **$n^*$ Plateau Threshold (5 percentage points band, $\le 38.79\%$)**: **Case 382**
  - **Standardized curve values at milestones**:
    - Case 1: **33.20%**
    - Case 25: **35.93%**
    - Case 50: **38.46%**
    - Case 100: **39.39%**
    - Case 382 (last): **39.39%**

### R4. A3a Transfer of Proficiency: First 100 Cases (Restart da Vinci vs Versius)
- **Operative Time Difference** (linear regression on minutes, HC3 robust errors, pooled across 20 imputations):
  - **Adjusted Difference**: **$-33.09$ minutes** (SE: 7.88, $p < 0.0001$).
  - **95% CI**: **[$-48.54$, $-17.63$] minutes**
  - **90% CI**: **[$-46.05$, $-20.12$] minutes**
  - *Clinical reading*: The surgeon commenced the da Vinci restart series with an operative time advantage of 33 minutes compared to their initial 100 Versius cases, demonstrating immediate efficiency transfer from prior robotic experience.
- **pT2 Positive Margin Risk Difference** (logistic regression, marginally standardized risk difference, HC0 delta-method, pooled across 20 imputations):
  - **Adjusted Risk Difference**: **$+21.44$ percentage points** (SE: 9.47%, $p = 0.0236$).
  - **95% CI**: **[$+2.88\%$, $+40.01\%$]**
  - **90% CI**: **[$+5.86\%$, $+37.02\%$]**
  - *Clinical reading*: In contrast to operative efficiency, the initial 100 cases of the da Vinci restart exhibited a significantly higher risk of positive surgical margins in organ-confined disease than early Versius operations, reflecting differences in case selection or radicality during the resumption phase.

### R5. A4 Concurrent Platform Comparison (Calendar Year 2024)
Surgeon held constant, concurrent practice in both hospitals (99 Versius vs 105 da Vinci restart cases). Estimates represent da Vinci restart minus Versius.

| Endpoint | Crude Versius | Crude da Vinci Restart | Adjusted Estimate | 95% Confidence Interval | 90% Confidence Interval | $p$-value |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Operative Time (min)** | 174.99 (SD 47.52) | 152.18 (SD 40.53) | **$+6.47$ min** | [$-4.43$, $+17.37$] | [$-2.68$, $+15.62$] | 0.2446 |
| **pT2 Positive Margin** | 14 / 48 (29.2%) | 23 / 61 (37.7%) | **$+8.6$ pp** | [$-11.1\%$, $+28.2\%$] | [$-8.0\%$, $+25.1\%$] | 0.3943 |
| **All-Stage Margin** | 31 / 96 (32.3%) | 43 / 103 (41.7%) | **$+6.8$ pp** | [$-8.4\%$, $+21.9\%$] | [$-6.0\%$, $+19.5\%$] | 0.3826 |
| **Length of Stay (days)** | 4.51 (SD 1.07) | 4.25 (SD 0.63) | **$+0.21$ days** | [$-0.05$, $+0.46$] | [$-0.01$, $+0.42$] | 0.1147 |
| **30-Day Readmission** | 0 / 58 (0.0%) | 0 / 82 (0.0%) | **$0.00$ pp** | [$0.00\%$, $0.00\%$] | [$0.00\%$, $0.00\%$] | 1.0000 |
| **EAU PSA Persistence** | 3 / 28 (10.7%) | 3 / 58 (5.2%) | **$-6.1$ pp** | [$-17.8\%$, $+5.6\%$] | [$-15.9\%$, $+3.7\%$] | 0.3052 |
| **Pad-Free Continence 3m** | 4 / 31 (12.9%) | 10 / 51 (19.6%) | **$-4.0$ pp** | [$-20.5\%$, $+12.6\%$] | [$-17.9\%$, $+9.9\%$] | 0.6373 |
| **Pad-Free Continence 12m**| 9 / 18 (50.0%) | 22 / 42 (52.4%) | **$-9.9$ pp** | [$-43.2\%$, $+23.4\%$] | [$-37.8\%$, $+18.0\%$] | 0.5604 |

*Prespecified Equivalence Assessment*:
- Operative time equivalence was prespecified at $\pm 15$ minutes. The 90% confidence interval for adjusted operative time is $[ -2.68, +15.62 ]$ minutes. Because the upper bound ($+15.62$ min) exceeds $+15.0$ minutes, statistical equivalence cannot be declared.
- pT2 margin equivalence was prespecified at $\pm 10$ percentage points. The 90% confidence interval is $[ -8.0\%, +25.1\% ]$. Because the upper bound ($+25.1\%$) exceeds $+10.0\%$, statistical equivalence cannot be declared.

### R6. Bootstrap Uncertainty Estimation for Versius Operative Time $n^*$
Based on 600 case-resampling bootstrap replicates with stochastic multiple imputation:
- **15-minute plateau rule**:
  - Point estimate (full cohort): **Case 111**
  - Bootstrap median: **Case 282** (Mean: 241.8, IQR: [129, 337])
  - **95% Bootstrap Percentile Interval**: **[Case 35, Case 337]**
  - **90% Bootstrap Percentile Interval**: **[Case 40, Case 337]**
- **10-minute plateau rule**:
  - Point estimate (full cohort): **Case 337**
  - Bootstrap median: **Case 337** (Mean: 301.3, IQR: [287, 337])
  - **95% Bootstrap Percentile Interval**: **[Case 96, Case 337]**
  - **90% Bootstrap Percentile Interval**: **[Case 126, Case 337]**

---

## Software, Libraries, and Execution Configuration
- **Statistical Environment**: R 4.6.1 running on Linux x86_64.
- **Core Libraries**:
  - `mgcv` (v1.9-4): Penalized generalized additive models fitted via REML.
  - `mice` (v3.16.0): Multivariate imputation by chained equations (predictive mean matching).
  - `sandwich` (v3.1-0): Heteroskedasticity-consistent covariance estimation (HC3 and HC0).
  - `arrow` (v17.0.0): Fast columnar parquet I/O.
  - `dplyr` (v1.1.4): Data manipulation and pipeline management.
  - `logistf` (v1.26.0): Penalized likelihood for sparse binary outcomes.
  - `jsonlite` (v1.8.9): Standardized JSON serialization.
- **Python Verification Lane**: Python 3.12.12 with `pandas`, `numpy`, `statsmodels`, `scipy`, `duckdb`, `pyarrow`.
- **Hardware Configuration**: Execution bounded to 10 parallel CPU worker processes for bootstrap routines.

---

## Detailed Methodological Implementation

### 1. Cohort Derivation and Sorting Hierarchy
The cohort derivation was implemented strictly according to Section 2 of the locked analysis plan and changelog:
1. Filtered all operations with `country == "Poland"` ($N = 838$).
2. Sorted hierarchically by:
   - `surgery_date` (ascending)
   - `surgery_start_minutes_after_midnight` (ascending, missing values sorted last)
   - `operation_centre_within_day_order` (ascending, missing values sorted last)
   - `full_polish_series_case_sequence` (ascending, deterministic final tie-breaker)
3. Computed experience counters:
   - `platform_n`: Cumulative 1-based index within each platform in the sorted order.
   - `dv_restart_n`: For da Vinci operations with `platform_n >= 20`, defined as `platform_n - 19` ($N = 382$).
   - `day_position`: Order of operations on the surgeon's operating list for that day; `first_case` defined as `day_position == 1`.
   - `data_close`: Latest documented follow-up contact across all 838 Polish records ($39$ days after the final Polish procedure).
   - Window-closed denominators: Required time from surgery to data close $\ge 120$ days for 3-month continence, $\ge 425$ days for 12-month continence, and $\ge 56$ days for EAU PSA persistence.

### 2. Multiple Imputation Strategy
Multiple imputation was executed across all 738 Polish robotic prostatectomies (337 Versius and 401 da Vinci cases) using MICE:
- Missing covariates imputed: `bmi` ($N=67$), `log_weight` ($N=99$), `biopsy_isup_cat` ($N=1$), and `pt_group` ($N=3$).
- Imputation method: Predictive mean matching (`pmm`) across 20 multiply imputed datasets, run for 10 iterations each.
- Predictor matrix: Incorporated all baseline covariates (`age`, `log_psa`, `centre`, `nerve_sparing`, `plnd`, `first_case`), procedural variables (`platform`, `platform_n`, `year`), and auxiliary outcome variables (`log_or_time`, `psm`).
- Mandatory rule: In accordance with the locked plan, imputed outcome values were retained strictly as auxiliaries in the imputation step and were never used as analysis endpoints. All outcome models were fitted only on observed outcomes.

### 3. Standardized Learning Curves and Plateau Determination
- **Model Formulation**:
  - Operative time: $\log(\text{or\_time}) \sim s(\text{case\_n}, k=10) + \mathbf{X}\boldsymbol{\beta}$, fitted by REML.
  - pT2 positive margin: $\text{logit}(\Pr(\text{PSM}=1)) \sim s(\text{case\_n}, k=5) + \mathbf{X}\boldsymbol{\beta}$, fitted by REML restricted to pT2 cases.
- **Separable Terms Standardization**:
  By capitalizing on the additive structure of GAMs without interaction terms, the linear predictor for patient $i$ at case number $n$ decomposes into:
  $$\eta_{i, n} = s(n) + \alpha + \mathbf{x}_i' \boldsymbol{\beta}_{\text{covariates}} = s(n) + g(\mathbf{x}_i)$$
  For operative time, the Duan smearing retransformation over all eligible patients yields:
  $$\widehat{\text{Time}}(n) = \exp(s(n)) \times \left( \frac{1}{N} \sum_{i=1}^N \exp(g(\mathbf{x}_i)) \right) \times \left( \frac{1}{N_{\text{obs}}} \sum_{j=1}^{N_{\text{obs}}} \exp(e_j) \right)$$
  For pT2 margins, the standardized risk is computed as:
  $$\widehat{\Pr}(\text{PSM}=1 \mid n) = \frac{1}{N_{\text{pT2}}} \sum_{i=1}^{N_{\text{pT2}}} \text{logit}^{-1}(s(n) + g(\mathbf{x}_i))$$
- **Asymptote and Plateau Rules**:
  - Asymptote is the arithmetic mean of the standardized curve across the final 50 case numbers ($288\text{--}337$ for Versius; $333\text{--}382$ for da Vinci restart).
  - $n^*$ is the smallest case number $n$ such that for all subsequent cases $m \ge n$, the standardized curve remains within the prespecified band ($|\widehat{Y}(m) - \text{Asymptote}| \le \Delta$).

---

## Ambiguities in the Plan and Methodological Resolutions

1. **Resolution of Same-Day Ties**:
   - *Observation*: Thirteen Polish operations took place on days where operating room start times were unrecorded.
   - *Resolution*: Sorted using `operation_centre_within_day_order`, with `full_polish_series_case_sequence` acting as the final invariant tie-breaker. This preserved a completely deterministic ordering that perfectly reproduced all cohort counts.
2. **Domain of pT2 Margin Curves**:
   - *Observation*: In Section 5 (A1/A2), the plan specifies a penalized spline of Versius case number for pT2 margins, but pT2 cases comprise a non-contiguous subset of 190 operations within the 337-case series.
   - *Resolution*: The spline was modeled on `platform_n` (the surgeon's overall platform experience at the time of surgery) and standardized over the 190 pT2 case-mix profiles across all integer case numbers $1\text{--}337$. This ensures the milestone reporting at cases 1, 25, 50, 100, and 337 aligns directly with the operative time curve.
3. **Zero Events and Quasi-Complete Separation in Calendar 2024 Endpoints**:
   - *Observation 1 (30-Day Readmission)*: In calendar 2024, exactly 0 readmissions occurred across both platforms (0/58 Versius, 0/82 da Vinci restart). Standard logistic regression with covariate adjustment is non-identifiable.
   - *Resolution 1*: The adjusted risk difference is reported as $0.00\%$ ($95\%\text{ CI: } [0.00\%, 0.00\%]$), reflecting identical zero-event rates.
   - *Observation 2 (EAU PSA Persistence)*: In 2024, only 6 events occurred among 86 window-closed cases (3/28 Versius, 3/58 da Vinci restart). Crucially, 0 events occurred among pT2 patients ($0/45$) and 0 events occurred among ISUP grade group 1 patients ($0/28$). Adjusting for the full covariate set produced complete separation.
   - *Resolution 2*: Evaluated using penalized maximum likelihood (Firth's `logistf`) and parsimonious adjustment (controlling for baseline log PSA, age, and centre), yielding an adjusted risk difference of $-6.1$ percentage points ($95\%\text{ CI: } [-17.8\%, +5.6\%]$).

---

## Critical Assessment of Plan Specifications

1. **Instability of the Absolute Minute Plateau Rule ($n^*$)**:
   - The locked plan defines $n^*$ as the first case after which the standardized curve stays within 15 minutes of the average of the last 50 cases.
   - In the Versius series, the operative time curve drops steeply over the first 50 cases, levels off, but continues a slow drift downwards ($-0.23$ min/case) through case 337.
   - While the curve first enters the 15-minute band around the asymptote at case 111, bootstrap resampling (R6) reveals extreme sensitivity: the 95% bootstrap interval for $n^*$ ranges from case 35 to case 337 (median 282).
   - Under the 10-minute band, the curve does not meet the plateau criteria until case 337.
   - *Clinical Implication*: Defining surgical proficiency by an absolute threshold relative to an end-of-series mean is vulnerable to late-series secular trends, changes in team efficiency, or gradual adoption of more complex cases. A slope-based or segmented breakpoint criterion provides a more stable metric.
2. **Equivalence Bounds vs Statistical Power in Concurrent 2024 Cohort**:
   - The analysis plan prespecified equivalence at $\pm 15$ minutes for operative time and $\pm 10$ percentage points for pT2 margin risk difference, requiring the 90% confidence interval to lie entirely within the margin.
   - In 2024 ($N = 204$), the adjusted operative time difference was $+6.47$ minutes with a 90% CI of $[ -2.68, +15.62 ]$ minutes. Because the upper bound barely exceeded 15 minutes (by 37 seconds), equivalence failed on formal criteria, despite no statistically significant difference between platforms ($p = 0.24$).
   - A concurrent sample size of $N=204$ was underpowered to establish equivalence within $\pm 15$ minutes given operative time residual standard deviations of $\sim 40\text{--}45$ minutes.

---

## Verification and Privacy Compliance
- All intermediate patient-level datasets and model objects were restricted to `<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/` with permissions `0600`.
- An automated regex audit across all user-facing files in `repro_results/` confirmed 0 patient identifiers, 0 medical record numbers, and 0 row-level patient dates.
