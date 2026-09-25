# Manuscript6 statistical analysis plan, version 2

Status. Locked on 25 September 2026 before any outcome was examined against case order. Aggregate case mix, completeness and documentation patterns were examined first (scripts 00 to 02). Version 2 integrates the plan discussion by Grok (grok-4.7), AGY (Gemini 3.8 Flash) and muse-free (Muse Spark 1.3). Every change from version 1 and every disagreement not adopted is recorded in `plan/SAP_CHANGELOG.md`. Where this file and the changelog differ, the changelog's version 2 rules govern.

## 1. Design and population

Retrospective analysis of a prospectively maintained registry and a locked, source-verified research database (release v1.1, patient_wide SHA-256 recorded in `evidence/DERIVATION_SUMMARY.json`). One surgeon performed every operation.

- Primary analysis set: 738 consecutive robot-assisted radical prostatectomies (337 Versius, 401 da Vinci) at Salve Medica, Łódź, and the county hospital in Bełchatów, January 2021 to February 2026.
- Descriptive context: 100 laparoscopic prostatectomies by the same surgeon in the same period. They are not part of any learning-curve model.
- No exclusions. Missing outcomes are reported with denominators and never coded as absent events.

## 2. Case order and experience variables

Case order is rebuilt from the verified surgery date, then the operating-room start time, then the hospital's within-day order, then the released sequence as the last tie-breaker (`scripts/01_derive.py`). The released `full_polish_series_case_sequence` is not used because eight registry date corrections were not propagated into it. Thirteen same-day cases without a start time remain unresolved ties and are handled in sensitivity analysis S-ties.

| Variable | Definition |
|---|---|
| `platform_n` | Surgeon's case number on that platform (Versius 1 to 337, da Vinci 1 to 401) |
| `robotic_n` | Surgeon's cumulative robotic case number, both platforms (1 to 738) |
| `team_platform_n` | That hospital's case number on that platform |
| `prior_v`, `prior_dv` | Surgeon's prior cases on each platform |
| `dv_phase` | da Vinci before the first Versius case (13), during the Versius programme (116), after the last Versius case (272) |
| da Vinci restart series | da Vinci cases 20 to 401 (382 cases, 2024 to 2026), numbered `dv_restart_n` 1 to 382. Continuous da Vinci work resumed at case 20 after 240 Versius cases. Cases 1 to 13 (2021, before any Versius case) and 14 to 19 (2022 to 2023, sporadic, gaps of 120 to 175 days) are described separately. |
| 2024 concurrent period | Calendar 2024: Versius cases 235 to 333 (99) and da Vinci restart cases 1 to 105 (105). |
| Data close | The latest dated postoperative contact in the Polish series, 39 days after the last operation. Time-windowed outcomes use window-closed denominators (3-month continence 120 days, 12-month continence 425 days, EAU PSA persistence 56 days). |

## 3. Endpoints

Primary efficiency endpoint. Operative time, registry start-to-end clock minutes (734 of 738 from the clock record). Analysed on the log scale; results back-transformed to minutes.

Co-primary quality endpoint. Positive surgical margin in pT2 disease (190 Versius, 253 da Vinci pT2 cases).

Secondary endpoints.
- Positive surgical margin overall, adjusted for pT group.
- Major complication, Clavien-Dindo grade III or higher, highest documented grade at any time; reoperation; readmission within 30 days; urine leak. Reported with denominators of cases with a stated grade or state. The registry recorded Clavien-Dindo grade I for almost every case through 2023 and recorded "none" or left the grade blank from 2024. Grade I and II rates are therefore not comparable across eras and are not analysed as endpoints.
- Postoperative length of stay.
- Estimated blood loss (the ten implausible repeated 2400 mL values are set to missing). Blood-loss estimation practice may differ between platforms and eras, so it is descriptive and secondary.
- Lymph-node yield among dissections.
- PSA persistence (EAU definition, PSA at or above 0.1 ng/mL on days 28 to 56) on observed denominators.
- Pad-free continence at 3 and 12 months on observed denominators.

Margin-length reporting format differs by hospital and era (Versius-era Salve Medica reports mostly give an inequality). Margin length above 3 mm and margin site are descriptive only.

## 4. Covariates

One prespecified adjustment set, used in every model unless stated: age, BMI, log PSA, biopsy ISUP grade group (1, 2, 3, 4 to 5), log specimen weight, pT group (pT2, pT3a, pT3b or higher; omitted from pT2-restricted models), nerve sparing, pelvic lymph-node dissection, hospital. Operative-time models also include day position (first case of the day versus later). Clinical T stage is not used because its documentation changed between eras (cT1c 15.7% of Versius versus 52.1% of da Vinci cases), which reflects recording practice, not disease.

The pT2 margin models use age, BMI, log PSA, biopsy ISUP grade group, log specimen weight and hospital (no nerve sparing or dissection terms). Missing covariates (BMI 67, specimen weight 99, biopsy ISUP 1, pT 3 among robotic cases) are handled by multiple imputation with chained equations (predictive mean matching, 20 imputations, outcome and experience variables in the imputation model). Estimates are pooled with Rubin's rules. Complete-case analysis is a sensitivity analysis.

## 5. Primary analyses

A1. Versius learning curve. A generalized additive model of the outcome on a penalized regression spline of Versius case number (mgcv, REML, basis dimension 10 for operative time and 5 for pT2 margins), adjusted for the covariate set. Operative time uses a Gaussian model on log minutes with smearing retransformation; pT2 margins use a binomial model. The curve is standardized: at each case number it is the mean prediction over all Versius cases eligible for that outcome. Pointwise 95% intervals come from posterior simulation pooled over imputations.

A2. Plateau. The asymptote is the mean standardized prediction over the series' final 50 cases. The plateau case number n* is the smallest case number after which the curve stays within 15 minutes of the asymptote for every later case (10 minutes in sensitivity), or within 5 percentage points for pT2 margins. A 95% interval for n* comes from 2000 case-resampling bootstrap replicates, each with one stochastic imputation and the rule applied inside the replicate. If the curve never leaves the band, the plateau is reached from case 1. The slope over the final 100 cases is reported as a check.

A3. Transfer of proficiency to the da Vinci restart.
- A3a. Initial phase. Adjusted difference between da Vinci restart cases 1 to 100 and Versius cases 1 to 100 in operative time (minutes, linear model with HC3 errors) and pT2 margin (marginally standardized risk difference). Intercept and slope over the first 50 cases of each series from a model with a series by case-number interaction.
- A3b. Second learning curve. A1 and A2 applied to the restart series with `dv_restart_n`. The smooth term's approximate test against a constant is reported.
- A3c. Platform-specific versus cumulative robotic experience. Over all robotic cases, GAMs with (i) a smooth of cumulative robotic case number plus platform, (ii) platform-specific smooths of platform case number, and (iii) both, compared by AIC and F tests.

A4. Concurrent platform comparison, calendar 2024. Adjusted differences between da Vinci restart cases and Versius cases operated in 2024 (surgeon held constant, same calendar year, both hospitals): operative time (minutes), pT2 margin and all-stage margin adjusted for pT (standardized risk differences), length of stay, 30-day readmission, EAU PSA persistence, and pad-free continence at 3 and 12 months (window-closed, observed denominators). Equivalence is prespecified at ±15 minutes for operative time and ±10 percentage points for pT2 margin, declared only if the 90% interval lies within the margin. Major complications, reoperation and transfusion are reported with exact binomial intervals. Sensitivity adds a natural spline of calendar month (three degrees of freedom). This analysis describes outcomes with the surgeon and year held constant. It does not estimate a causal platform effect.

## 6. Secondary analyses

B1. All-era mature-phase comparison (sensitivity to A4). Versius cases after the Versius plateau against da Vinci restart cases after the restart plateau (all restart cases if the plateau is reached from case 1), same outcomes and models as A4. Era-confounded by design and labelled as such.

B2. Surgeon versus team learning, operative time only. Per platform, a GAM with smooths of the surgeon's platform case number and the hospital's platform case number plus the covariate set. The team term's contribution and the correlation between the two counts are reported. Above a correlation of 0.8 the analysis is shown only as hospital-stratified curves.

B3. Platform alternation cost. Among robotic cases in 2024, whether operative time differs when the surgeon's previous robotic case within seven days was on the other platform, adjusted for the covariate set and each platform's own case number.

B4. Complications (grade III or higher, reoperation, readmission, transfusion), length of stay, blood loss, lymph-node yield among dissections, and margin site and length, descriptive by platform, era and case-order quartile with denominators.

B5. Continence and PSA persistence against case order, exploratory, window-closed observed denominators, with inverse-probability-of-observation weighting (age, hospital, year, pT group, nerve sparing; weights truncated at the 1st and 99th percentiles) as sensitivity.

B6. da Vinci cases 1 to 19 and the 100 laparoscopic cases, descriptive only.

## 7. Sensitivity analyses

S-rcs. Restricted cubic splines (four knots at the 5th, 35th, 65th and 95th percentiles for operative time, three knots for pT2 margins) in place of the penalized spline.
S-10. Plateau margin of 10 minutes.
S-seg. Segmented regression of log operative time on case number, one breakpoint, Davies test, bootstrap interval for the breakpoint.
S-cusum. CUSUM of operative time per series (deviation from the series mean), descriptive.
S-racusum. Risk-adjusted CUSUM (Steiner log-likelihood score) for positive margins at all stages. The expected risk comes from a logistic model over all robotic cases on pT group, final ISUP grade group 3 or higher, log PSA, log specimen weight and hospital; its discrimination and calibration are reported.
S-lccusum. LC-CUSUM for pT2 margins, acceptable 10% and unacceptable 20%, alpha and beta 0.05, threshold set by simulation. Reported whatever the result.
S-ties. 200 random permutations of the 13 unresolved same-day ties; the range of n* is reported.
S-cc. Complete-case analysis.
S-noNSLND. Operative-time models without nerve sparing and lymph-node dissection.
S-NS. pT2 margin models with nerve sparing added.
S-dvall. The da Vinci curve over all 401 da Vinci cases by platform case number.

## 8. Multiplicity and reporting

Two primary endpoints (operative time and pT2 margin) and four primary questions (A1 and A2, A3a, A3b, A4). Every other result is secondary and is reported with intervals, not as a confirmatory test. Reporting follows STROBE and the IDEAL recommendations for new surgical devices, with case numbers at every milestone.

## 9. Assurance

Tier A for A1 to A4 (publication-defining). Objective reproduction by an independent implementation from this plan alone (different model family, different software), with numerical reconciliation of n*, the A3 contrasts and the A4 differences. Methods and implementation review by two model families using code and aggregate output only.

## 10. Privacy

Patient-level data stay in `analysis/restricted/` (mode 0700 and 0600). Every script prints aggregates only and scans its serialized output for date and identifier patterns before writing reader-facing files.
