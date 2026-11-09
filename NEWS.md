# drMultiOutITR — changes

## 0.1.0 (first public release)

- First public source release of the mixed-outcome individualized treatment
  ranking package, with two or more treatment choices.
- Includes a three-treatment worked example with continuous, ordinal, and
  right-censored survival outcomes, plus an optional RHC local-data vignette.
- Provides fitted-model prediction, alternative outcome weights, documentation,
  automated tests, and a GitHub-ready software citation.
- Numerical sensitivity settings retain their established defaults; no
  statistical algorithms were changed for this release.

## 0.0.3.9009 (pre-release development)

- Added step-by-step examples for three treatments with continuous, ordinal,
  and right-censored survival outcomes, plus an optional local-data RHC guide.
- Extended input validation and tests for fitted-model serialization,
  numerical stability, and mixed-outcome treatment ranking.
- Added optional validation utilities for paired survival-distribution and
  factorial simulation comparisons. These are diagnostic tools, not required
  for routine fitting or prediction.
- Clarified the `optimal_ranking` encoding, the interpretation of `theta_hat`,
  and the assumptions needed for observational analyses.

## 0.0.3.9008

- Made simulation validation helpers available to tests in both source and
  installed-package contexts.

## 0.0.3.9006

- Matched the reference implementation's first-grid tie convention when a
  predicted continuous-outcome CDF is constant across all grid points.

## 0.0.3.9005

- Added optional rounding of the forest second-stage training response for
  numerical sensitivity analyses. The default remains unrounded (`NULL`).
  Original AIPW pseudo-outcomes are retained without rounding.

## 0.0.3.9004

- Added options for propensity-score truncation and reproducible seed schedules.
- Added treatment-specific random survival forests as an optional survival
  nuisance learner.

## 0.0.3.9003–0.0.3.9000

- Expanded synthetic tests for AIPW calculations, propensity fitting,
  continuous quantile inversion, survival censoring adjustment, and ranking.
- Added checks comparing the package implementation with the original research
  calculations on matched numerical inputs.

## 0.0.2.9001

- Corrected continuous-outcome grid handling when `grid_length` is supplied
  without explicit thresholds.

## 0.0.2.9000

- Introduced the integrated `fit_dr_itr()` and `predict.dr_itr()` workflow for
  multiple treatments and mixed-type outcomes, with alternative outcome weights
  at prediction time and out-of-bag covariance scaling for forest fits.

## 0.0.1.9001

- Improved documentation links and numerical test stability.
