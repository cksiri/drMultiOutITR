#' Estimate doubly robust conditional favorable-response probabilities
#'
#' Uses the conditional augmented inverse-propensity-weighted (AIPW)
#' pseudo-outcome, followed by a second-stage regression over all training
#' profiles. For binary and ordinal outcomes the first-stage learner models
#' treatment-specific indicator probabilities. Second-stage regression must
#' also be estimated consistently for conditional summary consistency.
#'
#' @param data Observational training data without missing analysis columns.
#' @param treatment Treatment column name.
#' @param outcome Response column name.
#' @param covariates Baseline covariate column names.
#' @param propensity A `dr_propensity` object from [fit_dr_propensity()], or
#'   an n-by-K propensity matrix with named treatment columns.
#' @param newdata Data frame for prediction, defaulting to training profiles.
#' @param favorable Favorable observed value for the binary response.
#' @param outcome_learner Binary indicator learner: `"ranger"` or `"glm"`.
#' @param second_stage Pseudo-outcome regression: `"ranger"` or `"lm"`.
#' @param fold_id Optional matching treatment-stratified nuisance fold vector.
#'   Automatically inherited from a `dr_propensity` object. Both the
#'   propensity and outcome nuisance predictions should be out of fold.
#' @param seed Seed for the stochastic learner.
#' @param num.trees Forest size when using `ranger`.
#' @param min.node.size Forest minimum terminal-node size.
#' @param ids Optional unique identifiers for rows in `newdata`.
#' @param return_models Return final regressions and pseudo-outcomes for audit.
#' @param seed_schedule `"package"` or `"research"`; the latter follows
#'   the component/arm forest seed offsets of files 02 and 03.
#' @return `dr_summary_estimate` with `theta_array` and `theta_long`. These
#'   are estimated conditional counterfactual summaries, not individually
#'   observed causal effects.
#' @param second_stage_round_digits `NULL` (default, unchanged original AIPW
#'   training responses) or an integer 0--15. When using a `ranger` second
#'   stage, only the response passed to that forest is rounded to this many
#'   decimal places. The exact AIPW pseudo-outcomes remain stored unmodified.
#'   Intended as an optional reproducibility diagnostic, not a DR assumption;
#'   use sensitivity analyses before adopting for scientific results.
#'   Not supported with `second_stage = "lm"`.
#' @export
estimate_dr_binary <- function(data, treatment, outcome, covariates, propensity,
                               newdata = data, favorable = 1,
                               outcome_learner = c('ranger', 'glm'),
                               second_stage = c('ranger', 'lm'),
                               fold_id = NULL, seed = 2026L,
                               num.trees = 500L, min.node.size = 10L,
                               ids = NULL, return_models = FALSE,
                               seed_schedule = c('package', 'research'),
                               second_stage_round_digits = NULL) {
  outcome_learner <- match.arg(outcome_learner)
  second_stage <- match.arg(second_stage)
  .dr_check_data(data, treatment, covariates, outcome)
  h <- as.numeric(as.character(data[[outcome]]) == as.character(favorable))
  .dr_components(data, treatment, covariates, propensity, newdata,
                 H = matrix(h, ncol = 1L), labels = 'Pr_favorable',
                 binary = TRUE, outcome_learner = outcome_learner,
                 second_stage = second_stage, fold_id = fold_id, seed = seed,
                 num.trees = num.trees, min.node.size = min.node.size,
                 monotonic = 'none', ids = ids, return_models = return_models,
                 seed_schedule = seed_schedule,
                 second_stage_round_digits = second_stage_round_digits)
}

#' Estimate doubly robust ordinal cumulative probabilities
#'
#' For higher-is-better outcomes estimates P(Y > c | X,A=k) at each ordered
#' cutpoint c; for lower-is-better outcomes estimates P(Y <= c | X,A=k).
#' An `ordered` factor or numeric response is required to ensure that the
#' category order is clinically specified rather than alphabetical.
#'
#' @inheritParams estimate_dr_binary
#' @param direction `"higher"` or `"lower"`, describing favorable categories.
#' @param levels Optional ordered category values for numeric/character data.
#'   For an ordered factor, its stored level order is respected automatically.
#' @export
estimate_dr_ordinal <- function(data, treatment, outcome, covariates, propensity,
                                newdata = data, direction = c('higher', 'lower'),
                                levels = NULL,
                                outcome_learner = c('ranger', 'glm'),
                                second_stage = c('ranger', 'lm'),
                                fold_id = NULL, seed = 2026L,
                                num.trees = 500L, min.node.size = 10L,
                                ids = NULL, return_models = FALSE,
                               seed_schedule = c('package', 'research'),
                               second_stage_round_digits = NULL) {
  direction <- match.arg(direction)
  outcome_learner <- match.arg(outcome_learner)
  second_stage <- match.arg(second_stage)
  .dr_check_data(data, treatment, covariates, outcome)
  y <- data[[outcome]]
  if (!is.null(levels)) {
    if (anyDuplicated(levels) || length(levels) < 2L ||
        !all(as.character(y) %in% as.character(levels))) {
      stop('levels must give all ordered category values in increasing order.')
    }
    category <- match(as.character(y), as.character(levels))
    category_labels <- as.character(levels)
  } else if (is.ordered(y)) {
    category <- as.integer(y)
    category_labels <- base::levels(y)
  } else if (is.numeric(y)) {
    category_labels <- as.character(sort(unique(y)))
    category <- match(y, sort(unique(y)))
  } else {
    stop('Ordinal outcomes must be numeric, ordered factors, or have explicit levels.')
  }
  if (length(category_labels) < 2L) stop('Ordinal response needs at least two categories.')
  cut <- seq_len(length(category_labels) - 1L)
  h <- vapply(cut, function(c) {
    as.numeric(if (direction == 'higher') category > c else category <= c)
  }, numeric(length(category)))
  if (is.null(dim(h))) h <- matrix(h, ncol = length(cut))
  labels <- paste0(if (direction == 'higher') 'P_gt_' else 'P_le_',
                   category_labels[cut])
  out <- .dr_components(data, treatment, covariates, propensity, newdata,
                        H = h, labels = labels, binary = TRUE,
                        outcome_learner = outcome_learner, second_stage = second_stage,
                        fold_id = fold_id, seed = seed, num.trees = num.trees,
                        min.node.size = min.node.size,
                        monotonic = if (direction == 'higher') 'decreasing' else 'increasing',
                        ids = ids, return_models = return_models,
                        seed_schedule = seed_schedule,
                 second_stage_round_digits = second_stage_round_digits)
  out$thresholds <- category_labels[cut]
  out$direction <- direction
  out
}

.dr_quantile_inverse <- function(cdf, grid, probs, interpolation = 'grid') {
  if (interpolation == 'grid') {
    return(vapply(probs, function(p) {
      idx <- which(cdf >= p)
      if (length(idx)) grid[idx[1L]] else grid[length(grid)]
    }, numeric(1)))
  }
  keep <- !duplicated(cdf)
  cdf_u <- cdf[keep]
  y_u <- grid[keep]
  # A constant fitted CDF does not identify a unique quantile. Use the
  # first grid point on ties for compatibility with the study implementation.
  if (length(cdf_u) < 2L) return(rep(grid[1L], length(probs)))
  as.numeric(stats::approx(cdf_u, y_u, xout = probs, rule = 2,
                           ties = 'ordered')$y)
}

#' Estimate doubly robust conditional continuous-outcome quantiles
#'
#' Fits treatment-specific conditional CDF indicators on an outcome grid,
#' constructs conditional AIPW pseudo-outcomes, regresses these on X, clips
#' projected CDFs to the unit interval, enforces monotonicity and inverts the CDF.
#' Values outside the grid cannot be recovered; grid endpoints may produce
#' boundary quantiles. No outcome-level policy value is estimated.
#'
#' @inheritParams estimate_dr_binary
#' @param probs Quantile levels strictly between 0 and 1.
#' @param grid Optional increasing set of outcome thresholds. By default a
#'   sequence between training-response 2nd and 98th percentiles.
#' @param grid_length Threshold count when `grid` is missing.
#' @param interpolation `"grid"` or `"interpolate"` quantile inversion.
#' @param return_cdf Return the clipped and monotone conditional CDF array.
#' @export
estimate_dr_continuous <- function(data, treatment, outcome, covariates, propensity,
                                   newdata = data, probs = c(.25, .5, .75),
                                   grid = NULL, grid_length = 41L,
                                   interpolation = c('grid', 'interpolate'),
                                   outcome_learner = c('ranger', 'glm'),
                                   second_stage = c('ranger', 'lm'),
                                   fold_id = NULL, seed = 2026L,
                                   num.trees = 500L, min.node.size = 10L,
                                   ids = NULL, return_cdf = TRUE,
                                   return_models = FALSE,
                                   seed_schedule = c('package', 'research'),
                               second_stage_round_digits = NULL) {
  interpolation <- match.arg(interpolation)
  outcome_learner <- match.arg(outcome_learner)
  second_stage <- match.arg(second_stage)
  .dr_check_data(data, treatment, covariates, outcome)
  y <- data[[outcome]]
  if (!is.numeric(y) || any(!is.finite(y))) stop('Continuous outcome must be finite numeric.')
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1) || is.unsorted(probs, strictly = TRUE)) {
    stop('probs must be strictly increasing probabilities between 0 and 1.')
  }
  if (is.null(grid)) {
    if (length(grid_length) != 1L || grid_length < 3L) stop('grid_length must be >=3.')
    endpoints <- as.numeric(stats::quantile(y, c(.02, .98)))
    if (endpoints[1L] == endpoints[2L]) endpoints <- range(y)
    if (endpoints[1L] == endpoints[2L]) stop('A nondegenerate outcome grid is required.')
    grid <- seq(endpoints[1L], endpoints[2L], length.out = grid_length)
  }
  grid <- sort(unique(grid))
  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid))) {
    stop('grid must contain at least two distinct finite numeric thresholds.')
  }
  H <- vapply(grid, function(v) as.numeric(y <= v), numeric(length(y)))
  if (is.null(dim(H))) H <- matrix(H, ncol = length(grid))
  out <- .dr_components(data, treatment, covariates, propensity, newdata,
                        H, labels = paste0('F_', seq_along(grid)),
                        binary = TRUE, outcome_learner = outcome_learner,
                        second_stage = second_stage, fold_id = fold_id,
                        seed = seed, num.trees = num.trees,
                        min.node.size = min.node.size,
                        monotonic = 'increasing', ids = ids,
                        return_models = return_models,
                        seed_schedule = seed_schedule,
                 second_stage_round_digits = second_stage_round_digits)
  cdf <- out$theta_array
  res <- array(NA_real_, dim = c(dim(cdf)[1L], dim(cdf)[2L], length(probs)),
               dimnames = list(new_id = dimnames(cdf)[[1L]],
                               treatment = dimnames(cdf)[[2L]],
                               summary = paste0('q', format(probs, trim = TRUE))))
  for (i in seq_len(dim(cdf)[1L])) for (k in seq_len(dim(cdf)[2L])) {
    res[i, k, ] <- .dr_quantile_inverse(as.numeric(cdf[i, k, ]), grid,
                                       probs, interpolation)
  }
  ids_out <- dimnames(res)[[1L]]
  tl <- do.call(rbind, lapply(seq_len(dim(res)[1L]), function(i) {
    do.call(rbind, lapply(seq_along(probs), function(s) {
      data.frame(new_id = ids_out[i], summary_id = dimnames(res)[[3L]][s],
                 treatment = dimnames(res)[[2L]],
                 theta_hat = as.numeric(res[i, , s]), stringsAsFactors = FALSE)
    }))
  }))
  rownames(tl) <- NULL
  out$theta_array <- res
  out$theta_long <- tl
  out$summary_ids <- dimnames(res)[[3L]]
  out$probs <- probs
  out$grid <- grid
  out$interpolation <- interpolation
  if (return_cdf) out$cdf_array <- cdf
  out
}
