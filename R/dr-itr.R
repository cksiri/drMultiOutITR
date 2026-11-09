# Fit conditional outcome summaries and retain second-stage models for
# treatment ranking on new baseline profiles.

.dr_itr_specs <- function(outcomes, data) {
  if (!is.list(outcomes) || !length(outcomes) || is.null(names(outcomes)) ||
      anyNA(names(outcomes)) || any(!nzchar(names(outcomes))) ||
      anyDuplicated(names(outcomes))) {
    stop('outcomes must be a nonempty, uniquely named list of specifications.',
         call. = FALSE)
  }
  types <- c('binary', 'ordinal', 'continuous', 'survival')
  result <- outcomes
  for (id in names(result)) {
    s <- result[[id]]
    if (!is.list(s) || !is.character(s$type) || length(s$type) != 1L ||
        is.na(s$type) || !s$type %in% types) {
      stop('Outcome ', id, ' requires type binary, ordinal, continuous or survival.',
           call. = FALSE)
    }
    if (s$type == 'survival') {
      if (!is.character(s$time) || length(s$time) != 1L ||
          !is.character(s$status) || length(s$status) != 1L ||
          anyNA(c(s$time, s$status)) ||
          !all(c(s$time, s$status) %in% names(data))) {
        stop('Survival outcome ', id, ' requires time and status column names.',
             call. = FALSE)
      }
      if (is.null(s$times)) stop('Survival outcome ', id, ' requires times.')
      if (is.null(s$censoring)) s$censoring <- 'pooled_km'
      if (!s$censoring %in% c('pooled_km', 'cox')) {
        stop('Invalid censoring choice for ', id)
      }
      if (is.null(s$min_G)) s$min_G <- 0.05
    } else {
      if (!is.character(s$column) || length(s$column) != 1L ||
          is.na(s$column) || !s$column %in% names(data)) {
        stop('Outcome ', id, ' requires an existing column name.', call. = FALSE)
      }
      if (s$type == 'binary' && is.null(s$favorable)) s$favorable <- 1
      if (s$type == 'ordinal') {
        if (is.null(s$direction)) s$direction <- 'higher'
        if (!s$direction %in% c('higher', 'lower')) {
          stop('Ordinal direction must be higher or lower for ', id)
        }
      }
      if (s$type == 'continuous') {
        if (is.null(s$direction)) s$direction <- 'higher'
        if (!s$direction %in% c('higher', 'lower')) {
          stop('Continuous direction must be higher or lower for ', id)
        }
        if (is.null(s$probs)) s$probs <- c(.25, .5, .75)
        if (is.null(s$grid_length)) s$grid_length <- 41L
        if (is.null(s$interpolation)) s$interpolation <- 'grid'
      }
    }
    if (is.null(s$rho)) s$rho <- 0.5
    if (is.null(s$summary_weight)) s$summary_weight <- 1
    if (is.null(s$summary_scale)) s$summary_scale <- 1
    for (nm in c('rho', 'summary_scale', 'summary_weight')) {
      val <- s[[nm]]
      if (!is.numeric(val) || length(val) != 1L || !is.finite(val) ||
          (nm != 'summary_weight' && val <= 0) ||
          (nm == 'summary_weight' && val < 0)) {
        stop('Outcome ', id, ' has invalid ', nm, '.', call. = FALSE)
      }
    }
    result[[id]] <- s
  }
  result
}

.dr_itr_weights <- function(weights, outcome_ids) {
  if (is.null(weights)) {
    return(stats::setNames(rep(1 / length(outcome_ids), length(outcome_ids)),
                           outcome_ids))
  }
  if (!is.numeric(weights) || length(weights) != length(outcome_ids) ||
      any(!is.finite(weights)) || any(weights < 0) || all(weights == 0)) {
    stop('weights must be finite, nonnegative and have one entry per outcome; ',
         'at least one must be positive.', call. = FALSE)
  }
  if (is.null(names(weights))) {
    names(weights) <- outcome_ids
  } else {
    if (anyDuplicated(names(weights)) || !setequal(names(weights), outcome_ids)) {
      stop('Named weights must use exactly the outcome specification names.',
           call. = FALSE)
    }
    weights <- weights[outcome_ids]
  }
  weights
}

.dr_itr_covariate_schema <- function(data, covariates) {
  lapply(data[, covariates, drop = FALSE], function(x) {
    if (is.factor(x)) list(kind = 'factor', levels = levels(x), ordered = is.ordered(x))
    else if (is.character(x)) list(kind = 'factor', levels = unique(x), ordered = FALSE)
    else if (is.numeric(x)) list(kind = 'numeric')
    else if (is.logical(x)) list(kind = 'logical')
    else stop('Covariates must be numeric, logical, character or factor.')
  })
}

.dr_itr_cast <- function(data, schema) {
  if (!is.data.frame(data) || nrow(data) < 1L ||
      any(!names(schema) %in% names(data))) {
    stop('Data must contain every fitted baseline covariate.', call. = FALSE)
  }
  if (anyNA(data[, names(schema), drop = FALSE])) {
    stop('Missing baseline covariate values; preprocess before fitting/prediction.',
         call. = FALSE)
  }
  for (v in names(schema)) {
    info <- schema[[v]]
    x <- data[[v]]
    if (info$kind == 'factor') {
      unexpected <- setdiff(unique(as.character(x)), info$levels)
      if (length(unexpected)) {
        stop('Unknown level in covariate ', v, ': ', paste(unexpected, collapse = ', '),
             call. = FALSE)
      }
      data[[v]] <- factor(as.character(x), levels = info$levels,
                          ordered = isTRUE(info$ordered))
    } else if (info$kind == 'numeric') {
      if (!is.numeric(x) || any(!is.finite(x))) {
        stop('Covariate ', v, ' must be finite numeric.', call. = FALSE)
      }
    } else if (info$kind == 'logical') {
      if (!is.logical(x)) stop('Covariate ', v, ' must be logical.', call. = FALSE)
    }
  }
  data
}

.dr_itr_long <- function(array, ids) {
  lev <- dimnames(array)[[2L]]
  labs <- dimnames(array)[[3L]]
  rows <- vector('list', length(ids) * length(labs))
  ix <- 0L
  for (i in seq_along(ids)) for (s in seq_along(labs)) {
    ix <- ix + 1L
    rows[[ix]] <- data.frame(new_id = as.character(ids[i]),
                              summary_id = labs[s], treatment = lev,
                              theta_hat = as.numeric(array[i, , s]),
                              stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  tab
}

.dr_itr_model_predictions <- function(est, x, covariates, ids,
                                       use_oob = FALSE) {
  if (!inherits(est, 'dr_summary_estimate') || is.null(est$models)) {
    stop('Internal error: missing fitted second-stage models.', call. = FALSE)
  }
  K <- length(est$treatment_levels)
  # Continuous estimators store one model for each CDF grid threshold,
  # whereas their public theta_array holds the *inverted* quantiles.
  is_cdf <- !is.null(est$grid) && !is.null(est$probs)
  J <- if (is_cdf) length(est$grid) else length(est$summary_ids)
  raw_labels <- if (is_cdf) paste0('F_', seq_len(J)) else est$summary_ids
  raw <- array(NA_real_, c(nrow(x), K, J),
               dimnames = list(new_id = as.character(ids),
                               treatment = est$treatment_levels,
                               summary = raw_labels))
  for (s in seq_len(J)) for (k in seq_len(K)) {
    fit <- est$models[[(s - 1L) * K + k]]$final_fit
    if (is.null(fit)) stop('Internal error: second-stage model absent.')
    if (use_oob && fit$kind == 'ranger') {
      # For ranger regression fits this is the training-row OOB prediction.
      y <- fit$model$predictions
      if (!is.numeric(y) || length(y) != nrow(x) || any(!is.finite(y))) {
        stop('Non-finite or missing ranger OOB predictions. Increase trees or ',
             'check training size.', call. = FALSE)
      }
      raw[, k, s] <- as.numeric(y)
    } else if (use_oob && fit$kind != 'constant') {
      stop('OOB covariance scaling requires a ranger second-stage learner; ',
           'choose covariance_source="in_sample" for lm.', call. = FALSE)
    } else {
      raw[, k, s] <- .dr_predict_model(fit, x, covariates)
    }
  }
  raw <- .dr_clip(raw)
  orientation <- if (is_cdf) 'increasing' else if (!is.null(est$direction)) {
    if (est$direction == 'higher') 'decreasing' else 'increasing'
  } else if (!is.null(est$times)) 'decreasing' else 'none'
  if (orientation != 'none' && J > 1L) {
    for (i in seq_len(dim(raw)[1L])) for (k in seq_len(K)) {
      raw[i, k, ] <- if (orientation == 'increasing') {
        cummax(raw[i, k, ])
      } else cummin(raw[i, k, ])
    }
  }
  result <- raw
  if (is_cdf) {
    result <- array(NA_real_, c(nrow(x), K, length(est$probs)),
                    dimnames = list(new_id = as.character(ids),
                                    treatment = est$treatment_levels,
                                    summary = est$summary_ids))
    for (i in seq_len(nrow(x))) for (k in seq_len(K)) {
      result[i, k, ] <- .dr_quantile_inverse(
        as.numeric(raw[i, k, ]), est$grid, est$probs, est$interpolation)
    }
  }
  out <- list(theta_array = result, theta_long = .dr_itr_long(result, ids),
              treatment_levels = est$treatment_levels,
              summary_ids = dimnames(result)[[3L]])
  if (is_cdf) out$cdf_array <- raw
  class(out) <- 'dr_summary_estimate'
  out
}

.dr_itr_table <- function(ests, specs, x, covariates, ids, weights,
                          use_oob = FALSE) {
  tabs <- lapply(names(specs), function(id) {
    spec <- specs[[id]]
    est <- .dr_itr_model_predictions(ests[[id]], x, covariates, ids,
                                      use_oob = use_oob)
    as_dr_summary_table(est, outcome_id = id,
      larger_is_better = spec$type != 'continuous' || spec$direction == 'higher',
      outcome_weight = as.numeric(weights[id]), rho = spec$rho,
      summary_weight = spec$summary_weight,
      summary_scale = spec$summary_scale)
  })
  do.call(rbind, tabs)
}

#' Fit a doubly robust multivariate individualized treatment ranking rule
#'
#' End-to-end interface for observational data with two or more treatment arms
#' and mixed binary, ordinal, continuous and right-censored outcomes. The fit
#' combines treatment-propensity nuisance models, outcome-specific conditional
#' AIPW summary estimates, and a candidate-specific covariance-scaled ranking
#' rule. All outcome learners use the **same treatment-stratified nuisance
#' folds**. The second-stage fits are retained so that [predict()] can apply
#' the rule to new baseline profiles without refitting nuisance functions.
#'
#' For forest-based second-stage regression, the default covariance matrices
#' are fitted on out-of-bag (OOB) predictions for the training profiles.
#' `covariance_source="in_sample"` is available for comparison or for an `lm`
#' second stage; this option does not provide OOB separation. An *outer* holdout
#' split is required to evaluate recommendations independently of fitting.
#'
#' @param data Complete observational training data. Missing baseline predictors
#'   or outcome values must be handled before calling this function.
#' @param treatment Column name for observed treatment assignment.
#' @param covariates Character vector of pretreatment adjustment columns.
#' @param outcomes Named list of outcome specifications. Each entry has a `type`:
#'   `"binary"` with `column` and optional `favorable` (default `1`);
#'   `"ordinal"` with `column`, optional `direction` (`"higher"`/`"lower"`),
#'   and optional category `levels`; `"continuous"` with `column`, optional
#'   `direction`, `probs`, `grid`, `grid_length`, and `interpolation`;
#'   or `"survival"` with `time`, `status`, required evaluation `times`, and
#'   optional `censoring` (`"pooled_km"` or `"cox"`) and `min_G`. Per-outcome
#'   optional fields include `rho` (default 0.5), `summary_weight` (default 1)
#'   and `summary_scale` (default 1). Ordinal probabilities are oriented so
#'   that higher is always favorable; lower-is-better continuous quantiles
#'   retain their original measurement scale.
#' @param weights Optional nonnegative vector of outcome weights, named using
#'   exactly the names in `outcomes`, or unnamed in outcome-list order. Equal
#'   weights are used by default.
#' @param propensity_method `"ranger"` (default) or `"multinomial"`.
#' @param nuisance_folds Treatment-stratified cross-fitting folds; 1 disables
#'   nuisance cross-fitting. The fold assignment is shared across outcomes.
#' @param propensity_eps Minimum treatment-assignment probability; fixed
#'   truncation can introduce bias and compromise propensity-only robustness.
#' @param outcome_learner Outcome-nuisance learner: `"ranger"`, `"glm"`
#'   for non-survival binary indicators; `"ranger"` or `"lm"` for survival.
#' @param second_stage `"ranger"` (default) or `"lm"` for regression of
#'   AIPW pseudo-outcomes on baseline covariates.
#' @param covariance_source `"oob"` (default) or `"in_sample"`. OOB requires
#'   `second_stage="ranger"`; in-sample scaling is less separated from fitting.
#' @param seed Random seed used for folds and learners.
#' @param num.trees Number of forest trees used for outcome and second stages.
#' @param propensity_trees Number of propensity-forest trees.
#' @param min.node.size Minimum forest terminal node size.
#' @param propensity_min_node_size Optional propensity forest node size;
#'   defaults to `min.node.size`. File 30 used 50 for propensities.
#' @param propensity_seed Optional propensity forest seed; file 30 used
#'   master seed plus 9000.
#' @param propensity_truncation Probability adjustment: `"simplex"` or
#'   original file-01 `"floor_renorm"`.
#' @param seed_schedule `"package"` keeps historical offsets, while
#'   `"research"` uses the seed schedules of files 02, 03 and 29.
#' @param survival_outcome_learner `"transformed"` fits IPCW-response
#'   regression; `"rsf"` uses treatment-specific randomForestSRC event
#'   survival forests for first-stage nuisances as in manuscript file 30.
#' @param survival_num.trees,survival_min_node_size Optional separate survival
#'   forest size and final-stage node size, default to the common settings.
#'   RSF first-stage nodesize matches `survival_min_node_size`.
#' @param second_stage_round_digits `NULL` (default) keeps the original AIPW
#'   response for every second-stage forest. An integer from 0 to 15 rounds
#'   the forest training response while retaining exact pseudo-outcomes.
#'   This is optional numerical stabilization, not a change to the DR formula.
#' @param survival_second_stage_round_digits Optional survival-specific override;
#'   defaults to `second_stage_round_digits`. Use `12L` to stabilize survival
#'   forest responses without changing continuous/ordinal/binary learners.
#'   These settings require `second_stage = "ranger"`.
#' @param ridge Relative ridge regularization for candidate covariance matrices.
#' @param eigen_floor Minimum relative covariance eigenvalue.
#' @param magnitude_scale Rank engine magnitude scaling choice.
#' @param ties.method Method for tied outcome ranks passed to [base::rank()].
#' @return An object of class `dr_itr` containing fitted conditional models,
#'   a frozen `dr_rank_rule`, outcome specifications, and training diagnostics.
#'   The recommended treatment is a rank-aggregation output under causal
#'   identification assumptions, *not* an externally validated causal policy.
#' @export
#' @examples
#' # See inst/examples/two_treatment_example.R for a complete synthetic worked example.
fit_dr_itr <- function(data, treatment, covariates, outcomes, weights = NULL,
                       propensity_method = c('ranger', 'multinomial'),
                       nuisance_folds = 5L, propensity_eps = .05,
                       outcome_learner = 'ranger',
                       second_stage = c('ranger', 'lm'),
                       covariance_source = c('oob', 'in_sample'),
                       seed = 2026L, num.trees = 300L,
                       propensity_trees = 300L, min.node.size = 10L,
                       propensity_min_node_size = min.node.size,
                       propensity_seed = NULL,
                       propensity_truncation = c('simplex', 'floor_renorm'),
                       seed_schedule = c('package', 'research'),
                       survival_outcome_learner = c('transformed', 'rsf'),
                       survival_num.trees = num.trees,
                       survival_min_node_size = min.node.size,
                       second_stage_round_digits = NULL,
                       survival_second_stage_round_digits = second_stage_round_digits,
                       ridge = 1e-6, eigen_floor = 1e-8,
                       magnitude_scale = c('none', 'summary_scale'),
                       ties.method = 'average') {
  propensity_method <- match.arg(propensity_method)
  second_stage <- match.arg(second_stage)
  covariance_source <- match.arg(covariance_source)
  magnitude_scale <- match.arg(magnitude_scale)
  propensity_truncation <- match.arg(propensity_truncation)
  seed_schedule <- match.arg(seed_schedule)
  survival_outcome_learner <- match.arg(survival_outcome_learner)
  second_stage_round_digits <- .dr_validate_forest_digits(
    second_stage_round_digits, second_stage)
  survival_second_stage_round_digits <- .dr_validate_forest_digits(
    survival_second_stage_round_digits, second_stage)
  if (is.null(propensity_seed)) {
    propensity_seed <- if (seed_schedule == 'research') seed + 9000L else seed
  }
  if (covariance_source == 'oob' && second_stage != 'ranger') {
    stop('covariance_source="oob" requires second_stage="ranger".', call. = FALSE)
  }
  specs <- .dr_itr_specs(outcomes, data)
  outcome_ids <- names(specs)
  weights <- .dr_itr_weights(weights, outcome_ids)
  outcome_columns <- unique(unlist(lapply(specs, function(s) {
    if (s$type == 'survival') c(s$time, s$status) else s$column
  }), use.names = FALSE))
  levels <- .dr_check_data(data, treatment, covariates, outcome_columns)
  if (propensity_eps >= 1 / length(levels)) {
    stop('propensity_eps must be less than 1/K.', call. = FALSE)
  }
  # R forests handle factors directly; preserve levels for later validation.
  schema <- .dr_itr_covariate_schema(data, covariates)
  data <- .dr_itr_cast(data, schema)
  if (!is.numeric(num.trees) || length(num.trees) != 1L ||
      !is.finite(num.trees) || num.trees < 1 || num.trees != as.integer(num.trees)) {
    stop('num.trees must be a positive integer.', call. = FALSE)
  }
  if (!is.numeric(propensity_trees) || length(propensity_trees) != 1L ||
      !is.finite(propensity_trees) || propensity_trees < 1 ||
      propensity_trees != as.integer(propensity_trees)) {
    stop('propensity_trees must be a positive integer.', call. = FALSE)
  }
  node_settings <- list(min.node.size = min.node.size,
                        propensity_min_node_size = propensity_min_node_size,
                        survival_min_node_size = survival_min_node_size)
  for (nm in names(node_settings)) {
    x <- node_settings[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x != as.integer(x)) {
      stop(nm, ' must be a positive integer.', call. = FALSE)
    }
  }
  if (!is.numeric(survival_num.trees) || length(survival_num.trees) != 1L ||
      !is.finite(survival_num.trees) || survival_num.trees < 1 ||
      survival_num.trees != as.integer(survival_num.trees)) {
    stop('survival_num.trees must be a positive integer.', call. = FALSE)
  }
  if (!is.numeric(propensity_seed) || length(propensity_seed) != 1L ||
      !is.finite(propensity_seed)) stop('propensity_seed must be finite scalar.')
  ps <- fit_dr_propensity(data, treatment, covariates,
                          method = propensity_method, folds = nuisance_folds,
                          eps = propensity_eps, num.trees = propensity_trees,
                          min.node.size = propensity_min_node_size,
                          seed = propensity_seed,
                          truncation = propensity_truncation,
                          seed_schedule = seed_schedule)
  estimates <- vector('list', length(specs))
  names(estimates) <- outcome_ids
  train_ids <- as.character(seq_len(nrow(data)))
  for (i in seq_along(specs)) {
    s <- specs[[i]]
    research_offset <- switch(s$type,
                              ordinal = 100L, continuous = 200L,
                              survival = 300L, binary = 400L + i)
    outcome_seed <- if (seed_schedule == 'research') {
      seed + research_offset
    } else seed + 100000L * i
    common <- list(data = data, treatment = treatment,
                   covariates = covariates, propensity = ps,
                   newdata = data, outcome_learner = outcome_learner,
                   second_stage = second_stage, seed = outcome_seed,
                   num.trees = num.trees, min.node.size = min.node.size,
                   ids = train_ids, return_models = TRUE,
                   seed_schedule = seed_schedule,
                   second_stage_round_digits = if (s$type == 'survival') {
                     survival_second_stage_round_digits
                   } else second_stage_round_digits)
    if (s$type == 'binary') {
      args <- c(common, list(outcome = s$column, favorable = s$favorable))
      estimates[[i]] <- do.call(estimate_dr_binary, args)
    } else if (s$type == 'ordinal') {
      args <- c(common, list(outcome = s$column,
                             direction = s$direction, levels = s$levels))
      estimates[[i]] <- do.call(estimate_dr_ordinal, args)
    } else if (s$type == 'continuous') {
      # Use exact list indexing: $grid can partially match grid_length when
      # no explicit threshold vector was supplied. Passing the count as a
      # singleton grid causes estimate_dr_continuous() to reject it.
      args <- c(common, list(outcome = s$column, probs = s$probs,
                             grid = s[["grid", exact = TRUE]],
                             grid_length = s$grid_length,
                             interpolation = s$interpolation,
                             return_cdf = FALSE))
      estimates[[i]] <- do.call(estimate_dr_continuous, args)
    } else {
      if (survival_outcome_learner == 'transformed' &&
          !outcome_learner %in% c('ranger', 'lm')) {
        stop('Transformed survival outcomes require outcome_learner="ranger" or "lm".')
      }
      args <- common
      args$outcome_learner <- if (survival_outcome_learner == 'rsf') 'rsf' else outcome_learner
      args$num.trees <- survival_num.trees
      args$min.node.size <- survival_min_node_size
      args <- c(args, list(time = s$time, status = s$status,
                           times = s$times, censoring = s$censoring,
                           min_G = s$min_G,
                           rsf_num.trees = survival_num.trees,
                           rsf_nodesize = survival_min_node_size))
      estimates[[i]] <- do.call(estimate_dr_survival, args)
    }
  }
  train_table <- .dr_itr_table(estimates, specs, data, covariates,
                               train_ids, weights,
                               use_oob = identical(covariance_source, 'oob'))
  rank_fit <- fit_rank_rule(train_table, magnitude_scale = magnitude_scale,
                             ridge = ridge, eigen_floor = eigen_floor,
                             ties.method = ties.method)
  # Baseline data are intentionally not stored as a separate copy.
  structure(list(treatment = treatment, covariates = covariates,
                 treatments = levels, outcomes = specs, weights = weights,
                 propensity = ps, estimates = estimates,
                 rank_rule = rank_fit, schema = schema,
                 covariance_source = covariance_source,
                 n_train = nrow(data), nuisance_folds = nuisance_folds,
                 num.trees = num.trees, seed = seed,
                 seed_schedule = seed_schedule,
                 propensity_seed = propensity_seed,
                 propensity_truncation = propensity_truncation,
                 propensity_min_node_size = propensity_min_node_size,
                 survival_outcome_learner = survival_outcome_learner,
                 survival_num.trees = survival_num.trees,
                 survival_min_node_size = survival_min_node_size,
                 second_stage_round_digits = second_stage_round_digits,
                 survival_second_stage_round_digits = survival_second_stage_round_digits),
            class = 'dr_itr')
}

#' Predict treatments for new patient profiles from a fitted DR rule
#'
#' Applies the stored conditional DR second-stage regressions to each new
#' baseline profile, constructs treatment-specific outcome summaries, and
#' applies the **frozen** covariance-scaled ranking rule. No nuisance models
#' or covariance matrices are refitted. Input should include baseline
#' covariates only; neither assigned treatment nor outcomes are required.
#'
#' @param object A fitted `dr_itr` returned by [fit_dr_itr()].
#' @param newdata Data frame of new baseline covariate profiles.
#' @param weights Optional alternative outcome weights; named or in the
#'   original outcome order. Covariance matrices are *not* recomputed.
#' @param ids Optional distinct identifiers for each prediction profile.
#' @param details Whether to return candidate-level ranking scores.
#' @param return_summaries Whether to return estimated conditional outcome
#'   summaries in standard long format alongside the recommendations.
#' @param ... Reserved; additional arguments are not currently supported.
#' @return By default a `data.frame` of treatment selections. If either
#'   `details` or `return_summaries` is `TRUE`, returns a list with `selection`,
#'   and optionally `details` and `summaries`.
#' @export
#' @method predict dr_itr
#' @importFrom stats predict
predict.dr_itr <- function(object, newdata, weights = NULL, ids = NULL,
                           details = FALSE, return_summaries = FALSE, ...) {
  if (!inherits(object, 'dr_itr')) stop('object must be a dr_itr fit.', call. = FALSE)
  if (length(list(...))) stop('Unused arguments passed to predict.dr_itr.', call. = FALSE)
  newdata <- .dr_itr_cast(newdata, object$schema)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(newdata)))
  if (length(ids) != nrow(newdata) || anyNA(ids) ||
      anyDuplicated(as.character(ids)) || any(!nzchar(as.character(ids)))) {
    stop('ids must contain one distinct, nonmissing value per profile.',
         call. = FALSE)
  }
  ids <- as.character(ids)
  w <- .dr_itr_weights(weights, names(object$outcomes))
  if (is.null(weights)) w <- object$weights
  tabs <- .dr_itr_table(object$estimates, object$outcomes, newdata,
                        object$covariates, ids, object$weights,
                        use_oob = FALSE)
  result <- stats::predict(object$rank_rule, tabs, weights = w, details = details)
  if (!details && !return_summaries) return(result)
  if (!details) result <- list(selection = result)
  if (return_summaries) result$summaries <- tabs
  result
}

#' Print an integrated DR ITR model
#'
#' @param x A fitted `dr_itr` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.dr_itr <- function(x, ...) {
  cat('Doubly robust individualized treatment ranking model\n')
  cat('  Training profiles: ', x$n_train, '\n', sep = '')
  cat('  Treatments: ', paste(x$treatments, collapse = ', '), '\n', sep = '')
  cat('  Outcomes: ', paste(names(x$outcomes), collapse = ', '), '\n', sep = '')
  cat('  Nuisance folds: ', x$nuisance_folds, '\n', sep = '')
  cat('  Covariance scaling: ', x$covariance_source, '\n', sep = '')
  invisible(x)
}
