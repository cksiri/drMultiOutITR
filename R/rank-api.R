# Public entry points for the covariance-scaled rank aggregation layer.
# This is a summary-prediction-based interface; it does NOT estimate DR nuisances.

# Create a canonical description of all outcome-summary components.
.rank_signature <- function(tab) {
  required <- c('outcome_id', 'summary_id', 'larger_is_better',
                'rho', 'summary_weight', 'summary_scale')
  if (!is.logical(tab$larger_is_better) || anyNA(tab$larger_is_better)) {
    stop('larger_is_better must be a nonmissing logical (TRUE/FALSE) column.')
  }
  keys <- paste(as.character(tab$outcome_id), as.character(tab$summary_id), sep = '\r')
  groups <- split(seq_len(nrow(tab)), keys)
  records <- lapply(groups, function(idx) {
    z <- tab[idx, required, drop = FALSE]
    if (any(vapply(z, function(column) length(unique(column)) != 1L, logical(1)))) {
      stop('Each outcome_id/summary_id must have a single direction, rho, ',
           'summary_weight and summary_scale across profiles and treatments.')
    }
    z[1L, , drop = FALSE]
  })
  result <- do.call(rbind, records)
  result <- result[order(as.character(result$outcome_id),
                         as.character(result$summary_id)), , drop = FALSE]
  rownames(result) <- NULL
  result$outcome_id <- as.character(result$outcome_id)
  result$summary_id <- as.character(result$summary_id)
  result
}

# Original shape check is local to each profile. Also enforce a common
# component set across profiles, necessary for frozen fit/predict semantics.
.check_rank_table <- function(tab, treatment_order = NULL, outcome_order = NULL) {
  tab <- validate_matrix_summary_table(tab)
  check_matrix_summary_structure(tab, treatment_order, outcome_order)
  sig <- .rank_signature(tab)
  expected <- paste(sig$outcome_id, sig$summary_id, sep = '\r')
  ids <- unique(tab$new_id)
  for (id in ids) {
    this <- tab[tab$new_id == id, , drop = FALSE]
    actual <- unique(paste(this$outcome_id, this$summary_id, sep = '\r'))
    if (!setequal(expected, actual)) {
      stop('Profile ', id, ' does not have exactly the expected summary components.')
    }
  }
  tab
}

#' Fit the covariance-scaled treatment ranking rule
#'
#' Fit the rank-aggregation stage from already estimated counterfactual
#' summary predictions. The covariance scaling matrices are estimated
#' **only** from the supplied training profiles, separately for each candidate
#' treatment ranking. The function does not fit propensity scores, censoring
#' models, or doubly robust outcome regressions.
#'
#' @param training_summaries Long-format `data.frame` with one row per
#'   profile/outcome/summary/treatment. Required columns are `new_id`,
#'   `outcome_id`, `summary_id`, `treatment`, `theta_hat`,
#'   `larger_is_better` (logical), `outcome_weight`, and `rho`. Optional
#'   columns `summary_weight` and `summary_scale` default to 1.
#' @param magnitude_scale Either `"none"` (manuscript default) or
#'   `"summary_scale"` for optional within-summary magnitude scaling.
#' @param ridge Nonnegative relative diagonal ridge applied to each covariance.
#' @param eigen_floor Positive relative eigenvalue lower bound.
#' @param ties.method Ranking tie method accepted by [base::rank()].
#' @return An object of class `dr_rank_rule` with frozen training-derived
#'   candidate-specific covariance matrices; use [predict()] to apply it.
#' @export
#' @examples
#' # See inst/examples/ranking_from_summaries.R for reproducible synthetic inputs.
fit_rank_rule <- function(training_summaries,
                          magnitude_scale = c('none', 'summary_scale'),
                          ridge = 1e-6, eigen_floor = 1e-8,
                          ties.method = 'average') {
  magnitude_scale <- match.arg(magnitude_scale)
  training_summaries <- .check_rank_table(training_summaries)
  if (length(unique(training_summaries$new_id)) < 2L) {
    stop('At least two training profiles are needed for covariance scaling.')
  }
  k <- length(get_matrix_treatment_order(training_summaries))
  if (k < 2L) stop('At least two treatment alternatives are required.')
  if (k > 6L) warning('Full permutation enumeration scales as K!; this may be costly.')
  fit <- estimate_smmr_scaling_matrices(
    scaling_summary_table = training_summaries,
    magnitude_scale = magnitude_scale,
    ridge = ridge, eigen_floor = eigen_floor,
    ties.method = ties.method, check_structure = FALSE
  )
  structure(
    list(scaling_fit = fit,
         signature = .rank_signature(training_summaries),
         treatment_order = fit$treatment_order,
         outcome_order = fit$outcome_order,
         training_outcome_weights = fit$outcome_weights,
         ties.method = ties.method),
    class = 'dr_rank_rule'
  )
}

#' Predict individualized treatment rankings from a fitted rule
#'
#' Apply a previously fitted covariance-scaled ranking rule to held-out
#' counterfactual summary predictions. This function does not refit the
#' covariance matrices. Changing clinical outcome weights only changes the
#' final aggregation, not the frozen covariance estimates.
#'
#' @param object A `dr_rank_rule` object returned by [fit_rank_rule()].
#' @param newdata Long-format summary predictions with the same outcome,
#'   summary, and treatment structure as the training data.
#' @param weights Optional **named** nonnegative numeric vector of outcome
#'   weights. Names must match all outcome IDs. By default, uses weights
#'   specified when fitting the ranking rule.
#' @param details If `TRUE`, return a list including per-candidate scores;
#'   otherwise return the recommendation data frame.
#' @param ... Reserved for future methods; not currently used.
#' @return By default a data frame with `new_id`, `best_treatment`,
#'   `optimal_ranking`, `min_W`, and `candidate_id` (the latter indexes the
#'   candidate permutation, not a treatment). If `details=TRUE`, a list with
#'   `selection` and `details` data frames.
#' @export
#' @method predict dr_rank_rule
#' @importFrom stats predict
predict.dr_rank_rule <- function(object, newdata, weights = NULL,
                                 details = FALSE, ...) {
  if (!inherits(object, 'dr_rank_rule')) stop('object must be a dr_rank_rule.')
  sf <- object$scaling_fit
  newdata <- .check_rank_table(newdata, sf$treatment_order, sf$outcome_order)
  if (!setequal(get_matrix_treatment_order(newdata), sf$treatment_order)) {
    stop('Treatment choices differ between fit and newdata.')
  }
  if (!setequal(as.character(unique(newdata$outcome_id)), sf$outcome_order)) {
    stop('Outcomes differ between fit and newdata.')
  }
  if (!identical(.rank_signature(newdata), object$signature)) {
    stop('newdata summary components or metadata differ from training data.')
  }
  if (is.null(weights)) {
    lambda <- object$training_outcome_weights
  } else {
    if (!is.numeric(weights) || is.null(names(weights)) ||
        anyDuplicated(names(weights)) ||
        !setequal(names(weights), sf$outcome_order) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop('weights must be a finite, nonnegative named vector for each outcome.')
    }
    lambda <- weights[match(sf$outcome_order, names(weights))]
  }
  if (all(lambda == 0)) stop('At least one outcome weight must be positive.')
  lambda_half <- sqrt(as.numeric(lambda))

  ids <- unique(newdata$new_id)
  candidates <- sf$candidates
  pick <- vector('list', length(ids))
  detail_rows <- if (isTRUE(details)) vector('list', length(ids) * nrow(candidates)) else NULL
  z <- 0L
  for (ii in seq_along(ids)) {
    id <- ids[ii]
    d <- newdata[newdata$new_id == id, , drop = FALSE]
    W <- numeric(nrow(candidates))
    for (cc in seq_len(nrow(candidates))) {
      delta <- as.integer(candidates[cc, ])
      pi <- compute_pi_one_patient(
        summary_table_one = d, delta = delta,
        treatment_order = sf$treatment_order,
        outcome_order = sf$outcome_order,
        magnitude_scale = sf$magnitude_scale,
        ties.method = object$ties.method
      )
      weighted_pi <- lambda_half * as.numeric(pi)
      Sigma <- sf$scaling[[cc]]$Sigma
      inv_times_pi <- tryCatch(solve(Sigma, weighted_pi),
                               error = function(e) qr.solve(Sigma, weighted_pi))
      score <- sum(weighted_pi * inv_times_pi)
      if (!is.finite(score)) score <- Inf
      W[cc] <- score
      if (isTRUE(details)) {
        z <- z + 1L
        row <- data.frame(
          new_id = id, candidate_id = cc,
          candidate_ranking = paste(delta, collapse = '-'),
          W = score,
          condition_number = sf$scaling[[cc]]$condition_number,
          stringsAsFactors = FALSE
        )
        for (jj in seq_along(sf$outcome_order)) {
          row[[paste0('pi_', sf$outcome_order[jj])]] <- pi[jj]
        }
        detail_rows[[z]] <- row
      }
    }
    if (all(!is.finite(W))) stop('All ranking scores are nonfinite for profile ', id)
    best <- which.min(W)
    delta <- as.integer(candidates[best, ])
    pick[[ii]] <- data.frame(
      new_id = id,
      best_treatment = sf$treatment_order[which.min(delta)],
      optimal_ranking = paste(delta, collapse = '-'),
      min_W = W[best], candidate_id = best,
      stringsAsFactors = FALSE
    )
  }
  selection <- do.call(rbind, pick)
  rownames(selection) <- NULL
  if (!isTRUE(details)) return(selection)
  detail_table <- do.call(rbind, detail_rows)
  rownames(detail_table) <- NULL
  list(selection = selection, details = detail_table)
}

#' Apply the original covariance-scaled ranking implementation
#'
#' A direct compatibility interface to the original manuscript rank engine.
#' For repeated held-out predictions, prefer [fit_rank_rule()] followed by
#' [predict()] to reuse the fitted covariance scaling matrices.
#'
#' @param summaries Long-format counterfactual summary predictions to rank.
#' @param training_summaries Long-format predictions on the training profiles,
#'   used exclusively for covariance scaling.
#' @param ... Additional parameters passed to the original rank engine,
#'   e.g. `magnitude_scale`, `ridge`, `eigen_floor`, `return_details`.
#' @return The original rank-engine list with `selection`, `scaling_fit`,
#'   and, by default, `details`.
#' @export
rank_treatments <- function(summaries, training_summaries, ...) {
  summaries <- .check_rank_table(summaries)
  training_summaries <- .check_rank_table(training_summaries)
  select_rank_based_treatment_matrix(
    summary_table = summaries,
    scaling_summary_table = training_summaries, ...
  )
}
