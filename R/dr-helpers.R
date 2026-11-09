# Utilities for conditional AIPW estimation across outcome types.

.dr_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install the optional package '", pkg, "' for this learner.", call. = FALSE)
  }
}

.dr_clip <- function(x) pmin(pmax(x, 0), 1)

# Optional rounding applies only to the second-stage ranger response.
# Stored AIPW pseudo-outcomes and first-stage predictions remain unrounded.
.dr_validate_forest_digits <- function(digits, second_stage) {
  if (is.null(digits)) return(NULL)
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0 || digits > 15 ||
      digits != as.integer(digits)) {
    stop('second_stage_round_digits must be NULL or an integer between 0 and 15.',
         call. = FALSE)
  }
  if (second_stage != 'ranger') {
    stop('second_stage_round_digits requires second_stage="ranger".', call. = FALSE)
  }
  as.integer(digits)
}

.dr_response_for_second_stage <- function(phi, digits) {
  if (is.null(digits)) phi else round(phi, digits = digits)
}

.dr_check_data <- function(data, treatment, covariates, outcomes = character()) {
  if (!is.data.frame(data) || nrow(data) < 2L) {
    stop('data must be a data frame with at least two observations.')
  }
  if (!is.character(treatment) || length(treatment) != 1L || is.na(treatment)) {
    stop('treatment must name one column.')
  }
  if (!length(covariates) || anyDuplicated(covariates) ||
      !is.character(covariates) || anyNA(covariates)) {
    stop('covariates must be a nonempty vector of distinct column names.')
  }
  if (any(covariates %in% c('.target', '.treatment', '.trt', '.time', '.censor'))) {
    stop('Internal reserved covariate names are not supported.')
  }
  needed <- unique(c(treatment, covariates, outcomes))
  if (any(!needed %in% names(data))) {
    stop('Missing columns: ', paste(setdiff(needed, names(data)), collapse = ', '))
  }
  if (anyNA(data[, needed, drop = FALSE])) {
    stop('Missing values in the treatment, covariates, or outcomes; preprocess before fitting.')
  }
  if (any(!stats::complete.cases(data[, needed, drop = FALSE])) ||
      any(vapply(data[, needed, drop = FALSE], function(v) {
        is.numeric(v) && any(!is.finite(v))
      }, logical(1)))) {
    stop('Non-finite or missing values in analysis columns.')
  }
  a <- as.character(data[[treatment]])
  levels <- if (is.factor(data[[treatment]])) levels(droplevels(data[[treatment]])) else sort(unique(a))
  if (length(levels) < 2L || any(!nzchar(levels)) || anyNA(a)) {
    stop('At least two nonmissing treatment levels are required.')
  }
  levels
}

.dr_check_newdata <- function(newdata, covariates) {
  if (!is.data.frame(newdata) || nrow(newdata) == 0L) {
    stop('newdata must be a nonempty data frame.')
  }
  if (any(!covariates %in% names(newdata)) ||
      anyNA(newdata[, covariates, drop = FALSE])) {
    stop('newdata must contain complete values for all covariates.')
  }
  invisible(TRUE)
}

.dr_make_folds <- function(treatment, nfolds, seed) {
  if (nfolds == 1L) return(rep(1L, length(treatment)))
  if (nfolds < 2L || nfolds > min(table(treatment))) {
    stop('folds must be 1 or between 2 and the smallest treatment-arm count.')
  }
  set.seed(seed)
  out <- integer(length(treatment))
  for (a in unique(as.character(treatment))) {
    idx <- which(as.character(treatment) == a)
    out[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  out
}

.dr_validate_folds <- function(fold_id, a) {
  if (length(fold_id) != length(a) || anyNA(fold_id)) {
    stop('fold_id must contain one nonmissing fold label per training row.')
  }
  fold_id <- as.integer(factor(fold_id))
  if (length(unique(fold_id)) == 1L) return(rep(1L, length(a)))
  for (f in unique(fold_id)) {
    if (!setequal(unique(a[fold_id != f]), unique(a))) {
      stop('Every nuisance-training split must contain all treatment arms.')
    }
  }
  fold_id
}

.dr_fit_model <- function(y, x, covariates, learner, seed, num.trees, min.node.size,
                          binary = FALSE) {
  if (length(unique(y)) == 1L) return(list(kind = 'constant', value = as.numeric(y[1])))
  df <- data.frame(.target = y, x[, covariates, drop = FALSE], check.names = FALSE)
  form <- stats::reformulate(covariates, response = '.target')
  if (learner == 'ranger') {
    .dr_require('ranger')
    if (binary) df$.target <- factor(df$.target, levels = c(0, 1))
    mod <- ranger::ranger(form, data = df, probability = binary,
                          num.trees = num.trees, min.node.size = min.node.size,
                          respect.unordered.factors = 'order', seed = seed)
  } else if (learner == 'glm' && binary) {
    mod <- stats::glm(form, family = stats::binomial(), data = df)
  } else if (learner == 'lm' && !binary) {
    mod <- stats::lm(form, data = df)
  } else {
    stop('Invalid learner for transformed outcome type.')
  }
  list(kind = learner, model = mod, binary = binary)
}

.dr_predict_model <- function(fit, newdata, covariates) {
  if (fit$kind == 'constant') return(rep(fit$value, nrow(newdata)))
  d <- newdata[, covariates, drop = FALSE]
  if (fit$kind == 'ranger') {
    p <- stats::predict(fit$model, data = d)$predictions
    if (isTRUE(fit$binary)) {
      p <- as.matrix(p)
      if (!'1' %in% colnames(p)) stop('Probability forest did not return the favorable class.')
      p <- p[, '1']
    }
  } else if (fit$kind == 'glm') {
    p <- stats::predict(fit$model, newdata = d, type = 'response')
  } else {
    p <- stats::predict(fit$model, newdata = d)
  }
  if (length(p) != nrow(newdata) || any(!is.finite(p))) {
    stop('Nuisance prediction failed; check factor levels, sparse folds and covariate support.')
  }
  as.numeric(p)
}

.dr_ps_input <- function(propensity, data, treatment, covariates) {
  folds <- NULL
  if (inherits(propensity, 'dr_propensity')) {
    if (propensity$treatment != treatment || !identical(propensity$covariates, covariates)) {
      stop('Propensity object uses different treatment or covariate columns.')
    }
    if (!identical(propensity$row_ids, rownames(data))) {
      stop('Propensity object training-row order differs from data.')
    }
    p <- propensity$ehat
    folds <- propensity$fold_id
  } else {
    p <- propensity
  }
  p <- as.matrix(p)
  if (!is.numeric(p) || nrow(p) != nrow(data) || any(!is.finite(p)) ||
      any(p <= 0) || any(p > 1) || is.null(colnames(p)) ||
      anyDuplicated(colnames(p))) {
    stop('propensity must be a finite, positive n-by-K matrix with treatment-column names.')
  }
  a <- as.character(data[[treatment]])
  if (!setequal(unique(a), colnames(p))) {
    stop('Propensity treatment column names must match all observed treatment levels.')
  }
  p <- p[, colnames(p), drop = FALSE]
  if (!is.null(folds)) folds <- .dr_validate_folds(folds, a)
  list(ehat = p, fold_id = folds, levels = colnames(p))
}

#' Construct the conditional AIPW pseudo-outcome for one treatment
#'
#' Implements m(X) + I(A=k) / e_k(X) * (H - m(X)). This conditional
#' pseudo-outcome is doubly robust with respect to treatment propensity
#' and outcome regression **if** exchangeability and other assumptions hold.
#' Survival analysis also requires valid censoring adjustment.
#'
#' @param transformed Transformed observed outcome (indicator or IPCW response).
#' @param received Character vector of treatments actually received.
#' @param target One treatment label.
#' @param propensity Positive propensity for the target treatment, per row.
#' @param outcome_prediction First-stage prediction m(X), per row.
#' @return Numeric pseudo-outcome vector, potentially outside the unit interval.
#' @export
#' @examples
#' dr_pseudo_outcome(c(1, 0), c('A', 'B'), 'A', c(.5, .5), c(.4, .4))
dr_pseudo_outcome <- function(transformed, received, target, propensity,
                              outcome_prediction) {
  n <- length(transformed)
  if (length(received) != n || length(propensity) != n ||
      length(outcome_prediction) != n || any(!is.finite(transformed)) ||
      any(!is.finite(outcome_prediction)) || any(!is.finite(propensity)) ||
      any(propensity <= 0 | propensity > 1) || anyNA(received) || anyNA(target)) {
    stop('Pseudo-outcome inputs must have equal lengths and finite predictions; propensity in (0, 1].')
  }
  outcome_prediction + (as.character(received) == as.character(target)) /
    propensity * (transformed - outcome_prediction)
}

.dr_components <- function(data, treatment, covariates, propensity, newdata,
                           H, labels, binary = TRUE, outcome_learner = 'ranger',
                           second_stage = 'ranger', fold_id = NULL, seed = 2026L,
                           num.trees = 500L, min.node.size = 10L, clip = TRUE,
                           monotonic = c('none', 'increasing', 'decreasing'),
                           ids = NULL, return_models = FALSE,
                           clip_first_stage = FALSE,
                           seed_schedule = c('package', 'research'),
                           first_stage_array = NULL,
                           second_stage_round_digits = NULL) {
  monotonic <- match.arg(monotonic)
  seed_schedule <- match.arg(seed_schedule)
  second_stage_round_digits <- .dr_validate_forest_digits(
    second_stage_round_digits, second_stage)
  .dr_check_data(data, treatment, covariates)
  .dr_check_newdata(newdata, covariates)
  ps <- .dr_ps_input(propensity, data, treatment, covariates)
  # Original standalone DR functions re-floor and renormalize the supplied
  # estimated propensity matrix. Preserve that behavior only when explicitly
  # requesting the file-30 research schedule and truncation convention.
  if (seed_schedule == 'research' && inherits(propensity, 'dr_propensity') &&
      identical(propensity$truncation, 'floor_renorm')) {
    ps$ehat <- pmax(ps$ehat, propensity$eps)
    ps$ehat <- ps$ehat / rowSums(ps$ehat)
  }
  if (is.null(fold_id)) fold_id <- ps$fold_id
  if (is.null(fold_id)) fold_id <- rep(1L, nrow(data))
  fold_id <- .dr_validate_folds(fold_id, as.character(data[[treatment]]))
  if (length(dim(H)) != 2L || nrow(H) != nrow(data) || ncol(H) != length(labels) ||
      any(!is.finite(H))) stop('H must be a finite training-row by component matrix.')
  if (binary && any(!(H %in% c(0, 1)))) {
    stop('Binary transformed responses must be 0 or 1.')
  }
  if (is.null(ids)) ids <- seq_len(nrow(newdata))
  if (length(ids) != nrow(newdata) || anyNA(ids) || anyDuplicated(ids)) {
    stop('ids must uniquely identify each newdata row.')
  }
  levels <- ps$levels
  n <- nrow(data)
  J <- ncol(H)
  K <- length(levels)
  if (!is.null(first_stage_array) &&
      (!is.array(first_stage_array) ||
       !identical(as.integer(dim(first_stage_array)), c(n, K, J)) ||
       any(!is.finite(first_stage_array)))) {
    stop('first_stage_array must be a finite training-row by treatment by component array.')
  }
  preds <- array(NA_real_, c(nrow(newdata), K, J),
                 dimnames = list(new_id = as.character(ids),
                                 treatment = levels, summary = labels))
  models <- if (return_models) vector('list', K * J) else NULL
  for (s in seq_len(J)) {
    response <- H[, s]
    for (k in seq_len(K)) {
      mhat <- rep(NA_real_, n)
      arm <- as.character(data[[treatment]]) == levels[k]
      if (!is.null(first_stage_array)) {
        mhat <- as.numeric(first_stage_array[, k, s])
      } else {
        for (f in unique(fold_id)) {
          training <- if (length(unique(fold_id)) == 1L) rep(TRUE, n) else fold_id != f
          testing <- fold_id == f
          first_seed <- if (seed_schedule == 'research') {
            seed + 1000L * s + k +
              (if (length(unique(fold_id)) == 1L) 0L else 100000L * f)
          } else seed + s * 1000L + k * 10L + f
          if (seed_schedule == 'research') set.seed(as.integer(first_seed))
          fit <- .dr_fit_model(response[training & arm],
                               data[training & arm, , drop = FALSE], covariates,
                               learner = outcome_learner, seed = first_seed,
                               num.trees = num.trees, min.node.size = min.node.size,
                               binary = binary)
          pred <- .dr_predict_model(fit, data[testing, , drop = FALSE], covariates)
          if (binary || clip_first_stage) pred <- .dr_clip(pred)
          mhat[testing] <- pred
        }
      }
      phi <- dr_pseudo_outcome(response, data[[treatment]], levels[k],
                               ps$ehat[, levels[k]], mhat)
      second_seed <- if (seed_schedule == 'research') {
        seed + 2000L * s + k
      } else seed + 20000L + s * 1000L + k
      # Retain exact AIPW values for audit; stabilize only forest training input.
      phi_fit <- .dr_response_for_second_stage(phi, second_stage_round_digits)
      if (seed_schedule == 'research') set.seed(as.integer(second_seed))
      final_fit <- .dr_fit_model(phi_fit, data, covariates, second_stage,
                                 second_seed, num.trees, min.node.size, binary = FALSE)
      value <- .dr_predict_model(final_fit, newdata, covariates)
      preds[, k, s] <- if (clip) .dr_clip(value) else value
      if (return_models) {
        models[[(s - 1L) * K + k]] <- list(final_fit = final_fit,
                                             pseudo_outcome = phi,
                                             pseudo_outcome_for_fit = if (!is.null(second_stage_round_digits)) phi_fit else NULL,
                                             first_stage_prediction = mhat,
                                             second_seed = second_seed)
        names(models)[(s - 1L) * K + k] <- paste(labels[s], levels[k], sep = '__')
      }
    }
  }
  if (monotonic != 'none' && J > 1L) {
    for (i in seq_len(dim(preds)[1L])) for (k in seq_len(K)) {
      if (monotonic == 'increasing') {
        preds[i, k, ] <- cummax(preds[i, k, ])
      } else {
        preds[i, k, ] <- cummin(preds[i, k, ])
      }
    }
  }
  long <- do.call(rbind, lapply(seq_len(dim(preds)[1L]), function(i) {
    do.call(rbind, lapply(seq_len(J), function(s) {
      data.frame(new_id = as.character(ids[i]), summary_id = labels[s],
                 treatment = levels, theta_hat = as.numeric(preds[i, , s]),
                 stringsAsFactors = FALSE)
    }))
  }))
  rownames(long) <- NULL
  result <- list(theta_array = preds, theta_long = long,
                 treatment_levels = levels, summary_ids = labels,
                 nuisance_folds = length(unique(fold_id)),
                 outcome_learner = outcome_learner, second_stage = second_stage,
                 second_stage_round_digits = second_stage_round_digits)
  if (return_models) result$models <- models
  class(result) <- 'dr_summary_estimate'
  result
}
