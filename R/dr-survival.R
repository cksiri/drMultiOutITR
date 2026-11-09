# Censoring models and treatment-specific survival nuisance estimation.

.dr_censor_fit <- function(data, treatment, covariates, time, status, method) {
  .dr_require('survival')
  if (!any(status == 0L)) return(list(kind = 'none'))
  if (method == 'pooled_km') {
    mod <- survival::survfit(survival::Surv(time, as.integer(1 - status)) ~ 1)
    return(list(kind = method, fit = mod))
  }
  df <- data.frame(.time = time, .censor = as.integer(1 - status),
                   .trt = factor(as.character(data[[treatment]])),
                   data[, covariates, drop = FALSE], check.names = FALSE)
  form <- stats::reformulate(c('.trt', covariates),
                            response = 'survival::Surv(.time, .censor)')
  fit <- survival::coxph(form, data = df, ties = 'breslow', model = TRUE)
  if (anyNA(stats::coef(fit))) {
    stop('Conditional censoring Cox model is singular; simplify predictors or use a suitable alternative.')
  }
  list(kind = 'cox', fit = fit,
       basehaz = survival::basehaz(fit, centered = FALSE))
}

.dr_censor_predict <- function(fit, newdata, treatment, covariates, times, min_G) {
  n <- nrow(newdata)
  if (fit$kind == 'none') return(matrix(1, n, length(times)))
  if (fit$kind == 'pooled_km') {
    G <- as.numeric(summary(fit$fit, times = times, extend = TRUE)$surv)
    if (length(G) != length(times)) stop('Could not predict pooled censoring survival.')
    return(matrix(rep(pmax(G, min_G), each = n), nrow = n))
  }
  b <- fit$basehaz
  j <- findInterval(times, b$time)
  haz <- numeric(length(times))
  haz[j > 0L] <- b$hazard[j[j > 0L]]
  df <- data.frame(.trt = factor(as.character(newdata[[treatment]]),
                                  levels = fit$fit$xlevels$.trt),
                   newdata[, covariates, drop = FALSE], check.names = FALSE)
  risk <- as.numeric(stats::predict(fit$fit, newdata = df, type = 'risk',
                                    reference = 'zero'))
  if (length(risk) != n || any(!is.finite(risk))) {
    stop('Conditional censoring probabilities could not be computed.')
  }
  pmax(exp(-outer(risk, haz)), min_G)
}

# Treatment-specific event-time forests are distinct from regressions of
# IPCW-transformed indicators.
.dr_rsf_first_stage <- function(data, treatment, time, status, covariates,
                                times, fold_id, seed, ntree, nodesize,
                                return_models = FALSE) {
  .dr_require('randomForestSRC')
  .dr_require('survival')
  levels <- unique(as.character(data[[treatment]]))
  # Use the same treatment ordering as the propensity prediction object.
  levels <- if (is.factor(data[[treatment]])) base::levels(data[[treatment]]) else sort(levels)
  n <- nrow(data); K <- length(levels); TT <- length(times)
  result <- array(NA_real_, c(n, K, TT),
                  dimnames = list(NULL, levels, paste0('t', times)))
  fits <- if (return_models) list() else NULL
  for (f in unique(fold_id)) {
    train <- if (length(unique(fold_id)) == 1L) rep(TRUE, n) else fold_id != f
    test <- fold_id == f
    for (k in seq_along(levels)) {
      in_arm <- train & as.character(data[[treatment]]) == levels[k]
      if (!any(in_arm)) stop('No observations for treatment ', levels[k], ' in RSF training fold.')
      pred <- matrix(1, sum(test), TT)
      if (any(status[in_arm] == 1L)) {
        d <- as.data.frame(data[in_arm, covariates, drop = FALSE])
        d$.dr_time <- as.numeric(time[in_arm])
        d$.dr_status <- as.integer(status[in_arm])
        form <- stats::as.formula('Surv(.dr_time, .dr_status) ~ .',
                                   env = asNamespace('survival'))
        forest_seed <- as.integer(seed + 1000L + k +
                        if (length(unique(fold_id)) == 1L) 0L else 100000L * f)
        model <- randomForestSRC::rfsrc(
          formula = form, data = d, ntree = as.integer(ntree),
          nodesize = as.integer(nodesize), forest = TRUE,
          seed = forest_seed, verbose = FALSE)
        forest_pred <- stats::predict(
          model, newdata = as.data.frame(data[test, covariates, drop = FALSE]))
        surv <- as.matrix(forest_pred$survival)
        grid <- forest_pred$time.interest
        if (is.null(grid)) grid <- model$time.interest
        if (is.null(grid) || nrow(surv) != sum(test) ||
            ncol(surv) != length(grid)) {
          stop('RSF returned incompatible survival curve and time grid.')
        }
        grid <- as.numeric(grid)
        for (tt in seq_along(times)) {
          j <- which(grid <= times[tt])
          if (length(j)) pred[, tt] <- surv[, max(j)]
        }
        if (return_models) fits[[paste0('fold', f, '__', levels[k])]] <- model
      } else if (return_models) {
        fits[[paste0('fold', f, '__', levels[k])]] <- list(kind = 'constant_one')
      }
      result[test, k, ] <- pmax(0, pmin(1, pred))
    }
  }
  if (any(!is.finite(result))) stop('RSF nuisance predictions are not finite.')
  list(mhat = result, models = fits)
}

#' Estimate doubly robust conditional survival probabilities
#'
#' Constructs IPCW transformations `I(U > t) / G(t | A,X)`, uses treatment-
#' specific regression of those transformed responses, computes treatment AIPW
#' pseudo-outcomes and fits a conditional second-stage regression. Conditional
#' independent censoring and censoring positivity are required in addition to
#' usual treatment identification assumptions. **Double robustness here covers
#' the treatment propensity and outcome regression, conditional on a valid
#' censoring estimator**; it does not protect against misspecified censoring.
#'
#' @inheritParams estimate_dr_binary
#' @param time Name of observed follow-up time column.
#' @param status Name of indicator column (1 = observed event, 0 = censored).
#' @param times Positive, increasing evaluation times.
#' @param censoring `"pooled_km"` requires censoring independent of covariates
#'   and treatment; `"cox"` models censoring conditional on treatment and X.
#' @param min_G Positive lower bound for estimated censoring survival, used
#'   for numerical stabilization (can introduce bias).
#' @param outcome_learner `"ranger"` or `"lm"` regress the IPCW transformed
#'   response in each treatment arm. `"rsf"` fits treatment-specific event-time
#'   random survival forests, as used in the final simulation (file 30).
#' @param second_stage Numeric pseudo-outcome learner: `"ranger"` or `"lm"`.
#' @param rsf_num.trees,rsf_nodesize Random-survival-forest first-stage settings;
#'   default to the corresponding numeric forest parameters.
#' @param seed_schedule `"research"` follows files 04/29 seed offsets;
#'   `"package"` retains the previous package behavior.
#' @return `dr_summary_estimate` with conditional survival summaries at `times`,
#'   censoring model information and the fraction of small G estimates.
#' @export
estimate_dr_survival <- function(data, treatment, time, status, covariates,
                                 propensity, newdata = data, times,
                                 censoring = c('pooled_km', 'cox'), min_G = .05,
                                 outcome_learner = c('ranger', 'lm', 'rsf'),
                                 second_stage = c('ranger', 'lm'),
                                 fold_id = NULL, seed = 2026L,
                                 num.trees = 500L, min.node.size = 10L,
                                 ids = NULL, return_models = FALSE,
                                 seed_schedule = c('package', 'research'),
                                 rsf_num.trees = num.trees,
                                 rsf_nodesize = min.node.size,
                                 second_stage_round_digits = NULL) {
  censoring <- match.arg(censoring)
  seed_schedule <- match.arg(seed_schedule)
  outcome_learner <- match.arg(outcome_learner)
  second_stage <- match.arg(second_stage)
  .dr_check_data(data, treatment, covariates, c(time, status))
  .dr_check_newdata(newdata, covariates)
  t_obs <- data[[time]]
  d <- data[[status]]
  if (!is.numeric(t_obs) || any(!is.finite(t_obs)) || any(t_obs < 0) ||
      any(!d %in% c(0, 1))) {
    stop('Observed follow-up must be nonnegative numeric and status must be 0/1.')
  }
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
      any(times <= 0) || is.unsorted(times, strictly = TRUE)) {
    stop('times must be strictly increasing and positive.')
  }
  if (length(min_G) != 1L || !is.finite(min_G) || min_G <= 0 || min_G > 1) {
    stop('min_G must lie in (0,1].')
  }
  ps <- .dr_ps_input(propensity, data, treatment, covariates)
  if (is.null(fold_id)) fold_id <- ps$fold_id
  if (is.null(fold_id)) fold_id <- rep(1L, nrow(data))
  fold_id <- .dr_validate_folds(fold_id, as.character(data[[treatment]]))
  G <- matrix(NA_real_, nrow(data), length(times))
  for (f in unique(fold_id)) {
    train <- if (length(unique(fold_id)) == 1L) rep(TRUE, nrow(data)) else fold_id != f
    test <- fold_id == f
    fit <- .dr_censor_fit(data[train, , drop = FALSE], treatment, covariates,
                          t_obs[train], d[train], method = censoring)
    G[test, ] <- .dr_censor_predict(fit, data[test, , drop = FALSE],
                                    treatment, covariates, times, min_G)
  }
  if (any(!is.finite(G))) stop('Nonfinite censoring survival predictions.')
  H <- vapply(seq_along(times), function(j) {
    as.numeric(t_obs > times[j]) / G[, j]
  }, numeric(nrow(data)))
  if (is.null(dim(H))) H <- matrix(H, ncol = length(times))
  rsf <- NULL
  if (outcome_learner == 'rsf') {
    rsf <- .dr_rsf_first_stage(data, treatment, t_obs, d, covariates,
                               times, fold_id, seed, rsf_num.trees,
                               rsf_nodesize, return_models = return_models)
    # Reorder to match propensity's treatment levels (which control AIPW).
    rsf$mhat <- rsf$mhat[, ps$levels, , drop = FALSE]
  }
  result <- .dr_components(data, treatment, covariates, propensity, newdata,
                           H, labels = paste0('S_t', format(times, trim = TRUE)),
                           binary = FALSE,
                           outcome_learner = if (is.null(rsf)) outcome_learner else 'ranger',
                           second_stage = second_stage, fold_id = fold_id,
                           seed = seed, num.trees = num.trees,
                           min.node.size = min.node.size,
                           monotonic = 'decreasing', ids = ids,
                           return_models = return_models,
                           clip_first_stage = TRUE,
                           seed_schedule = seed_schedule,
                           first_stage_array = if (is.null(rsf)) NULL else rsf$mhat,
                           second_stage_round_digits = second_stage_round_digits)
  if (!is.null(rsf)) {
    result$outcome_learner <- 'rsf'
    if (return_models) result$rsf_models <- rsf$models
  }
  result$times <- times
  result$censoring <- censoring
  result$min_G <- min_G
  result$censoring_G_train <- G
  result$G_at_floor_fraction <- mean(G <= min_G + .Machine$double.eps)
  result
}
