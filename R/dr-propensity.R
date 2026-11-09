#' Estimate generalized propensity scores for observational treatments
#'
#' Supports multiple treatment options. The default is a probability forest;
#' `"multinomial"` fits a binary logistic model (K=2) or multinomial logit
#' (`nnet` for K>2). With `folds>1`, the returned `ehat` matrix comprises
#' out-of-fold nuisance predictions, **not** the probabilities of a single
#' fitted full-sample model. Pass this object to the DR summary estimators to
#' use the same nuisance folds for treatment-specific outcome regressions.
#'
#' @param data Training data frame, without missing baseline predictors.
#' @param treatment Name of treatment column (two or more levels).
#' @param covariates Character vector of baseline adjustment columns.
#' @param method `"ranger"` probability forest or `"multinomial"` logit.
#' @param folds Number of stratified nuisance folds (1 means in-sample fits).
#' @param fold_id Optional fold identifier per row, overriding generated folds.
#' @param eps Lower probability bound; bounded scores are renormalized.
#' @param num.trees Number of forest trees when `method="ranger"`.
#' @param min.node.size Minimum forest terminal-node size.
#' @param seed Seed for fold assignment and random forest fitting.
#' @param truncation `"simplex"` (default) enforces the probability floor
#'   and sums to one; `"floor_renorm"` first applies `pmax(p, eps)`
#'   and then renormalizes (research files 01/30). The latter is not a
#'   strict post-renormalization lower bound.
#' @param seed_schedule `"package"` uses `seed + fold`; `"research"`
#'   uses `seed` for one-fold fits, matching file 01. For multiple folds
#'   it uses `seed + fold` to create distinct forests.
#' @return Object of class `dr_propensity` containing `$ehat`, a training-row
#'   by treatment matrix; `$fold_id`; learner settings; and the treatment levels.
#'   Probabilities are truncated, so propensity-only robustness may be affected.
#' @export
fit_dr_propensity <- function(data, treatment, covariates,
                              method = c('ranger', 'multinomial'),
                              folds = 1L, fold_id = NULL, eps = 0.05,
                              num.trees = 300L, min.node.size = 10L,
                              seed = 2026L,
                              truncation = c('simplex', 'floor_renorm'),
                              seed_schedule = c('package', 'research')) {
  method <- match.arg(method)
  truncation <- match.arg(truncation)
  seed_schedule <- match.arg(seed_schedule)
  levels <- .dr_check_data(data, treatment, covariates)
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) ||
      eps <= 0 || eps >= 1 / length(levels)) {
    stop('eps must lie strictly between 0 and 1/K.')
  }
  if (is.null(fold_id)) {
    if (length(folds) != 1L || !is.numeric(folds) ||
        !is.finite(folds) || folds != as.integer(folds)) stop('folds must be an integer.')
    fold_id <- .dr_make_folds(data[[treatment]], as.integer(folds), seed)
  }
  fold_id <- .dr_validate_folds(fold_id, as.character(data[[treatment]]))
  actual <- as.character(data[[treatment]])
  n <- nrow(data)
  preds <- matrix(NA_real_, nrow = n, ncol = length(levels),
                  dimnames = list(NULL, levels))
  for (f in unique(fold_id)) {
    training <- if (length(unique(fold_id)) == 1L) rep(TRUE, n) else fold_id != f
    testing <- fold_id == f
    df <- data.frame(.treatment = factor(actual[training], levels = levels),
                     data[training, covariates, drop = FALSE], check.names = FALSE)
    form <- stats::reformulate(covariates, response = '.treatment')
    if (method == 'ranger') {
      .dr_require('ranger')
      if (seed_schedule == 'research') set.seed(as.integer(seed))
      model <- ranger::ranger(form, data = df, probability = TRUE,
                              num.trees = num.trees, min.node.size = min.node.size,
                              respect.unordered.factors = 'order',
                              seed = if (seed_schedule == 'research' &&
                                         length(unique(fold_id)) == 1L) seed else seed + f)
      p <- stats::predict(model, data = data[testing, covariates, drop = FALSE])$predictions
    } else if (length(levels) == 2L) {
      df$.treatment <- as.integer(df$.treatment == levels[2L])
      model <- stats::glm(form, data = df, family = stats::binomial())
      p1 <- stats::predict(model, newdata = data[testing, covariates, drop = FALSE],
                            type = 'response')
      p <- cbind(1 - p1, p1)
      colnames(p) <- levels
    } else {
      .dr_require('nnet')
      model <- nnet::multinom(form, data = df, trace = FALSE)
      p <- stats::predict(model, newdata = data[testing, covariates, drop = FALSE],
                            type = 'probs')
    }
    # Multinomial predictions for a single validation row may be a vector.
    if (method == 'multinomial' && length(levels) > 2L && is.null(dim(p))) {
      p <- matrix(p, nrow = 1L, dimnames = list(NULL, levels))
    }
    p <- as.matrix(p)
    if (nrow(p) != sum(testing) || !setequal(colnames(p), levels) ||
        any(!is.finite(p))) {
      stop('Propensity estimation returned missing or malformed treatment probabilities.')
    }
    p <- p[, levels, drop = FALSE]
    # Research-compatible alternative: pmax followed by row renormalization.
    if (truncation == 'floor_renorm') {
      p <- pmax(p, eps)
      p <- p / rowSums(p)
      preds[testing, ] <- p
      next
    }
    # Project onto the probability simplex with a genuine eps floor,
    # preserving relative probabilities among the unconstrained arms.
    for (i in seq_len(nrow(p))) {
      raw <- p[i, ]
      if (any(raw < 0) || sum(raw) <= 0) stop('Invalid propensity probabilities.')
      raw <- raw / sum(raw)
      remaining <- seq_along(raw)
      adjusted <- numeric(length(raw))
      mass <- 1
      repeat {
        alloc <- if (sum(raw[remaining]) > 0) {
          mass * raw[remaining] / sum(raw[remaining])
        } else rep(mass / length(remaining), length(remaining))
        low <- alloc < eps
        if (!any(low)) {
          adjusted[remaining] <- alloc
          break
        }
        adjusted[remaining[low]] <- eps
        mass <- mass - eps * sum(low)
        remaining <- remaining[!low]
      }
      preds[which(testing)[i], ] <- adjusted
    }
  }
  structure(list(ehat = preds, fold_id = fold_id, treatment_levels = levels,
                 row_ids = rownames(data),
                 treatment = treatment, covariates = covariates,
                 method = method, eps = eps, truncation = truncation,
                 seed_schedule = seed_schedule, num.trees = num.trees,
                 min.node.size = min.node.size), class = 'dr_propensity')
}
