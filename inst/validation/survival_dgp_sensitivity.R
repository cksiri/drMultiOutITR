# Paired lognormal/Weibull survival-model sensitivity experiments.
# Requires the separate study reference scripts; sourcing defines functions
# but does not fit models or change package estimators.

.phase4b_require <- function() {
  needed <- c('.phase3_reference_environment', 'run_phase3_final_simulation',
              '.phase3_summary_array', '.phase3_metric',
              '.phase3c_refit_survival_only', '.phase4_aligned_choices')
  absent <- needed[!vapply(needed, exists, logical(1), mode = 'function')]
  if (length(absent)) {
    stop('Source the Phase 3, Phase 3C, and Phase 4 validation runners first. Missing: ',
         paste(absent, collapse = ', '), call. = FALSE)
  }
  if (!exists('fit_dr_itr', mode = 'function')) {
    stop('First run devtools::load_all().', call. = FALSE)
  }
}

# Build DGP cells; paired model rows use exactly the same simulation seed.
make_phase4b_survival_scenarios <- function(design = c('pilot', 'extended'),
                                            weibull_shape = 6) {
  design <- match.arg(design)
  if (length(weibull_shape) != 1L || !is.numeric(weibull_shape) ||
      !is.finite(weibull_shape) || weibull_shape <= 0) {
    stop('weibull_shape must be one positive finite number.')
  }
  n <- if (design == 'pilot') 300L else c(300L, 600L)
  censor <- if (design == 'pilot') c('25', '50') else c('none', '25', '50')
  z <- expand.grid(survival_model = c('lognormal', 'weibull'),
                   censoring = censor, n_total = n,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  z$sigma <- 1
  z$rho <- 0.3
  z$weibull_shape <- weibull_shape
  z$pair_id <- paste0('N', z$n_total, '_s1_rho03_c', z$censoring)
  z$scenario_id <- paste0(z$pair_id, '_', z$survival_model)
  z <- z[, c('scenario_id', 'pair_id', 'survival_model', 'n_total',
             'sigma', 'rho', 'censoring', 'weibull_shape')]
  rownames(z) <- NULL
  z
}

.phase4b_check_scenarios <- function(z) {
  cols <- c('scenario_id', 'pair_id', 'survival_model', 'n_total',
            'sigma', 'rho', 'censoring', 'weibull_shape')
  if (!is.data.frame(z) || !nrow(z) || !all(cols %in% names(z)) ||
      anyNA(z[, cols, drop = FALSE])) {
    stop('Scenarios must be a nonempty data.frame with: ',
         paste(cols, collapse = ', '))
  }
  if (anyDuplicated(z$scenario_id) || any(!nzchar(as.character(z$scenario_id))) ||
      any(!nzchar(as.character(z$pair_id))) ||
      any(!z$survival_model %in% c('lognormal', 'weibull')) ||
      any(!z$censoring %in% c('none', '25', '50')) ||
      !is.numeric(z$n_total) || any(!z$n_total %in% c(300L, 600L)) ||
      !is.numeric(z$sigma) || any(!is.finite(z$sigma) | z$sigma <= 0) ||
      !is.numeric(z$rho) || any(!is.finite(z$rho) | z$rho < -0.5 | z$rho >= 1) ||
      !is.numeric(z$weibull_shape) ||
      any(!is.finite(z$weibull_shape) | z$weibull_shape <= 0)) {
    stop('Invalid Phase 4B scenario values or duplicated scenario IDs.')
  }
  # A pairing means all design settings except the survival family agree.
  for (id in unique(z$pair_id)) {
    a <- z[z$pair_id == id, , drop = FALSE]
    if (anyDuplicated(a$survival_model) ||
        length(unique(a$n_total)) != 1L ||
        length(unique(a$sigma)) != 1L ||
        length(unique(a$rho)) != 1L ||
        length(unique(a$censoring)) != 1L ||
        length(unique(a$weibull_shape)) != 1L) {
      stop('Paired scenarios must share all settings except survival_model: ', id)
    }
  }
  z <- z[, cols, drop = FALSE]
  z$n_total <- as.integer(z$n_total)
  z$censoring <- as.character(z$censoring)
  z$survival_model <- as.character(z$survival_model)
  z$pair_id <- as.character(z$pair_id)
  z$scenario_id <- as.character(z$scenario_id)
  rownames(z) <- NULL
  z
}

.phase4b_jobs <- function(scenarios, B, seed_start) {
  z <- .phase4b_check_scenarios(scenarios)
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) ||
      B != as.integer(B) || B < 1L || B > 10L) {
    stop('B must be an integer from 1 to 10 in this pilot.')
  }
  if (length(seed_start) != 1L || !is.numeric(seed_start) ||
      !is.finite(seed_start) || seed_start != as.integer(seed_start) ||
      seed_start < 1L) {
    stop('seed_start must be a positive integer.')
  }
  pair_ids <- unique(z$pair_id)
  ans <- do.call(rbind, lapply(seq_len(nrow(z)), function(i) {
    data.frame(z[rep(i, B), , drop = FALSE], replicate = seq_len(B))
  }))
  seed_values <- as.double(seed_start) +
    (match(ans$pair_id, pair_ids) - 1) * 10000 + ans$replicate - 1L
  if (any(seed_values + 1000001 > .Machine$integer.max)) {
    stop('Seeds are too large for the independent censoring-calibration stream.')
  }
  ans$seed <- as.integer(seed_values)
  rownames(ans) <- NULL
  ans
}

# Compute the censoring rate solving E[1-exp(-rate*T)] = target on an
# INDEPENDENT calibration simulation (not the training event times).
.phase4b_censor_rate <- function(T_event, target) {
  if (!is.numeric(T_event) || !length(T_event) ||
      any(!is.finite(T_event) | T_event <= 0) ||
      length(target) != 1L || !is.finite(target) ||
      target <= 0 || target >= 1) {
    stop('T_event must be positive finite and target must lie in (0,1).')
  }
  fun <- function(rate) mean(-expm1(-rate * T_event)) - target
  upper <- 1
  while (fun(upper) <= 0) upper <- upper * 2
  stats::uniroot(fun, lower = 0, upper = upper, tol = 1e-11)$root
}

# Transform the SAME Gaussian latent event-time residual used in the
# original lognormal simulator into a Weibull margin. Because this is a
# monotone quantile transform of the original latent score, the dependence
# (Gaussian copula) with the continuous/ordinal outcome errors is preserved.
# Unlike drawing fresh Weibull residuals, it isolates the survival margin.
# Original log(T) mean = f3/d2; its marginal residual SD is analytic.
.phase4b_weibull_times <- function(ref, sim, shape) {
  if (length(shape) != 1L || !is.numeric(shape) ||
      !is.finite(shape) || shape <= 0 ||
      any(!is.finite(sim$T_event) | sim$T_event <= 0)) {
    stop('Invalid Weibull shape or original event times.')
  }
  covs <- grep('^X[0-9]+$', names(sim$dat), value = TRUE)
  x <- sim$dat[, covs, drop = FALSE]
  mu <- as.numeric(ref$smmr_f_mean(3L, sim$dat$A, x)) / sim$d2
  sd_log <- sqrt((sim$sigma / (sim$d1 * sim$d2))^2 + sim$sigma_L^2)
  if (!is.finite(sd_log) || sd_log <= 0) {
    stop('Cannot standardize original lognormal latent survival errors.')
  }
  z <- (log(sim$T_event) - mu) / sd_log
  u <- stats::pnorm(-z)
  u <- pmin(1 - .Machine$double.eps,
            pmax(.Machine$double.eps, u))
  times <- exp(mu) * (-log(u))^(1 / shape)
  if (any(!is.finite(times) | times <= 0)) stop('Invalid Weibull times.')
  times
}

.phase4b_true_weibull_survival <- function(ref, xnew, time_points,
                                           d2, treat_levels, shape) {
  index <- expand.grid(treatment = as.character(treat_levels),
                       time_index = seq_along(time_points),
                       new_id = seq_len(nrow(xnew)),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  xx <- xnew[index$new_id, , drop = FALSE]
  scale <- exp(ref$smmr_f_mean(3L, index$treatment, xx) / d2)
  s <- exp(- (time_points[index$time_index] / scale)^shape)
  data.frame(new_id = index$new_id, time_index = index$time_index,
             time_point = time_points[index$time_index],
             treatment = index$treatment, theta_hat = s,
             stringsAsFactors = FALSE)
}

# Install a modified simulation and (for Weibull) its EXACT analytic
# counterfactual survival oracle into a private sourced research environment.
# Original source code is never edited; no package R/ functions are changed.
.phase4b_install_survival_dgp <- function(ref, survival_model,
                                          censoring_mode, shape = 6,
                                          calibration_n = 5000L) {
  if (!is.environment(ref) ||
      !survival_model %in% c('lognormal', 'weibull') ||
      !censoring_mode %in% c('none', '25', '50') ||
      length(shape) != 1L || !is.finite(shape) || shape <= 0 ||
      length(calibration_n) != 1L ||
      calibration_n != as.integer(calibration_n) ||
      calibration_n < 200L) {
    stop('Invalid survival-DGP configuration.')
  }
  original_sim <- ref$simulate_smmr_original_like_data
  ref$simulate_smmr_original_like_data <- function(...) {
    args <- list(...)
    if (!is.null(args$randomized) && isTRUE(args$randomized)) {
      stop('Phase 4B requires observational, not randomized, treatment.')
    }
    if (!is.null(args$n_per_arm)) {
      stop('Phase 4B requires n_per_arm=NULL (observational sample).')
    }
    seed <- if (is.null(args$seed)) 2026L else as.integer(args$seed)
    sim <- do.call(original_sim, args)
    after_sim_rng <- get('.Random.seed', envir = .GlobalEnv)
    on.exit(assign('.Random.seed', after_sim_rng, envir = .GlobalEnv), add = TRUE)
    times <- if (survival_model == 'weibull') {
      .phase4b_weibull_times(ref, sim, shape)
    } else sim$T_event
    if (censoring_mode == 'none') {
      rate <- 0
      C <- rep(Inf, length(times))
    } else {
      # This independent calibration sample avoids tuning the censor rate
      # on the SAME subjects used for fitting the forest.
      cal_args <- args
      cal_args$n_total <- as.integer(calibration_n)
      cal_args$seed <- as.integer(seed + 1000001L)
      cal <- do.call(original_sim, cal_args)
      T_cal <- if (survival_model == 'weibull') {
        .phase4b_weibull_times(ref, cal, shape)
      } else cal$T_event
      rate <- .phase4b_censor_rate(T_cal, as.numeric(censoring_mode) / 100)
      # Keep the same underlying exponential censoring uniforms for the
      # two survival models; change only their calibrated scale.
      C <- sim$C * (sim$censor_rate / rate)
    }
    sim$T_event <- as.numeric(times)
    sim$C <- as.numeric(C)
    sim$censor_rate <- rate
    sim$censoring <- censoring_mode
    sim$dat$time <- pmin(sim$T_event, sim$C)
    sim$dat$status <- as.integer(sim$T_event <= sim$C)
    sim$time_points <- as.numeric(stats::quantile(
      sim$T_event, probs = c(.25, .5, .75), names = FALSE))
    sim$survival_model <- survival_model
    sim$weibull_shape <- if (survival_model == 'weibull') shape else NA_real_
    sim$target_censoring_fraction <- if (censoring_mode == 'none') 0 else
      as.numeric(censoring_mode) / 100
    sim
  }
  if (survival_model == 'weibull') {
    ref$true_smmr_survival <- function(xnew, time_points, sigma = 1,
                                        d1 = 10, d2 = 5, sigma_L = .2,
                                        treat_levels = c('1', '2', '3')) {
      .phase4b_true_weibull_survival(ref, xnew, time_points, d2,
                                     treat_levels, shape)
    }
  }
  invisible(ref)
}

#' Audit that the two simulated event-time mechanisms share the SAME
#' non-survival observations and base censoring draws, before fitting forests.
#' This is a cheap synthetic-DGP test; it does not run any model fits.
audit_phase4b_dgp_pair <- function(reference_directory, seed = 2026L,
                                   n_total = 300L, sigma = 1, rho = .3,
                                   censoring = c('25', '50', 'none'),
                                   weibull_shape = 6, calibration_n = 5000L) {
  if (!exists('.phase3_reference_environment', mode = 'function')) {
    stop('Source validation/run_final_simulation_phase3.R first.')
  }
  censoring <- match.arg(censoring)
  if (!n_total %in% c(300L, 600L) ||
      !is.finite(seed) || length(seed) != 1L ||
      !is.finite(sigma) || sigma <= 0 || !is.finite(rho) ||
      rho < -.5 || rho >= 1) stop('Invalid audit DGP settings.')
  args <- list(n_total=n_total, p=10L, sigma=sigma, rho=rho,
               censoring=if(censoring == 'none') '25' else censoring,
               seed=as.integer(seed), randomized=FALSE,
               confounding_strength=.75)
  ref_l <- .phase3_reference_environment(reference_directory)
  ref_w <- .phase3_reference_environment(reference_directory)
  .phase4b_install_survival_dgp(ref_l, 'lognormal', censoring,
                                 weibull_shape, calibration_n)
  .phase4b_install_survival_dgp(ref_w, 'weibull', censoring,
                                 weibull_shape, calibration_n)
  l <- do.call(ref_l$simulate_smmr_original_like_data, args)
  w <- do.call(ref_w$simulate_smmr_original_like_data, args)
  cols <- c('A', paste0('X',1:10), 'Ycont','Yord')
  invariant <- identical(l$dat[,cols,drop=FALSE],
                         w$dat[,cols,drop=FALSE])
  if (!invariant) stop('DGPs differ in nonsurvival observations.')
  xx <- l$dat[, paste0('X',1:10), drop=FALSE]
  mu <- as.numeric(ref_l$smmr_f_mean(3L, l$dat$A, xx)) / l$d2
  sd_log <- sqrt((l$sigma/(l$d1*l$d2))^2 + l$sigma_L^2)
  z <- (log(l$T_event) - mu)/sd_log
  implied_z <- stats::qnorm(exp(- (w$T_event/exp(mu))^weibull_shape),
                             lower.tail=FALSE)
  latent_score_max_abs_diff <- max(abs(implied_z - z))
  # Censoring draws may be rescaled but must retain their common uniforms.
  censor_base_diff <- if (censoring == 'none') 0 else
    max(abs(l$C * l$censor_rate - w$C * w$censor_rate))
  if (!is.finite(latent_score_max_abs_diff) ||
      latent_score_max_abs_diff > 1e-8 ||
      !is.finite(censor_base_diff) || censor_base_diff > 1e-8) {
    stop('DGP event-score/censoring coupling failed numerical audit.')
  }
  data.frame(seed=seed, n_total=n_total, censoring=censoring,
    baseline_outcomes_identical=invariant,
    latent_score_max_abs_diff=latent_score_max_abs_diff,
    shared_exponential_censor_draw_max_abs_diff=censor_base_diff,
    lognormal_realized_censoring=mean(l$dat$status == 0),
    weibull_realized_censoring=mean(w$dat$status == 0),
    stringsAsFactors=FALSE)
}

# Perturb ONLY the second-stage training responses at floating-point scale.
# The AIPW pseudo-outcomes are retained UNMODIFIED; compare the resulting
# fit against the saved raw model using its original seeds and nuisance fits.
.phase4b_perturb_survival_only <- function(fit, dat, epsilon_multiple = 2) {
  if (!is.finite(epsilon_multiple) || epsilon_multiple <= 0 ||
      length(epsilon_multiple) != 1L) stop('epsilon_multiple must be positive.')
  out <- fit
  est <- out$estimates$survival
  for (j in seq_along(est$models)) {
    item <- est$models[[j]]
    phi <- as.numeric(item$pseudo_outcome)
    signed <- ifelse((seq_along(phi) + j) %% 2L == 0L, 1, -1)
    phi_fit <- phi + epsilon_multiple * .Machine$double.eps *
      pmax(1, abs(phi)) * signed
    if (length(phi) != nrow(dat) || any(!is.finite(phi_fit))) {
      stop('Invalid perturbation inputs.')
    }
    item$final_fit <- drMultiOutITR:::.dr_fit_model(
      y = phi_fit, x = dat, covariates = out$covariates, learner = 'ranger',
      seed = item$second_seed, num.trees = out$survival_num.trees,
      min.node.size = out$survival_min_node_size, binary = FALSE)
    item$pseudo_outcome_for_fit <- phi_fit
    est$models[[j]] <- item
  }
  new_train <- drMultiOutITR:::.dr_itr_model_predictions(
    est, dat, out$covariates, as.character(seq_len(nrow(dat))),
    use_oob = FALSE)
  est$theta_array <- new_train$theta_array
  est$theta_long <- new_train$theta_long
  out$estimates$survival <- est
  training_table <- drMultiOutITR:::.dr_itr_table(
    out$estimates, out$outcomes, dat, out$covariates,
    as.character(seq_len(nrow(dat))), out$weights, use_oob = FALSE)
  out$rank_rule <- fit_rank_rule(training_table, magnitude_scale = 'none',
                                 ridge = 1e-6, eigen_floor = 1e-8)
  if (!identical(lapply(fit$estimates$survival$models,
                        function(x) x$pseudo_outcome),
                 lapply(out$estimates$survival$models,
                        function(x) x$pseudo_outcome))) {
    stop('Perturbation unexpectedly changed stored AIPW responses.')
  }
  out
}

.phase4b_one_job <- function(job, reference_directory, digits,
                              calibration_n, epsilon_multiple, mode) {
  ref <- .phase3_reference_environment(reference_directory)
  .phase4b_install_survival_dgp(ref, job$survival_model,
    job$censoring, job$weibull_shape, calibration_n)
  # Pass the modified data-generating model through the reference environment.
  one <- run_phase3_final_simulation(
    reference_directory = reference_directory, reference_override = ref,
    mode = mode, reference_survival = 'final_rsf', alignment = 'research',
    seed = job$seed, n_total = job$n_total, sigma = job$sigma,
    rho = job$rho, censoring = if (job$censoring == 'none') '25' else job$censoring)
  orig <- one$original_fit
  pkg <- one$package_fit
  stable <- .phase3c_refit_survival_only(pkg, orig$dat, digits)
  perturbed <- .phase4b_perturb_survival_only(pkg, orig$dat, epsilon_multiple)
  ids <- as.character(seq_len(nrow(orig$xnew)))
  raw <- one$package_predictions
  p_stable <- stats::predict(stable, orig$xnew, ids = ids,
                              details = TRUE, return_summaries = TRUE)
  p_perturb <- stats::predict(perturbed, orig$xnew, ids = ids,
                               details = TRUE, return_summaries = TRUE)
  trts <- as.character(orig$sim$treat_levels)
  a_raw <- .phase3_summary_array(raw$summaries, 'survival', ids, trts)
  a_round <- .phase3_summary_array(p_stable$summaries, 'survival', ids, trts)
  a_eps <- .phase3_summary_array(p_perturb$summaries, 'survival', ids, trts)
  d_round <- .phase3_metric('raw versus rounded', a_raw, a_round)
  d_eps <- .phase3_metric('raw versus machine-epsilon perturbed', a_raw, a_eps)
  ora <- orig$matrix_fits$true$selection
  o_orig <- .phase4_aligned_choices(ora, orig$matrix_fits$estimated$selection)
  o_raw <- .phase4_aligned_choices(ora, raw$selection)
  o_round <- .phase4_aligned_choices(ora, p_stable$selection)
  o_eps <- .phase4_aligned_choices(ora, p_perturb$selection)
  raw_round <- .phase4_aligned_choices(raw$selection, p_stable$selection)
  raw_eps <- .phase4_aligned_choices(raw$selection, p_perturb$selection)
  non_surv <- c(one$comparison$max_abs_diff[
      one$comparison$component %in% c('propensity (estimated)',
        'counterfactual summary: continuous', 'counterfactual summary: ordinal')],
    one$scaling_comparison$max_abs_diff[
      one$scaling_comparison$component %in% c('scaling summary: continuous',
        'scaling summary: ordinal')])
  if (length(non_surv) != 5L || any(!is.finite(non_surv)) ||
      max(non_surv) > 1e-8) {
    stop('Non-survival parity did not hold under paired DGP: ', max(non_surv))
  }
  if (!isTRUE(one$ranking_parity_on_identical_summaries)) {
    stop('Rank parity failed under identical counterfactual summaries.')
  }
  if (!identical(lapply(pkg$estimates$survival$models,
                        function(x) x$pseudo_outcome),
                 lapply(stable$estimates$survival$models,
                        function(x) x$pseudo_outcome))) {
    stop('Rounding modified stored AIPW pseudo-outcomes.')
  }
  data.frame(n_test = length(o_raw$same_best),
    realized_censoring = mean(orig$dat$status == 0),
    target_censoring = orig$sim$target_censoring_fraction,
    oracle_original_correct = sum(o_orig$same_best),
    oracle_raw_correct = sum(o_raw$same_best),
    oracle_rounded_correct = sum(o_round$same_best),
    oracle_perturbed_correct = sum(o_eps$same_best),
    raw_rounded_choice_switches = sum(!raw_round$same_best),
    raw_perturbed_choice_switches = sum(!raw_eps$same_best),
    raw_rounded_full_rank_switches = sum(!raw_round$same_rank),
    raw_perturbed_full_rank_switches = sum(!raw_eps$same_rank),
    raw_rounded_survival_max = d_round$max_abs_diff,
    raw_rounded_survival_mean = d_round$mean_abs_diff,
    raw_perturbed_survival_max = d_eps$max_abs_diff,
    raw_perturbed_survival_mean = d_eps$mean_abs_diff,
    non_survival_max = max(non_surv),
    rank_parity = one$ranking_parity_on_identical_summaries,
    stringsAsFactors = FALSE)
}

.phase4b_aggregate <- function(replicates, group = 'scenario_id') {
  if (!is.data.frame(replicates) || !nrow(replicates) ||
      !group %in% names(replicates)) stop('Invalid aggregate input.')
  ids <- unique(replicates[[group]])
  x <- do.call(rbind, lapply(ids, function(id) {
    subset <- replicates[replicates[[group]] == id, , drop = FALSE]
    ok <- subset[subset$success, , drop = FALSE]
    total <- if (nrow(ok)) sum(ok$n_test) else 0L
    acc <- function(key) if (total) sum(ok[[key]]) / total else NA_real_
    rate <- function(key) if (nrow(ok)) mean(ok[[key]]) else NA_real_
    max_safe <- function(key) if (nrow(ok)) max(ok[[key]]) else NA_real_
    data.frame(group_id = id, jobs = nrow(subset), completed = nrow(ok),
      failures = nrow(subset) - nrow(ok), n_test = total,
      mean_realized_censoring = rate('realized_censoring'),
      original_oracle_agreement = acc('oracle_original_correct'),
      raw_oracle_agreement = acc('oracle_raw_correct'),
      rounded_oracle_agreement = acc('oracle_rounded_correct'),
      perturbed_oracle_agreement = acc('oracle_perturbed_correct'),
      raw_round_switches = if(nrow(ok)) sum(ok$raw_rounded_choice_switches) else 0L,
      raw_epsilon_switches = if(nrow(ok)) sum(ok$raw_perturbed_choice_switches) else 0L,
      max_round_survival_change = max_safe('raw_rounded_survival_max'),
      max_epsilon_survival_change = max_safe('raw_perturbed_survival_max'),
      mean_epsilon_survival_change = rate('raw_perturbed_survival_mean'),
      max_non_survival_difference = max_safe('non_survival_max'),
      all_rank_parity = if(nrow(ok)) all(ok$rank_parity) else NA,
      stringsAsFactors = FALSE)
  }))
  names(x)[1L] <- group
  rownames(x) <- NULL
  x
}

# Compare matching seeds/parameters for the two DGPs. A pair with any
# failure is excluded, not silently treated as a zero-difference result.
.phase4b_paired_models <- function(replicates) {
  a <- replicates[replicates$success &
                    replicates$survival_model == 'lognormal', , drop = FALSE]
  b <- replicates[replicates$success &
                    replicates$survival_model == 'weibull', , drop = FALSE]
  keys <- c('pair_id', 'replicate', 'seed')
  cols <- c(keys, 'realized_censoring', 'raw_perturbed_survival_max',
            'raw_rounded_survival_max', 'raw_perturbed_choice_switches',
            'oracle_raw_correct')
  x <- merge(a[, cols, drop = FALSE], b[, cols, drop = FALSE],
             by = keys, suffixes = c('_lognormal', '_weibull'))
  x$epsilon_sensitivity_delta_weibull_minus_lognormal <-
    x$raw_perturbed_survival_max_weibull -
    x$raw_perturbed_survival_max_lognormal
  x$rounding_sensitivity_delta_weibull_minus_lognormal <-
    x$raw_rounded_survival_max_weibull -
    x$raw_rounded_survival_max_lognormal
  x
}

#' Run a paired survival-model / censoring / sample-size sensitivity pilot.
#'
#' Requires original Code_set, ranger, randomForestSRC, survival, and
#' devtools::load_all(). The oracle is rebuilt for each survival family.
#' No participant-level results are written or returned.
#'
#' @param reference_directory The extracted Code_set directory or its
#'   00_utils_validation.R path.
#' @param scenarios Output of make_phase4b_survival_scenarios().
#' @param B Replicates per cell (1..10); start with B=1 or B=2.
#' @param seed_start Base random seed, matched across survival mechanisms.
#' @param digits Rounding digits for sensitivity arm, default 11.
#' @param calibration_n Independent calibration sample size for censoring.
#' @param epsilon_multiple Size of deterministic perturbation, in units of
#'   machine epsilon times max(1,abs(phi)); default 2.
#' @param mode Phase3 mode ('pilot' or 'paper_one_replicate').
#' @param continue_on_error Save failed job metadata and continue.
#' @return Metrics-only per_replicate/by_scenario/by_model/paired_models,
#'   failures, and settings (no raw simulated data or fitted models).
run_phase4b_survival_sensitivity <- function(
    reference_directory,
    scenarios = make_phase4b_survival_scenarios('pilot'),
    B = 1L, seed_start = 2026L, digits = 11L,
    calibration_n = 5000L, epsilon_multiple = 2,
    mode = c('pilot', 'paper_one_replicate'),
    continue_on_error = TRUE) {
  .phase4b_require()
  mode <- match.arg(mode)
  if (!requireNamespace('randomForestSRC', quietly = TRUE) ||
      !requireNamespace('ranger', quietly = TRUE) ||
      !requireNamespace('survival', quietly = TRUE)) {
    stop('Install randomForestSRC, ranger, and survival to run Phase 4B.')
  }
  if (length(digits) != 1L || !is.numeric(digits) || !is.finite(digits) ||
      digits != as.integer(digits) || digits < 0 || digits > 15) {
    stop('digits must be an integer between 0 and 15.')
  }
  if (!is.logical(continue_on_error) || length(continue_on_error) != 1L ||
      is.na(continue_on_error)) {
    stop('continue_on_error must be TRUE or FALSE.')
  }
  if (!is.numeric(calibration_n) || length(calibration_n) != 1L ||
      !is.finite(calibration_n) || calibration_n != as.integer(calibration_n) ||
      calibration_n < 200L || calibration_n > 50000L) {
    stop('calibration_n must be an integer from 200 to 50000.')
  }
  if (length(epsilon_multiple) != 1L || !is.numeric(epsilon_multiple) ||
      !is.finite(epsilon_multiple) || epsilon_multiple <= 0) {
    stop('epsilon_multiple must be positive.')
  }
  if (length(reference_directory) != 1L ||
      !is.character(reference_directory) || is.na(reference_directory) ||
      !(file.exists(reference_directory) || dir.exists(reference_directory))) {
    stop('Choose an existing original Code_set directory or source file.')
  }
  jobs <- .phase4b_jobs(scenarios, B, seed_start)
  rows <- vector('list', nrow(jobs))
  cat('\n=== Phase 4B: paired survival DGP / machine-precision sensitivity ===\n')
  cat('Jobs: ', nrow(jobs), '; ', mode, ' mode, ', digits,
      '-decimal rounding, epsilon multiple ', epsilon_multiple,
      ', independent censor calibration n=', calibration_n, '.\n', sep = '')
  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, , drop = FALSE]
    cat('\n=== Job ', i, '/', nrow(jobs), ' | ', job$scenario_id,
        ' | replicate ', job$replicate, ' | seed ', job$seed, ' ===\n', sep = '')
    errmsg <- ''
    outcome <- tryCatch(
      .phase4b_one_job(job, reference_directory, digits,
                       calibration_n, epsilon_multiple, mode),
      error = function(e) { errmsg <<- conditionMessage(e); NULL })
    if (is.null(outcome)) {
      cat('FAILED: ', errmsg, '\n', sep = '')
      if (!continue_on_error) stop(errmsg, call. = FALSE)
      outcome <- data.frame(n_test=NA_integer_,
        realized_censoring=NA_real_, target_censoring=NA_real_,
        oracle_original_correct=NA_integer_, oracle_raw_correct=NA_integer_,
        oracle_rounded_correct=NA_integer_, oracle_perturbed_correct=NA_integer_,
        raw_rounded_choice_switches=NA_integer_,
        raw_perturbed_choice_switches=NA_integer_,
        raw_rounded_full_rank_switches=NA_integer_,
        raw_perturbed_full_rank_switches=NA_integer_,
        raw_rounded_survival_max=NA_real_, raw_rounded_survival_mean=NA_real_,
        raw_perturbed_survival_max=NA_real_,
        raw_perturbed_survival_mean=NA_real_, non_survival_max=NA_real_,
        rank_parity=NA, stringsAsFactors = FALSE)
    } else {
      cat('Realized censoring: ', round(outcome$realized_censoring, 3),
          '; epsilon max survival change: ',
          signif(outcome$raw_perturbed_survival_max, 4),
          '; rounding max survival change: ',
          signif(outcome$raw_rounded_survival_max, 4),
          '; raw/epsilon first-choice switches: ',
          outcome$raw_perturbed_choice_switches, '\n', sep = '')
    }
    rows[[i]] <- data.frame(job, success = !nzchar(errmsg),
                            error_message = errmsg, outcome,
                            stringsAsFactors = FALSE)
    invisible(gc(verbose = FALSE))
  }
  per_rep <- do.call(rbind, rows)
  rownames(per_rep) <- NULL
  by_scenario <- .phase4b_aggregate(per_rep, 'scenario_id')
  by_model <- .phase4b_aggregate(per_rep, 'survival_model')
  pairs <- .phase4b_paired_models(per_rep)
  failures <- per_rep[!per_rep$success,
                      c('scenario_id','replicate','seed','error_message'),
                      drop = FALSE]
  ans <- list(per_replicate = per_rep, by_scenario = by_scenario,
              by_model = by_model, paired_models = pairs,
              failures = failures,
              settings = list(B=B, digits=digits, mode=mode,
                 calibration_n=calibration_n, epsilon_multiple=epsilon_multiple))
  cat('\n=== Phase 4B by survival model (descriptive) ===\n')
  print(by_model, row.names = FALSE, digits = 5)
  cat('\n=== Matched lognormal / Weibull runs ===\n')
  print(pairs, row.names = FALSE, digits = 5)
  cat('\nLIMITATIONS: pilot only; dependent held-out profiles within each ',
      'replicate; the alternative is Weibull with the original Gaussian copula, ',
      'not the paper\'s lognormal DGP; censoring calibrated using an ',
      'independent synthetic sample; no policy-value inference.\n', sep = '')
  invisible(ans)
}
