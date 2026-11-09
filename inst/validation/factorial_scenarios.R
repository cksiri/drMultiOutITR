# Factorial simulation diagnostics for original-code comparisons.
# This script requires the study's separate reference simulation utilities.
# The simulated oracle ranks are not independent policy-value estimates.

#' Return four pilot cells or the complete original 24-cell factorial grid.
#'
#' The four-cell pilot varies sample size (300/600) and censoring (25/50)
#' with sigma=1 and rho=0.3. The factorial design additionally varies
#' sigma in {1,2,3} and rho in {0.3,0.8}. The original file 30 is the
#' authoritative source for these scenario definitions.
make_phase4_scenarios <- function(design = c('pilot', 'factorial')) {
  design <- match.arg(design)
  grid <- expand.grid(
    n_total = c(300L, 600L), sigma = c(1, 2, 3),
    rho = c(0.3, 0.8), censoring = c('25', '50'),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- grid[order(grid$n_total, grid$sigma, grid$rho,
                     as.integer(grid$censoring)), , drop = FALSE]
  if (design == 'pilot') {
    grid <- grid[grid$sigma == 1 & grid$rho == .3, , drop = FALSE]
  }
  grid$scenario_id <- paste0('N', grid$n_total, '_s', grid$sigma,
                             '_rho', ifelse(grid$rho == .3, '03', '08'),
                             '_c', grid$censoring)
  grid <- grid[, c('scenario_id', 'n_total', 'sigma', 'rho', 'censoring')]
  rownames(grid) <- NULL
  grid
}

.phase4_validate_scenarios <- function(scenarios) {
  if (!is.data.frame(scenarios) || !nrow(scenarios)) {
    stop('scenarios must be a nonempty data frame.')
  }
  required <- c('scenario_id', 'n_total', 'sigma', 'rho', 'censoring')
  if (!all(required %in% names(scenarios))) {
    stop('scenarios must contain: ', paste(required, collapse = ', '))
  }
  if (anyNA(scenarios[, required, drop = FALSE])) {
    stop('scenarios may not contain missing values.')
  }
  sid <- as.character(scenarios$scenario_id)
  if (any(!nzchar(sid)) || anyDuplicated(sid)) {
    stop('scenario_id must be nonempty and unique.')
  }
  if (!is.numeric(scenarios$n_total) ||
      any(!(scenarios$n_total %in% c(300L, 600L))) ||
      !is.numeric(scenarios$sigma) ||
      any(!is.finite(scenarios$sigma)) ||
      any(scenarios$sigma <= 0) ||
      !is.numeric(scenarios$rho) ||
      any(!is.finite(scenarios$rho)) ||
      any(scenarios$rho < -0.5 | scenarios$rho >= 1) ||
      any(!(as.character(scenarios$censoring) %in% c('25', '50')))) {
    stop('Scenario values must satisfy Phase 3 n_total/sigma/rho/censoring constraints.')
  }
  scenarios$scenario_id <- sid
  scenarios$n_total <- as.integer(scenarios$n_total)
  scenarios$sigma <- as.numeric(scenarios$sigma)
  scenarios$rho <- as.numeric(scenarios$rho)
  scenarios$censoring <- as.character(scenarios$censoring)
  scenarios <- scenarios[, required, drop = FALSE]
  rownames(scenarios) <- NULL
  scenarios
}

.phase4_make_jobs <- function(scenarios, B, seed_start) {
  scenarios <- .phase4_validate_scenarios(scenarios)
  if (length(B) != 1L || !is.numeric(B) || is.na(B) ||
      !is.finite(B) || B != as.integer(B) || B < 1L || B > 20L) {
    stop('B must be an integer between 1 and 20 for this validation pilot.')
  }
  if (length(seed_start) != 1L || !is.numeric(seed_start) ||
      is.na(seed_start) || !is.finite(seed_start) ||
      seed_start != as.integer(seed_start) || seed_start < 1) {
    stop('seed_start must be one positive integer.')
  }
  # Match file 30 seed schedule: start + (scenario_index-1)*100000 + replicate-1.
  jobs <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    data.frame(scenario_index = i, scenarios[rep(i, B), , drop = FALSE],
               replicate = seq_len(B), stringsAsFactors = FALSE)
  }))
  s <- as.double(seed_start) +
    (as.double(jobs$scenario_index) - 1) * 100000 + jobs$replicate - 1
  if (any(s > .Machine$integer.max)) {
    stop('Seed schedule exceeds R integer range. Use fewer scenarios/smaller seed_start.')
  }
  jobs$seed <- as.integer(s)
  rownames(jobs) <- NULL
  jobs
}

.phase4_aligned_choices <- function(reference, comparison, label = 'comparison') {
  required <- c('new_id', 'best_treatment', 'optimal_ranking')
  if (!is.data.frame(reference) || !is.data.frame(comparison) ||
      !all(required %in% names(reference)) ||
      !all(required %in% names(comparison))) {
    stop('Both selections need new_id, best_treatment, optimal_ranking: ', label)
  }
  a <- as.character(reference$new_id)
  b <- as.character(comparison$new_id)
  if (!length(a) || anyNA(a) || anyNA(b) || anyDuplicated(a) ||
      anyDuplicated(b) || length(a) != length(b) || !setequal(a, b)) {
    stop('Selections must have the same unique patient IDs: ', label)
  }
  comparison <- comparison[match(a, b), , drop = FALSE]
  original_best <- as.character(reference$best_treatment)
  new_best <- as.character(comparison$best_treatment)
  original_rank <- as.character(reference$optimal_ranking)
  new_rank <- as.character(comparison$optimal_ranking)
  if (anyNA(c(original_best, new_best, original_rank, new_rank))) {
    stop('Missing selection or ranking: ', label)
  }
  list(n = length(a), same_best = original_best == new_best,
       same_rank = original_rank == new_rank,
       original_best = original_best, comparison_best = new_best)
}

.phase4_aggregate <- function(replicates) {
  if (!is.data.frame(replicates) || !nrow(replicates) ||
      !all(c('scenario_id', 'success', 'n_test', 'raw_oracle_correct',
             'rounded_oracle_correct', 'original_oracle_correct',
             'raw_rounded_switches', 'raw_rounded_rank_switches',
             'rounded_improvements', 'rounded_deteriorations',
             'raw_oracle_rank_correct', 'rounded_oracle_rank_correct',
             'original_oracle_rank_correct',
             'raw_round_survival_test_max', 'raw_round_survival_scaling_max',
             'non_survival_max_diff', 'rank_parity') %in% names(replicates))) {
    stop('Missing fields in Phase 4 replicate results.')
  }
  ids <- unique(as.character(replicates$scenario_id))
  rows <- lapply(ids, function(id) {
    sub <- replicates[replicates$scenario_id == id, , drop = FALSE]
    ok <- sub[sub$success %in% TRUE, , drop = FALSE]
    nprof <- if (nrow(ok)) sum(ok$n_test) else 0L
    calc_accuracy <- function(col) {
      if (!nprof) NA_real_ else sum(ok[[col]]) / nprof
    }
    max_safe <- function(col) {
      if (!nrow(ok)) NA_real_ else max(ok[[col]])
    }
    data.frame(scenario_id = id, jobs = nrow(sub), completed = nrow(ok),
      failed = nrow(sub) - nrow(ok), heldout_profiles = nprof,
      original_oracle_agreement = calc_accuracy('original_oracle_correct'),
      raw_oracle_agreement = calc_accuracy('raw_oracle_correct'),
      rounded_oracle_agreement = calc_accuracy('rounded_oracle_correct'),
      original_oracle_full_rank_agreement = calc_accuracy('original_oracle_rank_correct'),
      raw_oracle_full_rank_agreement = calc_accuracy('raw_oracle_rank_correct'),
      rounded_oracle_full_rank_agreement = calc_accuracy('rounded_oracle_rank_correct'),
      raw_rounded_first_agreement = if (nprof) 1 - sum(ok$raw_rounded_switches)/nprof else NA_real_,
      raw_rounded_full_rank_agreement = if (nprof) 1 - sum(ok$raw_rounded_rank_switches)/nprof else NA_real_,
      raw_rounded_choice_switches = if (nrow(ok)) sum(ok$raw_rounded_switches) else 0L,
      raw_rounded_ranking_switches = if (nrow(ok)) sum(ok$raw_rounded_rank_switches) else 0L,
      rounding_improvements_vs_oracle = if (nrow(ok)) sum(ok$rounded_improvements) else 0L,
      rounding_deteriorations_vs_oracle = if (nrow(ok)) sum(ok$rounded_deteriorations) else 0L,
      max_raw_rounded_survival_test_difference = max_safe('raw_round_survival_test_max'),
      max_raw_rounded_survival_scaling_difference = max_safe('raw_round_survival_scaling_max'),
      max_non_survival_reference_difference = max_safe('non_survival_max_diff'),
      all_identical_input_rank_checks_pass = if (nrow(ok)) all(ok$rank_parity) else NA,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.phase4_one_job <- function(job, reference_directory, digits, mode) {
  one <- run_phase3_final_simulation(
    reference_directory = reference_directory,
    mode = mode, reference_survival = 'final_rsf', alignment = 'research',
    seed = job$seed, n_total = job$n_total, sigma = job$sigma,
    rho = job$rho, censoring = job$censoring)
  ref <- one$original_fit
  package_raw <- one$package_fit
  stable <- .phase3c_refit_survival_only(package_raw, ref$dat, digits)

  # Audit that the paired intervention does NOT alter the nuisance estimates
  # or the original (unrounded) AIPW training responses.
  if (!identical(package_raw$propensity$ehat, stable$propensity$ehat)) {
    stop('Propensity nuisance estimates changed during paired refit.')
  }
  if (!identical(lapply(package_raw$estimates$survival$models,
                        function(x) x$pseudo_outcome),
                 lapply(stable$estimates$survival$models,
                        function(x) x$pseudo_outcome))) {
    stop('Survival AIPW pseudo-outcomes changed during paired refit.')
  }
  if (!isTRUE(one$ranking_parity_on_identical_summaries)) {
    stop('Ranking parity on identical summary inputs failed.')
  }
  xnew <- ref$xnew
  ids <- as.character(seq_len(nrow(xnew)))
  rounded <- stats::predict(stable, xnew, ids = ids,
                            details = TRUE, return_summaries = TRUE)
  raw <- one$package_predictions
  trts <- as.character(ref$sim$treat_levels)
  compare_test <- .phase3_metric('raw vs rounded: survival heldout',
    .phase3_summary_array(raw$summaries, 'survival', ids, trts),
    .phase3_summary_array(rounded$summaries, 'survival', ids, trts))
  compare_reference <- .phase3_metric('original vs rounded: survival heldout',
    .phase3_summary_array(ref$estimated_summary_tables$test, 'survival', ids, trts),
    .phase3_summary_array(rounded$summaries, 'survival', ids, trts))
  training_ids <- as.character(seq_len(nrow(ref$dat)))
  raw_train <- drMultiOutITR:::.dr_itr_table(
    package_raw$estimates, package_raw$outcomes, ref$dat,
    package_raw$covariates, training_ids, package_raw$weights,
    use_oob = FALSE)
  rounded_train <- drMultiOutITR:::.dr_itr_table(
    stable$estimates, stable$outcomes, ref$dat,
    stable$covariates, training_ids, stable$weights, use_oob = FALSE)
  compare_scaling <- .phase3_metric('raw vs rounded: survival training',
    .phase3_summary_array(raw_train, 'survival', training_ids, trts),
    .phase3_summary_array(rounded_train, 'survival', training_ids, trts))

  oracle <- ref$matrix_fits$true$selection
  original <- ref$matrix_fits$estimated$selection
  raw_select <- raw$selection
  round_select <- rounded$selection
  oracle_raw <- .phase4_aligned_choices(oracle, raw_select, 'oracle/raw')
  oracle_round <- .phase4_aligned_choices(oracle, round_select, 'oracle/rounded')
  oracle_orig <- .phase4_aligned_choices(oracle, original, 'oracle/original')
  raw_round <- .phase4_aligned_choices(raw_select, round_select, 'raw/rounded')
  orig_raw <- .phase4_aligned_choices(original, raw_select, 'original/raw')
  orig_round <- .phase4_aligned_choices(original, round_select, 'original/rounded')
  if (any(!is.finite(unlist(one$comparison$max_abs_diff))) ||
      any(!is.finite(unlist(one$scaling_comparison$max_abs_diff)))) {
    stop('Nonfinite numerical comparison in Phase 3 result.')
  }
  non_survival <- c(one$comparison$max_abs_diff[
    one$comparison$component %in% c('propensity (estimated)',
    'counterfactual summary: continuous', 'counterfactual summary: ordinal')],
    one$scaling_comparison$max_abs_diff[
      one$scaling_comparison$component %in%
        c('scaling summary: continuous', 'scaling summary: ordinal')])
  if (length(non_survival) != 5L) {
    stop('Expected five non-survival reference parity metrics.')
  }
  n <- oracle_raw$n
  if (n != nrow(xnew) || oracle_round$n != n ||
      oracle_orig$n != n || raw_round$n != n) {
    stop('Selection comparison has unexpected number of profiles.')
  }
  # Restore all vectors to oracle ID ordering for improvement counts.
  raw_matches <- oracle_raw$same_best
  rounded_matches <- oracle_round$same_best
  data.frame(
    n_test = n, realized_censoring = mean(ref$dat$status == 0),
    raw_oracle_correct = sum(raw_matches),
    rounded_oracle_correct = sum(rounded_matches),
    original_oracle_correct = sum(oracle_orig$same_best),
    raw_oracle_rank_correct = sum(oracle_raw$same_rank),
    rounded_oracle_rank_correct = sum(oracle_round$same_rank),
    original_oracle_rank_correct = sum(oracle_orig$same_rank),
    raw_rounded_switches = sum(!raw_round$same_best),
    raw_rounded_rank_switches = sum(!raw_round$same_rank),
    rounded_improvements = sum(!raw_matches & rounded_matches),
    rounded_deteriorations = sum(raw_matches & !rounded_matches),
    original_raw_first_agreement = mean(orig_raw$same_best),
    original_rounded_first_agreement = mean(orig_round$same_best),
    raw_round_survival_test_max = compare_test$max_abs_diff,
    raw_round_survival_test_mean = compare_test$mean_abs_diff,
    raw_round_survival_scaling_max = compare_scaling$max_abs_diff,
    original_round_survival_test_max = compare_reference$max_abs_diff,
    non_survival_max_diff = max(non_survival),
    rank_parity = one$ranking_parity_on_identical_summaries,
    stringsAsFactors = FALSE)
}

#' Run paired raw/11-decimal experiments across original simulation scenarios.
#'
#' Each job uses the same generated data, first-stage models, and AIPW
#' pseudo-outcomes in its raw and rounded PACKAGE fits. Both are compared
#' to the original study's Monte-Carlo oracle rule. Runs sequentially.
#'
#' @param reference_directory Original Code_set folder or 00_utils_validation.R.
#' @param scenarios Output of make_phase4_scenarios() or subset of its rows.
#' @param B Replicates per scenario, 1-20. Start with B=2.
#' @param seed_start First simulation seed, default 2026; matches original
#'   file-30 seed indexing over scenarios passed to this function.
#' @param digits Optional response rounding precision, default 11L. The
#'   unrounded estimator remains the package default.
#' @param mode 'pilot': 120 trees, 12 new profiles, 120 oracle scaling profiles;
#'   'paper_one_replicate': 600 trees, 100 new profiles, 1000 scaling profiles.
#' @param continue_on_error If TRUE, capture per-job errors and proceed.
#' @param output_dir Optional existing folder in which to write *metrics only*
#'   CSVs; NULL writes nothing.
#' @return A list with jobs, per_replicate, per_scenario, overall,
#'   failures, settings; does not retain data or models.
run_phase4_scenario_validation <- function(
    reference_directory,
    scenarios = make_phase4_scenarios('pilot'), B = 2L,
    seed_start = 2026L, digits = 11L, mode = c('pilot', 'paper_one_replicate'),
    continue_on_error = TRUE, output_dir = NULL) {
  mode <- match.arg(mode)
  # Require the companion research comparators with compatible output fields.
  required <- c('run_phase3_final_simulation', '.phase3_metric',
                '.phase3_summary_array', '.phase3c_refit_survival_only')
  missing <- required[!vapply(required, exists, logical(1), mode = 'function')]
  if (length(missing)) {
    stop('First source validation/run_final_simulation_phase3.R and ',
         'validation/run_phase3c_multi_replicate.R; missing: ',
         paste(missing, collapse = ', '))
  }
  if (!exists('fit_dr_itr', mode = 'function')) {
    stop('First run devtools::load_all() on drMultiOutITR.')
  }
  if (!is.logical(continue_on_error) || length(continue_on_error) != 1L ||
      is.na(continue_on_error)) {
    stop('continue_on_error must be TRUE or FALSE.')
  }
  if (length(digits) != 1L || !is.numeric(digits) || is.na(digits) ||
      !is.finite(digits) || digits != as.integer(digits) ||
      digits < 0 || digits > 15) {
    stop('digits must be an integer from 0 to 15.')
  }
  if (!is.null(output_dir) &&
      (length(output_dir) != 1L || !is.character(output_dir) ||
       is.na(output_dir) || !dir.exists(output_dir))) {
    stop('output_dir must be NULL or an existing folder.')
  }
  jobs <- .phase4_make_jobs(scenarios, B, seed_start)
  if (length(reference_directory) != 1L ||
      !is.character(reference_directory) || is.na(reference_directory) ||
      !(dir.exists(reference_directory) || file.exists(reference_directory))) {
    stop('reference_directory must be an existing Code_set directory or file.')
  }
  if (!requireNamespace('randomForestSRC', quietly = TRUE)) {
    stop('Install randomForestSRC for the manuscript-style survival learner.')
  }
  results <- vector('list', nrow(jobs))
  cat('\n=== Phase 4: scenario-based ORACLE agreement & rounding sensitivity ===\n')
  cat('Jobs: ', nrow(jobs), '; mode: ', mode, '; rounding: ', digits,
      ' decimals; research configuration; equal weights.\n', sep = '')
  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, , drop = FALSE]
    cat('\n=== Job ', i, '/', nrow(jobs), ' | ', job$scenario_id,
        ' | replicate ', job$replicate, ' | seed ', job$seed,
        ' ===\n', sep = '')
    error_text <- ''
    ans <- tryCatch(
      .phase4_one_job(job, reference_directory, digits, mode),
      error = function(e) {
        error_text <<- conditionMessage(e)
        NULL
      })
    if (is.null(ans)) {
      cat('FAILED: ', error_text, '\n', sep = '')
      if (!continue_on_error) stop(error_text, call. = FALSE)
      results[[i]] <- data.frame(
        job, success = FALSE, error_message = error_text,
        n_test = NA_integer_, realized_censoring = NA_real_,
        raw_oracle_correct = NA_integer_,
        rounded_oracle_correct = NA_integer_,
        original_oracle_correct = NA_integer_,
        raw_oracle_rank_correct = NA_integer_,
        rounded_oracle_rank_correct = NA_integer_,
        original_oracle_rank_correct = NA_integer_,
        raw_rounded_switches = NA_integer_,
        raw_rounded_rank_switches = NA_integer_,
        rounded_improvements = NA_integer_,
        rounded_deteriorations = NA_integer_,
        original_raw_first_agreement = NA_real_,
        original_rounded_first_agreement = NA_real_,
        raw_round_survival_test_max = NA_real_,
        raw_round_survival_test_mean = NA_real_,
        raw_round_survival_scaling_max = NA_real_,
        original_round_survival_test_max = NA_real_,
        non_survival_max_diff = NA_real_,
        rank_parity = NA, stringsAsFactors = FALSE)
    } else {
      results[[i]] <- data.frame(job, success = TRUE,
                                  error_message = '', ans,
                                  stringsAsFactors = FALSE)
      cat('Oracle concordance (original/raw/rounded): ',
          ans$original_oracle_correct, '/', ans$n_test, ' / ',
          ans$raw_oracle_correct, '/', ans$n_test, ' / ',
          ans$rounded_oracle_correct, '/', ans$n_test, '\n', sep = '')
      cat('Raw vs rounded first-choice switches: ',
          ans$raw_rounded_switches,
          '; maximum survival probability change: ',
          format(ans$raw_round_survival_test_max, digits = 5), '\n', sep = '')
      if (ans$non_survival_max_diff > 1e-8) {
        cat('CHECK: non-survival reference differences remain: ',
            format(ans$non_survival_max_diff, digits = 5), '\n', sep = '')
      }
    }
    invisible(gc(verbose = FALSE))
  }
  per_rep <- do.call(rbind, results)
  rownames(per_rep) <- NULL
  per_scenario <- .phase4_aggregate(per_rep)
  ok <- per_rep[per_rep$success, , drop = FALSE]
  denom <- if (nrow(ok)) sum(ok$n_test) else 0L
  rate <- function(col) if (denom) sum(ok[[col]])/denom else NA_real_
  overall <- data.frame(
    scenario_cells = nrow(scenarios), jobs = nrow(jobs),
    completed = nrow(ok), failed = nrow(jobs) - nrow(ok),
    heldout_profiles = denom,
    original_oracle_agreement = rate('original_oracle_correct'),
    raw_oracle_agreement = rate('raw_oracle_correct'),
    rounded_oracle_agreement = rate('rounded_oracle_correct'),
    original_oracle_full_rank_agreement = rate('original_oracle_rank_correct'),
    raw_oracle_full_rank_agreement = rate('raw_oracle_rank_correct'),
    rounded_oracle_full_rank_agreement = rate('rounded_oracle_rank_correct'),
    raw_rounded_first_agreement = if (denom) 1 - sum(ok$raw_rounded_switches)/denom else NA_real_,
    raw_rounded_full_rank_agreement = if (denom) 1 - sum(ok$raw_rounded_rank_switches)/denom else NA_real_,
    raw_rounded_choice_switches = if (nrow(ok)) sum(ok$raw_rounded_switches) else 0L,
    raw_rounded_ranking_switches = if (nrow(ok)) sum(ok$raw_rounded_rank_switches) else 0L,
    rounding_improvements_vs_oracle = if (nrow(ok)) sum(ok$rounded_improvements) else 0L,
    rounding_deteriorations_vs_oracle = if (nrow(ok)) sum(ok$rounded_deteriorations) else 0L,
    max_raw_rounded_survival_test_difference = if (nrow(ok)) max(ok$raw_round_survival_test_max) else NA_real_,
    max_raw_rounded_survival_scaling_difference = if (nrow(ok)) max(ok$raw_round_survival_scaling_max) else NA_real_,
    max_non_survival_reference_difference = if (nrow(ok)) max(ok$non_survival_max_diff) else NA_real_,
    all_identical_input_rank_checks_pass = if (nrow(ok)) all(ok$rank_parity) else NA,
    stringsAsFactors = FALSE)
  failures <- per_rep[!per_rep$success,
                      c('scenario_id', 'replicate', 'seed', 'error_message'),
                      drop = FALSE]
  result <- list(
    per_replicate = per_rep, per_scenario = per_scenario,
    overall = overall, failures = failures, jobs = jobs,
    settings = list(scenarios = scenarios, B = B, seed_start = seed_start,
                    digits = digits, mode = mode,
                    research_alignment = 'final_rsf / research / equal_weights',
                    reference = 'original Code_set; estimated vs Monte Carlo DGP oracle',
                    package_default_rounding = NULL),
    versions = list(R = as.character(getRversion()),
                    ranger = as.character(utils::packageVersion('ranger')),
                    randomForestSRC = as.character(utils::packageVersion('randomForestSRC'))))
  cat('\n=== Phase 4 by-scenario results ===\n')
  print(per_scenario, row.names = FALSE, digits = 5)
  cat('\n=== Phase 4 overall pilot (descriptive, NOT independent-patient inference) ===\n')
  print(overall, row.names = FALSE, digits = 5)
  if (nrow(failures)) {
    cat('\nFailed jobs (inspect before interpreting completed jobs):\n')
    print(failures, row.names = FALSE)
  }
  cat('\nNOTE: held-out profiles within a fitted replicate are dependent. ',
      'Do not compute patient-level standard errors from these totals.\n', sep = '')
  cat('NOTE: This is a small, equal-weight pilot using a Monte Carlo DGP ',
      'oracle, not manuscript B=500 or estimated policy value.\n', sep = '')
  if (!is.null(output_dir)) {
    utils::write.csv(per_rep, file.path(output_dir,
                                       'phase4_per_replicate.csv'), row.names = FALSE)
    utils::write.csv(per_scenario, file.path(output_dir,
                                           'phase4_per_scenario.csv'), row.names = FALSE)
    utils::write.csv(overall, file.path(output_dir,
                                      'phase4_overall.csv'), row.names = FALSE)
    utils::write.csv(failures, file.path(output_dir,
                                       'phase4_failures.csv'), row.names = FALSE)
    cat('Wrote aggregate metrics CSVs to: ', output_dir, '\n', sep = '')
  }
  invisible(result)
}
