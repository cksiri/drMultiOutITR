# Load the diagnostic functions from the source tree or installed package.
phase4b_test_env <- new.env(parent = globalenv())
.phase4b_test_runner <- testthat::test_path('..', '..', 'inst', 'validation',
                                           'survival_dgp_sensitivity.R')
if (!file.exists(.phase4b_test_runner)) {
  .phase4b_test_runner <- system.file('validation',
    'survival_dgp_sensitivity.R', package = 'drMultiOutITR')
}
if (!nzchar(.phase4b_test_runner) || !file.exists(.phase4b_test_runner)) {
  stop('Phase 4B installed validation runner is missing from package.')
}
sys.source(.phase4b_test_runner, envir = phase4b_test_env)
rm(.phase4b_test_runner)

test_that('Survival sensitivity pilot and extended grids contain complete paired designs', {
  pilot <- phase4b_test_env$make_phase4b_survival_scenarios('pilot')
  full <- phase4b_test_env$make_phase4b_survival_scenarios('extended')
  expect_equal(nrow(pilot), 4L)
  expect_equal(nrow(full), 12L)
  expect_setequal(unique(pilot$survival_model), c('lognormal', 'weibull'))
  expect_setequal(unique(pilot$censoring), c('25', '50'))
  expect_setequal(unique(full$censoring), c('none', '25', '50'))
  expect_setequal(unique(full$n_total), c(300L, 600L))
  expect_equal(as.integer(table(pilot$pair_id)), rep(2L, 2))
  expect_equal(as.integer(table(full$pair_id)), rep(2L, 6))
  expect_equal(nrow(phase4b_test_env$.phase4b_check_scenarios(pilot)), 4L)
  expect_error(phase4b_test_env$make_phase4b_survival_scenarios('pilot', 0),
               'positive')
})

test_that('Survival sensitivity seeds are matched across DGP families and differ by replicate', {
  pilot <- phase4b_test_env$make_phase4b_survival_scenarios('pilot')
  jobs <- phase4b_test_env$.phase4b_jobs(pilot, 2L, 2026L)
  expect_equal(nrow(jobs), 8L)
  for (pair in unique(jobs$pair_id)) {
    sub <- jobs[jobs$pair_id == pair, , drop = FALSE]
    expect_equal(length(unique(sub$seed)), 2L)
    expect_identical(sub$seed[1:2], sub$seed[3:4])
  }
  expect_error(phase4b_test_env$.phase4b_jobs(pilot, 0L, 2026L), 'B must')
  expect_error(phase4b_test_env$.phase4b_jobs(pilot, 11L, 2026L), 'B must')
  bad <- pilot
  bad$scenario_id[2] <- bad$scenario_id[1]
  expect_error(phase4b_test_env$.phase4b_jobs(bad, 1L, 2026L), 'Invalid')
  bad <- pilot
  bad$n_total[2] <- 600L
  expect_error(phase4b_test_env$.phase4b_jobs(bad, 1L, 2026L), 'Paired')
})

test_that('independent censoring rate calibration meets target expectation', {
  event <- exp(seq(log(.2), log(4), length.out = 2000))
  r25 <- phase4b_test_env$.phase4b_censor_rate(event, .25)
  r50 <- phase4b_test_env$.phase4b_censor_rate(event, .5)
  expect_gt(r25, 0)
  expect_gt(r50, r25)
  expect_equal(mean(-expm1(-r25 * event)), .25, tolerance = 1e-9)
  expect_equal(mean(-expm1(-r50 * event)), .5, tolerance = 1e-9)
  expect_error(phase4b_test_env$.phase4b_censor_rate(event, 1), 'target')
  expect_error(phase4b_test_env$.phase4b_censor_rate(c(event, -1), .25),
               'T_event')
})

test_that('Weibull oracle uses conditional Weibull survival probabilities', {
  ref <- new.env(parent = globalenv())
  ref$smmr_f_mean <- function(j, k, x) {
    stopifnot(j == 3L)
    2 * as.numeric(k) + x$X1
  }
  xx <- data.frame(X1 = c(-.5, .8))
  times <- c(.7, 1.1, 2)
  long <- phase4b_test_env$.phase4b_true_weibull_survival(
    ref, xx, times, d2 = 5, treat_levels = c('1', '2', '3'), shape = 6)
  expect_equal(nrow(long), 18L)
  expect_equal(length(unique(long$new_id)), 2L)
  expect_equal(length(unique(long$treatment)), 3L)
  expect_true(all(long$theta_hat > 0 & long$theta_hat <= 1))
  row <- long[long$new_id == 1 & long$time_index == 2 &
                long$treatment == '2', ]
  expected <- exp(- (1.1 / exp((4 - .5) / 5))^6)
  expect_equal(row$theta_hat, expected, tolerance = 1e-13)
  for (i in unique(long$new_id)) for (k in unique(long$treatment)) {
    vals <- long$theta_hat[long$new_id == i & long$treatment == k]
    expect_true(all(diff(vals) <= 0))
  }
})

test_that('the Weibull quantile transform preserves the Gaussian event score', {
  ref <- new.env(parent = globalenv())
  ref$smmr_f_mean <- function(j, k, x) rep(1.2, nrow(x))
  z <- c(-1.4, 0, 1.3, -.2)
  sigma <- 1.2
  d1 <- 10
  d2 <- 5
  sigma_L <- .2
  sd_log <- sqrt((sigma / (d1 * d2))^2 + sigma_L^2)
  T_original <- exp(1.2 / d2 + z * sd_log)
  sim <- list(dat = data.frame(A = factor(c('1','2','3','1')),
                               X1 = rep(.1, 4)),
              T_event = T_original, sigma = sigma,
              d1 = d1, d2 = d2, sigma_L = sigma_L)
  T_weib <- phase4b_test_env$.phase4b_weibull_times(ref, sim, 6)
  expected <- exp(1.2 / d2) * (-log(stats::pnorm(-z)))^(1 / 6)
  expect_equal(T_weib, expected, tolerance = 1e-13)
  expect_identical(order(T_original), order(T_weib))
  expect_error(phase4b_test_env$.phase4b_weibull_times(ref, sim, 0),
               'Invalid Weibull')
})

test_that('model summary aggregates only successful jobs', {
  z <- data.frame(scenario_id = c('a', 'a', 'b'),
    survival_model = c('lognormal','lognormal','weibull'),
    success = c(TRUE, FALSE, TRUE), n_test = c(12L, NA, 12L),
    realized_censoring = c(.25, NA, .50),
    oracle_original_correct = c(9L, NA, 7L),
    oracle_raw_correct = c(9L, NA, 8L),
    oracle_rounded_correct = c(10L, NA, 7L),
    oracle_perturbed_correct = c(8L, NA, 8L),
    raw_rounded_choice_switches = c(1L, NA, 0L),
    raw_perturbed_choice_switches = c(2L, NA, 1L),
    raw_rounded_survival_max = c(.05, NA, .1),
    raw_perturbed_survival_max = c(.02, NA, .03),
    raw_perturbed_survival_mean = c(.001, NA, .002),
    non_survival_max = c(0, NA, 1e-12),
    rank_parity = c(TRUE, NA, TRUE))
  ans <- phase4b_test_env$.phase4b_aggregate(z)
  expect_equal(ans$jobs, c(2L, 1L))
  expect_equal(ans$completed, c(1L, 1L))
  expect_equal(ans$failures, c(1L, 0L))
  expect_equal(ans$raw_oracle_agreement, c(9/12, 8/12))
  expect_equal(ans$raw_epsilon_switches, c(2L, 1L))
  expect_equal(ans$max_epsilon_survival_change, c(.02, .03))
  expect_true(all(ans$all_rank_parity))
})
