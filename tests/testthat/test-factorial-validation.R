# Load the diagnostic functions from the source tree or installed package.
phase4_test_env <- new.env(parent = globalenv())
.phase4_test_runner <- testthat::test_path('..', '..', 'inst', 'validation',
                                         'factorial_scenarios.R')
if (!file.exists(.phase4_test_runner)) {
  .phase4_test_runner <- system.file('validation',
                                     'factorial_scenarios.R',
                                     package = 'drMultiOutITR')
}
if (!nzchar(.phase4_test_runner) || !file.exists(.phase4_test_runner)) {
  stop('Phase 4 installed validation runner was not packaged.')
}
sys.source(.phase4_test_runner, envir = phase4_test_env)
rm(.phase4_test_runner)

test_that('Factorial scenario scenario definitions mirror original file-30 factorial', {
  pilot <- phase4_test_env$make_phase4_scenarios('pilot')
  full <- phase4_test_env$make_phase4_scenarios('factorial')
  expect_equal(nrow(pilot), 4L)
  expect_equal(nrow(full), 24L)
  expect_identical(pilot$scenario_id,
    c('N300_s1_rho03_c25', 'N300_s1_rho03_c50',
      'N600_s1_rho03_c25', 'N600_s1_rho03_c50'))
  expect_setequal(unique(full$n_total), c(300L, 600L))
  expect_setequal(unique(full$sigma), c(1, 2, 3))
  expect_setequal(unique(full$rho), c(.3, .8))
  expect_setequal(unique(full$censoring), c('25', '50'))
  expect_identical(full, full[!duplicated(full$scenario_id), , drop = FALSE])
})

test_that('Factorial scenario job seeds and validation match original schedule', {
  pilot <- phase4_test_env$make_phase4_scenarios('pilot')
  jobs <- phase4_test_env$.phase4_make_jobs(pilot, 2L, 2026L)
  expect_equal(nrow(jobs), 8L)
  expect_identical(jobs$seed, as.integer(c(
    2026, 2027, 102026, 102027, 202026, 202027, 302026, 302027)))
  expect_error(phase4_test_env$.phase4_make_jobs(pilot, 0L, 2026L),
               'B must be')
  expect_error(phase4_test_env$.phase4_make_jobs(pilot, 21L, 2026L),
               'B must be')
  bad <- pilot
  bad$scenario_id[2L] <- bad$scenario_id[1L]
  expect_error(phase4_test_env$.phase4_make_jobs(bad, 2L, 2026L),
               'unique')
  bad <- pilot
  bad$censoring[1L] <- '75'
  expect_error(phase4_test_env$.phase4_make_jobs(bad, 2L, 2026L),
               'constraints')
})

test_that('oracle comparisons align unique held-out IDs, not row positions', {
  oracle <- data.frame(new_id = c('x', 'y', 'z'),
    best_treatment = c('A', 'B', 'C'),
    optimal_ranking = c('A>B>C', 'B>A>C', 'C>B>A'))
  pkg <- oracle[c(3L, 1L, 2L), , drop = FALSE]
  pkg$best_treatment[pkg$new_id == 'y'] <- 'C'
  result <- phase4_test_env$.phase4_aligned_choices(oracle, pkg)
  expect_identical(result$same_best, c(TRUE, FALSE, TRUE))
  expect_true(all(result$same_rank))
  expect_equal(result$n, 3L)
  pkg$new_id[2L] <- 'z'
  expect_error(phase4_test_env$.phase4_aligned_choices(oracle, pkg),
               'unique patient IDs')
})

test_that('Factorial scenario aggregates are keyed per scenario and exclude failures', {
  dat <- data.frame(
    scenario_id = c('a', 'a', 'b'),
    success = c(TRUE, FALSE, TRUE), n_test = c(12L, NA, 12L),
    raw_oracle_correct = c(8L, NA, 9L),
    rounded_oracle_correct = c(9L, NA, 8L),
    original_oracle_correct = c(10L, NA, 10L),
    raw_oracle_rank_correct = c(7L, NA, 8L),
    rounded_oracle_rank_correct = c(8L, NA, 7L),
    original_oracle_rank_correct = c(9L, NA, 9L),
    raw_rounded_switches = c(1L, NA, 1L),
    raw_rounded_rank_switches = c(2L, NA, 3L),
    rounded_improvements = c(1L, NA, 0L),
    rounded_deteriorations = c(0L, NA, 1L),
    raw_round_survival_test_max = c(.05, NA, .07),
    raw_round_survival_scaling_max = c(.06, NA, .08),
    non_survival_max_diff = c(0, NA, 1e-12),
    rank_parity = c(TRUE, NA, TRUE))
  result <- phase4_test_env$.phase4_aggregate(dat)
  expect_equal(nrow(result), 2L)
  expect_equal(result$failed, c(1L, 0L))
  expect_equal(result$heldout_profiles, c(12L, 12L))
  expect_equal(result$raw_oracle_agreement, c(8 / 12, 9 / 12))
  expect_equal(result$rounded_oracle_agreement, c(9 / 12, 8 / 12))
  expect_equal(result$max_raw_rounded_survival_scaling_difference, c(.06, .08))
  expect_true(all(result$all_identical_input_rank_checks_pass))
  expect_equal(result$raw_rounded_first_agreement, c(11/12, 11/12))
  expect_equal(result$raw_oracle_full_rank_agreement, c(7/12, 8/12))
})
