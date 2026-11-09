test_that('optional forest rounding is opt-in and rejects invalid specifications', {
  helper <- drMultiOutITR:::.dr_validate_forest_digits
  expect_null(helper(NULL, 'ranger'))
  expect_identical(helper(12L, 'ranger'), 12L)
  for (bad in list(-1L, 16L, 12.5, NaN, NA_real_, Inf, '12', numeric())) {
    expect_error(helper(bad, 'ranger'), 'integer between 0 and 15')
  }
  expect_error(helper(12L, 'lm'), 'requires second_stage="ranger"')
  phi <- c(-0.123456789012345, 0.999999999999999, 1.23456789123456)
  expect_identical(drMultiOutITR:::.dr_response_for_second_stage(phi, NULL), phi)
  expect_identical(drMultiOutITR:::.dr_response_for_second_stage(phi, 12L),
                   round(phi, 12L))
})

test_that('rounded identical forest responses produce identical refits', {
  skip_if_not_installed('ranger')
  set.seed(143L)
  x <- data.frame(z = rnorm(90), x = runif(90))
  phi <- round(sin(x$z) + x$x, 7L) + 0.00000037
  perturbed <- phi + rep(c(0, 4.440892098500626e-16), length.out = 90L)
  expect_true(max(abs(phi - perturbed)) > 0)
  y1 <- drMultiOutITR:::.dr_response_for_second_stage(phi, 12L)
  y2 <- drMultiOutITR:::.dr_response_for_second_stage(perturbed, 12L)
  expect_identical(y1, y2)
  fit1 <- drMultiOutITR:::.dr_fit_model(
    y1, x, c('x', 'z'), 'ranger', 201L, 35L, 5L)
  fit2 <- drMultiOutITR:::.dr_fit_model(
    y2, x, c('x', 'z'), 'ranger', 201L, 35L, 5L)
  expect_identical(drMultiOutITR:::.dr_predict_model(fit1, x, c('x', 'z')),
                   drMultiOutITR:::.dr_predict_model(fit2, x, c('x', 'z')))
})

test_that('survival rounding changes only the second-stage forest response', {
  skip_if_not_installed('survival')
  skip_if_not_installed('ranger')
  d <- make_dr_test_data(n = 180L)
  ps <- fit_dr_propensity(d, 'a', c('x', 'z'), method = 'multinomial',
                          folds = 1L, eps = .01)
  arguments <- list(data = d, treatment = 'a', time = 'time', status = 'status',
                    covariates = c('x', 'z'), propensity = ps,
                    times = c(.75, 1.5), outcome_learner = 'lm',
                    second_stage = 'ranger', num.trees = 25L,
                    min.node.size = 5L, seed = 850L, return_models = TRUE)
  raw <- do.call(estimate_dr_survival, arguments)
  stable <- do.call(estimate_dr_survival,
                    c(arguments, list(second_stage_round_digits = 12L)))
  expect_null(raw$second_stage_round_digits)
  expect_identical(stable$second_stage_round_digits, 12L)
  for (j in seq_along(raw$models)) {
    expect_identical(raw$models[[j]]$pseudo_outcome,
                     stable$models[[j]]$pseudo_outcome)
    expect_identical(raw$models[[j]]$first_stage_prediction,
                     stable$models[[j]]$first_stage_prediction)
    expect_null(raw$models[[j]]$pseudo_outcome_for_fit)
    expect_identical(stable$models[[j]]$pseudo_outcome_for_fit,
                     round(stable$models[[j]]$pseudo_outcome, 12L))
  }
  expect_true(all(is.finite(stable$theta_array)))
})

test_that('integrated survival-only option preserves other outcome forests', {
  skip_if_not_installed('survival')
  skip_if_not_installed('ranger')
  d <- make_dr_test_data(n = 160L)
  specs <- list(nonuse = list(type = 'binary', column = 'ybin'),
                survival = list(type = 'survival', time = 'time',
                                status = 'status', times = c(.75, 1.5)))
  args <- list(data = d, treatment = 'a', covariates = c('x', 'z'),
               outcomes = specs, propensity_method = 'multinomial',
               nuisance_folds = 1L, outcome_learner = 'ranger',
               covariance_source = 'in_sample', num.trees = 25L,
               min.node.size = 5L, seed = 624L)
  baseline <- do.call(fit_dr_itr, args)
  stable <- do.call(fit_dr_itr,
                    c(args, list(survival_second_stage_round_digits = 12L)))
  expect_identical(baseline$estimates$nonuse$theta_array,
                   stable$estimates$nonuse$theta_array)
  expect_identical(baseline$propensity$ehat, stable$propensity$ehat)
  expect_null(baseline$estimates$survival$second_stage_round_digits)
  expect_identical(stable$estimates$survival$second_stage_round_digits, 12L)
  for (j in seq_along(stable$estimates$survival$models)) {
    expect_identical(baseline$estimates$survival$models[[j]]$pseudo_outcome,
                     stable$estimates$survival$models[[j]]$pseudo_outcome)
    expect_identical(stable$estimates$survival$models[[j]]$pseudo_outcome_for_fit,
                     round(stable$estimates$survival$models[[j]]$pseudo_outcome, 12L))
  }
  expect_s3_class(stable$rank_rule, 'dr_rank_rule')
  expect_error(do.call(fit_dr_itr,
                       c(args, list(second_stage_round_digits = 16L))),
               'integer between 0 and 15')
})
