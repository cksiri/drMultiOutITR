test_that('research propensity option matches original floor-then-renormalize rule', {
  skip_if_not_installed('ranger')
  d <- make_dr_test_data(n = 240L)
  d$a <- factor(rep(c('LTR', 'STI', 'ODF'), length.out = nrow(d)),
                levels = c('LTR', 'STI', 'ODF'))
  ps <- fit_dr_propensity(d, 'a', c('x', 'z'), folds = 1L,
                          eps = .05, num.trees = 40L, min.node.size = 12L,
                          seed = 702L, truncation = 'floor_renorm',
                          seed_schedule = 'research')
  df <- data.frame(.treatment = d$a, x = d$x, z = d$z)
  set.seed(702L)
  fit <- ranger::ranger(.treatment ~ x + z, data = df, probability = TRUE,
                        num.trees = 40L, min.node.size = 12L,
                        respect.unordered.factors = 'order', seed = 702L)
  raw <- stats::predict(fit, data = d[, c('x', 'z')])$predictions
  expected <- pmax(raw[, levels(d$a), drop = FALSE], .05)
  expected <- expected / rowSums(expected)
  expect_equal(unname(ps$ehat), unname(expected), tolerance = 1e-12)
  expect_equal(rowSums(ps$ehat), rep(1, nrow(d)), tolerance = 1e-12)
  expect_identical(ps$truncation, 'floor_renorm')
  expect_identical(ps$seed_schedule, 'research')
})

test_that('research seed schedule gives file-02 and file-03 second-stage offsets', {
  skip_if_not_installed('ranger')
  d <- make_dr_test_data(n = 240L)
  ps <- fit_dr_propensity(d, 'a', c('x', 'z'), method = 'multinomial',
                          folds = 1L, eps = .01)
  ord <- estimate_dr_ordinal(d, 'a', 'yord', c('x', 'z'), ps,
                             seed = 300L, seed_schedule = 'research',
                             num.trees = 25L, min.node.size = 8L,
                             return_models = TRUE)
  cont <- estimate_dr_continuous(d, 'a', 'ycont', c('x', 'z'), ps,
                                 probs = c(.5), grid_length = 5L,
                                 seed = 400L, seed_schedule = 'research',
                                 num.trees = 25L, min.node.size = 8L,
                                 return_models = TRUE)
  expect_equal(ord$models[[1L]]$second_seed, 300L + 2000L + 1L)
  expect_equal(ord$models[[3L]]$second_seed, 300L + 4000L + 1L)
  expect_equal(cont$models[[1L]]$second_seed, 400L + 2000L + 1L)
  expect_true(all(is.finite(cont$theta_array)))
})

test_that('research RSF fits event-survival nuisances and the correct AIPW formula', {
  skip_if_not_installed('survival')
  skip_if_not_installed('randomForestSRC')
  d <- make_dr_test_data(n = 180L)
  ps <- fit_dr_propensity(d, 'a', c('x','z'), method = 'multinomial',
                          folds = 1L, eps = .01)
  fit <- estimate_dr_survival(d, 'a', 'time', 'status', c('x', 'z'), ps,
                              times = c(.75, 1.5), outcome_learner = 'rsf',
                              second_stage = 'ranger', seed_schedule = 'research',
                              seed = 901L, num.trees = 25L,
                              min.node.size = 5L, rsf_num.trees = 25L,
                              rsf_nodesize = 5L, return_models = TRUE)
  expect_equal(fit$outcome_learner, 'rsf')
  expect_true(all(is.finite(fit$theta_array)))
  expect_true(all(fit$theta_array >= 0 & fit$theta_array <= 1))
  expect_true(all(fit$theta_array[,,1] >= fit$theta_array[,,2]))
  expect_true(length(fit$rsf_models) >= 2L)
  first <- fit$models[[1L]]
  H <- as.numeric(d$time > .75) / fit$censoring_G_train[,1]
  expected <- dr_pseudo_outcome(H, d$a, ps$treatment_levels[1L],
                                ps$ehat[, ps$treatment_levels[1L]],
                                first$first_stage_prediction)
  expect_equal(first$pseudo_outcome, expected, tolerance = 1e-10)
  expect_equal(first$second_seed, 901L + 2000L + 1L)
})
