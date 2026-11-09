test_that('conditional AIPW pseudo-outcome implements equation correctly', {
  p <- dr_pseudo_outcome(transformed = c(1, 0, 1),
                         received = c('A', 'B', 'A'), target = 'A',
                         propensity = c(.5, .4, .25),
                         outcome_prediction = c(.2, .4, .6))
  expect_equal(p, c(1.8, .4, 2.2))
  expect_error(dr_pseudo_outcome(1, 'A', 'A', 0, .2), 'propensity')
})

test_that('propensity fits and out-of-fold probabilities retain named arms', {
  d <- make_dr_test_data()
  p <- fit_dr_propensity(d, 'a', c('x','z'), method = 'multinomial',
                         folds = 3, eps = .05)
  expect_s3_class(p, 'dr_propensity')
  expect_identical(colnames(p$ehat), c('A','B'))
  expect_equal(rowSums(p$ehat), rep(1, nrow(d)), tolerance = 1e-12)
  expect_equal(sort(unique(p$fold_id)), 1:3)
  expect_true(all(p$ehat >= .05 - 1e-12))
})

test_that('binary, ordinal, continuous estimators return coherent summaries', {
  # More per-arm/per-fold observations reduce accidental logistic separation
  # at the tails of the continuous-outcome CDF grid.
  d <- make_dr_test_data(n = 360L)
  ps <- fit_dr_propensity(d, 'a', c('x','z'), method = 'multinomial', folds = 3)
  args <- list(data = d, treatment = 'a', covariates = c('x','z'),
               propensity = ps, newdata = d[1:6, , drop = FALSE],
               outcome_learner = 'glm', second_stage = 'lm')
  b <- do.call(estimate_dr_binary, c(args, list(outcome = 'ybin')))
  o <- do.call(estimate_dr_ordinal, c(args, list(outcome = 'yord')))
  q <- do.call(estimate_dr_continuous, c(args, list(outcome = 'ycont',
                       probs = c(.25,.5,.75),
                       grid = as.numeric(stats::quantile(d$ycont,
                                     probs = seq(.15,.85,length.out=11))))))
  expect_equal(dim(b$theta_array), c(6L, 2L, 1L))
  expect_equal(dim(o$theta_array), c(6L, 2L, 3L))
  expect_equal(dim(q$theta_array), c(6L, 2L, 3L))
  expect_equal(dim(q$cdf_array), c(6L, 2L, 11L))
  expect_true(all(b$theta_array >= 0 & b$theta_array <= 1))
  expect_true(all(o$theta_array >= 0 & o$theta_array <= 1))
  for (i in 1:6) for (k in 1:2) {
    expect_true(all(diff(o$theta_array[i,k,]) <= 1e-12))
    expect_true(all(diff(q$cdf_array[i,k,]) >= -1e-12))
    expect_true(all(diff(q$theta_array[i,k,]) >= -1e-12))
  }
  expect_equal(b$nuisance_folds, 3L)
  expect_equal(nrow(as_dr_summary_table(b, 'binary')), 12L)
})

test_that('lower-is-better ordinal thresholds are nondecreasing', {
  d <- make_dr_test_data()
  p <- fit_dr_propensity(d, 'a', c('x','z'), method = 'multinomial')
  o <- estimate_dr_ordinal(d, 'a', 'yord', c('x','z'), p,
                           newdata = d[1:3,], direction = 'lower',
                           outcome_learner = 'glm', second_stage = 'lm')
  expect_true(all(diff(o$theta_array[1,1,]) >= -1e-12))
})

test_that('survival IPCW and conditional summaries obey probability constraints', {
  skip_if_not_installed('survival')
  d <- make_dr_test_data()
  p <- fit_dr_propensity(d, 'a', c('x','z'), method = 'multinomial', folds = 3)
  s <- estimate_dr_survival(d, 'a', 'time', 'status', c('x','z'), p,
                            newdata = d[1:5,], times = c(1,2),
                            censoring = 'pooled_km',
                            outcome_learner = 'lm', second_stage = 'lm')
  expect_equal(dim(s$theta_array), c(5L, 2L, 2L))
  expect_true(all(s$theta_array >= 0 & s$theta_array <= 1))
  expect_true(all(s$theta_array[1,1,1] >= s$theta_array[1,1,2]))
  expect_equal(nrow(as_dr_summary_table(s, 'survival')), 20L)
  sc <- estimate_dr_survival(d, 'a', 'time', 'status', c('x','z'), p,
                             newdata = d[1:3,], times = c(1,2),
                             censoring = 'cox', outcome_learner = 'lm',
                             second_stage = 'lm')
  expect_true(all(is.finite(sc$theta_array)))
})

test_that('DR summaries feed directly into the existing rank engine', {
  # Keep GLM indicator fits well populated for each threshold and fold.
  d <- make_dr_test_data(n = 360L)
  ps <- fit_dr_propensity(d, 'a', c('x','z'), method = 'multinomial', folds = 3)
  b <- estimate_dr_binary(d, 'a', 'ybin', c('x','z'), ps,
                          outcome_learner = 'glm', second_stage = 'lm')
  q <- estimate_dr_continuous(d, 'a', 'ycont', c('x','z'), ps,
                              grid = as.numeric(stats::quantile(d$ycont,
                                  probs = seq(.15,.85,length.out=9))),
                              outcome_learner = 'glm',
                              second_stage = 'lm')
  tab <- rbind(as_dr_summary_table(b, 'nonuse', outcome_weight = .5),
               as_dr_summary_table(q, 'earnings', outcome_weight = .5))
  training <- tab[tab$new_id %in% as.character(1:300), , drop=FALSE]
  new <- tab[tab$new_id %in% as.character(301:360), , drop=FALSE]
  rank <- fit_rank_rule(training)
  result <- predict(rank, new)
  expect_equal(nrow(result), 60L)
  expect_true(all(is.finite(result$min_W)))
  expect_true(all(result$best_treatment %in% c('A','B')))
})
