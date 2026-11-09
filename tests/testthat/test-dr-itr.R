test_that('continuous grid arguments distinguish a length from thresholds', {
  d <- make_dr_test_data(n = 60L)
  count_spec <- drMultiOutITR:::.dr_itr_specs(
    list(lab = list(type = 'continuous', column = 'ycont', grid_length = 7L)),
    d
  )
  expect_null(count_spec$lab[["grid", exact = TRUE]])
  expect_identical(count_spec$lab$grid_length, 7L)

  explicit_grid <- c(-2, 0, 2)
  vector_spec <- drMultiOutITR:::.dr_itr_specs(
    list(lab = list(type = 'continuous', column = 'ycont',
                    grid = explicit_grid, grid_length = 7L)), d
  )
  expect_identical(vector_spec$lab[["grid", exact = TRUE]], explicit_grid)
})

test_that('integrated forest fit predicts mixed outcomes with frozen scaling', {
  skip_if_not_installed('ranger')
  d <- make_dr_test_data(n = 180L)
  d$sex <- factor(rep(c('F', 'M'), length.out = nrow(d)))
  specs <- list(
    response = list(type = 'binary', column = 'ybin', favorable = 1),
    health = list(type = 'ordinal', column = 'yord', direction = 'higher'),
    lab = list(type = 'continuous', column = 'ycont',
               direction = 'lower', probs = c(.25, .5, .75), grid_length = 7L)
  )
  fit <- fit_dr_itr(d, treatment = 'a', covariates = c('x', 'z', 'sex'),
                    outcomes = specs, weights = c(.2, .3, .5),
                    propensity_method = 'ranger', nuisance_folds = 2L,
                    num.trees = 50L, propensity_trees = 55L,
                    min.node.size = 6L, seed = 91L)
  expect_s3_class(fit, 'dr_itr')
  expect_s3_class(fit$rank_rule, 'dr_rank_rule')
  expect_identical(fit$covariance_source, 'oob')
  expect_equal(fit$n_train, 180L)
  expect_equal(fit$weights, c(response = .2, health = .3, lab = .5))

  new <- d[1:4, c('x', 'z', 'sex'), drop = FALSE]
  pred <- predict(fit, new, ids = paste0('person-', 1:4),
                  details = TRUE, return_summaries = TRUE)
  expect_equal(nrow(pred$selection), 4L)
  manual_selection <- stats::predict(fit$rank_rule, pred$summaries)
  expect_equal(manual_selection, pred$selection)
  expect_true(all(pred$selection$best_treatment %in% c('A', 'B')))
  expect_true(all(is.finite(pred$selection$min_W)))
  expect_true(all(c('person-1', 'person-4') %in% pred$selection$new_id))
  expect_equal(nrow(pred$summaries), 4L * 2L * (1L + 3L + 3L))
  expect_true(all(pred$summaries$theta_hat[ pred$summaries$outcome_id == 'response'] >= 0))
  expect_true(all(pred$summaries$theta_hat[ pred$summaries$outcome_id == 'response'] <= 1))
  expect_false(unique(pred$summaries$larger_is_better[
    pred$summaries$outcome_id == 'lab']))
  expect_identical(pred$selection$best_treatment,
                   predict(fit, new, ids = paste0('person-', 1:4))$best_treatment)
  expect_equal(pred$selection$min_W,
               predict(fit, new, weights = c(lab = .5, health = .3,
                                              response = .2))$min_W,
               tolerance = 1e-8)
  expect_error(predict(fit, transform(new, sex = 'unknown')), 'Unknown level')
  expect_error(predict(fit, new, ids = rep('same', 4)), 'ids')

  # Saved models must preserve predictions and alternative-weight rankings.
  tmp_model <- tempfile(pattern = 'dr_itr_private_', fileext = '.rds')
  saveRDS(fit, tmp_model)
  restored <- readRDS(tmp_model)
  unlink(tmp_model)
  restored_prediction <- predict(restored, new, ids = paste0('person-', 1:4),
                                 details = TRUE, return_summaries = TRUE)
  expect_identical(restored_prediction$selection, pred$selection)
  expect_identical(restored_prediction$details, pred$details)
  expect_identical(restored_prediction$summaries, pred$summaries)
  original_weighted <- predict(fit, new, weights = c(response=.8, health=.1, lab=.1))
  restored_weighted <- predict(restored, new,
                              weights = c(response=.8, health=.1, lab=.1))
  expect_identical(restored_weighted, original_weighted)
})

test_that('integrated linear second stage requires explicit in-sample scaling', {
  d <- make_dr_test_data(n = 360L)
  specs <- list(nonuse = list(type = 'binary', column = 'ybin'),
                health = list(type = 'ordinal', column = 'yord'))
  expect_error(
    fit_dr_itr(d, 'a', c('x', 'z'), specs, second_stage = 'lm',
               propensity_method = 'multinomial', nuisance_folds = 2L,
               outcome_learner = 'glm'),
    'requires second_stage')
  fit <- fit_dr_itr(d, 'a', c('x', 'z'), specs,
                    weights = c(nonuse = .6, health = .4),
                    second_stage = 'lm', covariance_source = 'in_sample',
                    propensity_method = 'multinomial', nuisance_folds = 2L,
                    outcome_learner = 'glm')
  res <- predict(fit, d[1:3, c('x','z')], return_summaries = TRUE)
  expect_equal(nrow(res$selection), 3L)
  stored <- drMultiOutITR:::.dr_itr_model_predictions(
    fit$estimates$nonuse, d, c('x', 'z'), seq_len(nrow(d)))
  expect_equal(stored$theta_array, fit$estimates$nonuse$theta_array, tolerance = 1e-10)
  expect_equal(nrow(res$summaries), 3L * 2L * (1L + 3L))
  expect_error(predict(fit, d[1:2, c('x','z')], weights = c(wrong=.2, health=.8)),
               'Named weights')
})

test_that('survival summaries use stored second stages at prediction time', {
  skip_if_not_installed('survival')
  skip_if_not_installed('ranger')
  d <- make_dr_test_data(n = 170L)
  specs <- list(survival = list(type = 'survival', time = 'time',
                                status = 'status', times = c(.5, 1.5),
                                censoring = 'pooled_km'),
                nonuse = list(type = 'binary', column = 'ybin'))
  fit <- fit_dr_itr(d, 'a', c('x','z'), specs,
                    propensity_method = 'multinomial', nuisance_folds = 2L,
                    outcome_learner = 'ranger', second_stage = 'lm',
                    covariance_source = 'in_sample',
                    num.trees = 40L, min.node.size = 5L)
  res <- predict(fit, d[1:3, c('x','z')], return_summaries = TRUE)
  expect_equal(nrow(res$summaries), 3L * 2L * (2L + 1L))
  s <- res$summaries[res$summaries$outcome_id == 'survival', ]
  expect_true(all(s$theta_hat >= 0 & s$theta_hat <= 1))
  for (id in unique(s$new_id)) for (a in unique(s$treatment)) {
    v <- s$theta_hat[s$new_id == id & s$treatment == a]
    expect_true(all(diff(v) <= 1e-12))
  }
})

test_that('input validation does not silently accept invalid clinical choices', {
  d <- make_dr_test_data(n = 60L)
  expect_error(fit_dr_itr(d, 'a', c('x','z'), list()), 'outcomes')
  expect_error(fit_dr_itr(d, 'a', c('x','z'), list(a = list(type='other', column='ybin'))),
               'type')
  expect_error(fit_dr_itr(d, 'a', c('x','z'), list(a = list(type='binary', column='missing'))),
               'existing column')
  expect_error(drMultiOutITR:::.dr_itr_weights(c(-1,2), c('a','b')), 'nonnegative')
  expect_equal(drMultiOutITR:::.dr_itr_weights(NULL, c('a','b')), c(a=.5,b=.5))
  valid_spec <- list(nonuse = list(type = 'binary', column = 'ybin'))
  expect_error(fit_dr_itr(d, 'a', c('x', 'not_a_covariate'), valid_spec),
               'Missing columns')
  expect_error(fit_dr_itr(d, 'not_a_treatment', c('x', 'z'), valid_spec),
               'Missing columns')
  missing_baseline <- d
  missing_baseline$x[1L] <- NA_real_
  expect_error(fit_dr_itr(missing_baseline, 'a', c('x', 'z'), valid_spec),
               'Missing values')
})


test_that('three-treatment fit enumerates all candidate permutations', {
  skip_if_not_installed('ranger')
  d <- make_dr_test_data(n = 240L)
  set.seed(86)
  d$a <- sample(c('A', 'B', 'C'), nrow(d), replace = TRUE)
  specs <- list(nonuse = list(type = 'binary', column = 'ybin'),
                health = list(type = 'ordinal', column = 'yord'))
  fit <- fit_dr_itr(d, 'a', c('x','z'), specs,
                    propensity_method = 'ranger', nuisance_folds = 2L,
                    num.trees = 45L, propensity_trees = 50L,
                    min.node.size = 6L)
  expect_equal(length(fit$treatments), 3L)
  expect_equal(nrow(fit$rank_rule$scaling_fit$candidates), 6L)
  p <- predict(fit, d[1:4, c('x','z')], details = TRUE)
  expect_equal(nrow(p$selection), 4L)
  expect_equal(nrow(p$details), 4L * 6L)
  expect_true(all(p$selection$best_treatment %in% c('A','B','C')))
})
