test_that('candidate magnitude is indexed by rank position, not treatment label', {
  tab <- data.frame(new_id = 'p', outcome_id = 'benefit', summary_id = 's',
                    treatment = c('A','B','C'), theta_hat = c(3,2,1),
                    larger_is_better = TRUE, outcome_weight = 1,
                    rho = 0.5, summary_weight = 1, summary_scale = 1)
  best <- drMultiOutITR:::compute_pi_one_patient(tab, c(1,2,3),
                                                  c('A','B','C'), 'benefit')
  worst <- drMultiOutITR:::compute_pi_one_patient(tab, c(3,2,1),
                                                   c('A','B','C'), 'benefit')
  expect_equal(as.numeric(best), 0)
  expect_equal(as.numeric(worst), 4 * sqrt(2), tolerance = 1e-12)
})

test_that('frozen fit/predict exactly match original file-19 engine', {
  dat <- make_synthetic_summaries()
  tr <- dat[dat$new_id %in% unique(dat$new_id)[1:18], , drop = FALSE]
  te <- dat[dat$new_id %in% unique(dat$new_id)[19:24], , drop = FALSE]
  fit <- fit_rank_rule(tr)
  pred <- predict(fit, te, details = TRUE)
  ref <- rank_treatments(te, tr)
  expect_equal(nrow(fit$scaling_fit$candidates), 6L)
  expect_equal(pred$selection, ref$selection, tolerance = 1e-10)
  expect_equal(pred$details, ref$details, tolerance = 1e-10)
  expect_true(all(is.finite(pred$selection$min_W)))
  expect_true(all(pred$selection$best_treatment %in% c('LTR','STI','ODF')))
})

test_that('reweighting retains frozen covariance matrices', {
  dat <- make_synthetic_summaries()
  fit <- fit_rank_rule(dat)
  before <- fit$scaling_fit$scaling[[1]]$Sigma
  w <- c(earnings = 0.5, health = 0.25, nonuse = 0.25)
  pred <- predict(fit, dat, weights = w)
  expect_equal(nrow(pred), 24L)
  expect_identical(fit$scaling_fit$scaling[[1]]$Sigma, before)
})

test_that('malformed summary tables are rejected', {
  dat <- make_synthetic_summaries()
  bad <- dat[-1, , drop = FALSE]
  expect_error(fit_rank_rule(bad), 'incomplete')
  fit <- fit_rank_rule(dat)
  unknown <- dat
  unknown$summary_id[1] <- 'different'
  expect_error(predict(fit, unknown))
  expect_error(predict(fit, dat, weights = c(earnings=1, nonuse=1)))
})
