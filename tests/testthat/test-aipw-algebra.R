test_that("three-arm AIPW calculation is independently correct", {
  trt <- factor(c("LTR", "STI", "ODF", "LTR", "STI", "ODF"),
                levels = c("LTR", "STI", "ODF"))
  h <- c(1, 0, 1, 0, 1, 0)
  m <- c(.3, .4, .5, .2, .8, .6)
  e <- c(.4, .3, .3, .5, .2, .25)
  for (k in levels(trt)) {
    observed <- dr_pseudo_outcome(h, trt, k, e, m)
    independent <- numeric(length(h))
    for (i in seq_along(h)) {
      independent[i] <- m[i]
      if (as.character(trt[i]) == k) {
        independent[i] <- m[i] + (h[i] - m[i]) / e[i]
      }
    }
    expect_equal(observed, independent, tolerance = 1e-14)
    expect_equal(observed[as.character(trt) != k],
                 m[as.character(trt) != k])
  }
})

test_that("three-arm multinomial nuisance folds remain valid", {
  skip_if_not_installed("nnet")
  set.seed(4532)
  n <- 270L
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$trt <- factor(sample(c("LTR", "STI", "ODF"), n, replace = TRUE),
                  levels = c("LTR", "STI", "ODF"))
  p <- fit_dr_propensity(d, "trt", c("x", "z"),
                         method = "multinomial", folds = 3L,
                         eps = 0.05, seed = 2026)
  expect_identical(colnames(p$ehat), levels(d$trt))
  expect_equal(rowSums(p$ehat), rep(1, n), tolerance = 1e-12)
  expect_true(all(p$ehat >= .05 - 1e-12))
  expect_setequal(unique(p$fold_id), seq_len(3L))
  for (f in unique(p$fold_id)) {
    expect_setequal(unique(as.character(d$trt[p$fold_id != f])),
                    levels(d$trt))
  }
})

test_that("conditional CDF inversion agrees with hand calculations", {
  cdf <- c(.1, .2, .75, .75, 1)
  y_grid <- c(-2, -1, 0, 1, 2)
  probs <- c(.25, .50, .90)
  expect_equal(drMultiOutITR:::.dr_quantile_inverse(cdf, y_grid, probs, "grid"), c(0, 0, 2))
  interp <- drMultiOutITR:::.dr_quantile_inverse(cdf, y_grid, probs, "interpolate")
  expect_equal(length(interp), 3L)
  expect_true(all(is.finite(interp)))
  expect_true(all(diff(interp) >= 0))
})

test_that("pooled censoring G matches independent KM computation", {
  skip_if_not_installed("survival")
  d <- data.frame(trt = factor(rep(c("LTR", "STI", "ODF"), each = 4L)),
                  x = seq(-1, 1, length.out = 12L),
                  time = c(.5,.7,1.2,1.8,.6,1.1,1.5,2.1,.8,1.3,1.7,2.5),
                  status = c(1,0,1,0,1,1,0,1,0,1,1,0))
  times <- c(.9, 1.4)
  fit <- drMultiOutITR:::.dr_censor_fit(d, "trt", "x", d$time, d$status, "pooled_km")
  got <- drMultiOutITR:::.dr_censor_predict(fit, d, "trt", "x", times, min_G = .05)
  independent <- as.numeric(summary(
    survival::survfit(survival::Surv(d$time, 1 - d$status) ~ 1),
    times = times, extend = TRUE)$surv)
  independent <- pmax(independent, .05)
  expect_equal(got, matrix(rep(independent, each = nrow(d)),
                           nrow = nrow(d)), tolerance = 1e-12)
  ipcw <- as.numeric(d$time > times[1L]) / got[, 1L]
  manual <- ifelse(d$time > times[1L], 1 / independent[1L], 0)
  expect_equal(ipcw, manual, tolerance = 1e-12)
})


test_that("returned survival IDs locate all stored pseudo-outcomes", {
  skip_if_not_installed("survival")
  d <- make_dr_test_data(n = 180L)
  p <- fit_dr_propensity(d, "a", c("x", "z"),
                         method = "multinomial", folds = 1L)
  times <- c(0.75, 1.5, 2.5)
  fit <- estimate_dr_survival(
    d, "a", "time", "status", c("x", "z"), p,
    newdata = d[1:4, , drop = FALSE], times = times,
    censoring = "pooled_km", outcome_learner = "lm",
    second_stage = "lm", return_models = TRUE)
  expect_length(fit$summary_ids, length(times))
  expect_identical(fit$summary_ids, dimnames(fit$theta_array)[[3L]])
  expected_names <- unlist(lapply(fit$summary_ids, function(label) {
    paste0(label, "__", fit$treatment_levels)
  }), use.names = FALSE)
  expect_true(all(expected_names %in% names(fit$models)))
  expect_true(all(vapply(fit$models[expected_names], function(x) {
    is.numeric(x$pseudo_outcome) && length(x$pseudo_outcome) == nrow(d)
  }, logical(1L))))
})
