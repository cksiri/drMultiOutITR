test_that("three-arm candidate discordance follows the manuscript formula", {
  # Independent hand calculation for two one-component outcomes.
  d <- data.frame(
    new_id = rep("reference", 6L),
    outcome_id = rep(c("benefit", "harm"), each = 3L),
    summary_id = rep("component", 6L),
    treatment = rep(c("A", "B", "C"), 2L),
    theta_hat = c(4, 9, 1, 5, 2, 9),
    larger_is_better = rep(c(TRUE, FALSE), each = 3L),
    outcome_weight = rep(0.5, 6L),
    rho = rep(0.5, 6L),
    summary_weight = rep(1, 6L),
    summary_scale = rep(1, 6L)
  )
  # The observed order is B, A, C for BOTH outcomes:
  # rank vector (2,1,3), so that candidate has Pi=0 exactly.
  get_pi <- drMultiOutITR:::compute_pi_one_patient
  agreement <- get_pi(d, c(2L, 1L, 3L), c("A", "B", "C"),
                      c("benefit", "harm"))
  expect_equal(unname(agreement), c(0, 0), tolerance = 1e-12)

  # Candidate (1,2,3) swaps A and B:
  # higher-is-better: 2*sqrt(|4-9|)=2*sqrt(5)
  # lower-is-better:  2*sqrt(|5-2|)=2*sqrt(3)
  disagreement <- get_pi(d, c(1L, 2L, 3L), c("A", "B", "C"),
                         c("benefit", "harm"))
  expect_equal(unname(disagreement), c(2*sqrt(5), 2*sqrt(3)),
               tolerance = 1e-12)
})

test_that("frozen training covariance does not depend on heldout profiles", {
  make_table <- function(ids, shift = 0) {
    out <- vector("list", length(ids) * 3L * 2L)
    pos <- 0L
    for (i in seq_along(ids)) for (k in seq_along(c("A", "B", "C"))) {
      for (j in seq_along(c("benefit", "harm"))) {
        pos <- pos + 1L
        out[[pos]] <- data.frame(
          new_id = ids[i], outcome_id = c("benefit", "harm")[j],
          summary_id = "only", treatment = c("A", "B", "C")[k],
          theta_hat = if (j == 1L) sin(i / 3 + k) + 0.25*k + shift
                      else cos(i / 4 - k) + 0.3*k - shift,
          larger_is_better = (j == 1L), outcome_weight = .5,
          rho = .5, summary_weight = 1, summary_scale = 1
        )
      }
    }
    do.call(rbind, out)
  }
  tr <- make_table(paste0("tr", seq_len(24)))
  te1 <- make_table(paste0("te", seq_len(5)))
  te2 <- make_table(paste0("te", seq_len(5)), shift = 100)
  rule <- fit_rank_rule(tr)
  before <- lapply(rule$scaling_fit$scaling, function(x) x$Sigma)
  invisible(predict(rule, te1))
  invisible(predict(rule, te2))
  after <- lapply(rule$scaling_fit$scaling, function(x) x$Sigma)
  expect_equal(after, before, tolerance = 0)
  expect_equal(nrow(rule$scaling_fit$candidates), 6L)
  expect_equal(length(after), 6L)
})
