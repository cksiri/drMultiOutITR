# Fully synthetic examples; no RHC or DATOS participant records.
make_synthetic_summaries <- function(n = 24L) {
  ids <- sprintf('synthetic_%03d', seq_len(n))
  x <- seq(-1.3, 1.3, length.out = n)
  trt <- c('LTR', 'STI', 'ODF')
  rows <- vector('list', n * 3L * 3L)
  ii <- 0L
  for (j in seq_len(n)) for (k in seq_along(trt)) {
    earnings <- c(3500, 5000, 4000)[k] +
      c(2900, -2400, 1650)[k] * x[j] + 250 * sin((k + 1) * x[j])
    eta <- c(-0.3, 0.3, 0)[k] + c(0.7, -0.6, 0.15)[k] * x[j]
    pp <- c(earnings, plogis(eta - 0.8),
            plogis(c(0.3, -0.4, 0.15)[k] +
                   c(-0.5, 0.45, 0.9)[k] * x[j]))
    for (h in seq_along(pp)) {
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        new_id = ids[j], outcome_id = c('earnings','health','nonuse')[h],
        summary_id = c('median','excellent','no_use')[h],
        treatment = trt[k], theta_hat = pp[h],
        larger_is_better = TRUE, outcome_weight = 1 / 3, rho = 0.5,
        summary_weight = 1, summary_scale = 1,
        stringsAsFactors = FALSE
      )
    }
  }
  tab <- do.call(rbind, rows)
  tab$treatment <- factor(tab$treatment, levels = trt)
  tab
}
