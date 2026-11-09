test_that('constant CDF inversion follows ORIGINAL first-grid convention', {
  inv <- drMultiOutITR:::.dr_quantile_inverse
  grid <- seq(-3, 3, length.out = 41L)
  probs <- c(0.25, 0.5, 0.75)
  for (level in c(0, .5, 1)) {
    flat <- rep(level, length(grid))
    # Original 03_continuous_quantile_dr.R convention:
    # y_grid[which.min(abs(cdf - p_vec[1]))] on an all-equal CDF.
    original_value <- grid[which.min(abs(flat - probs[[1L]]))]
    expect_identical(inv(flat, grid, probs, 'interpolate'),
                     rep(original_value, length(probs)))
  }
})

test_that('nonconstant CDF interpolation and grid inversion remain unchanged', {
  inv <- drMultiOutITR:::.dr_quantile_inverse
  grid <- seq(-2, 2, length.out = 5L)
  probs <- c(.25, .50, .75)
  cdf <- c(0, .1, .4, .4, 1)
  keep <- !duplicated(cdf)
  original_interp <- as.numeric(stats::approx(
    x = cdf[keep], y = grid[keep], xout = probs,
    rule = 2, ties = 'ordered')$y)
  expect_equal(inv(cdf, grid, probs, 'interpolate'), original_interp,
               tolerance = 0)
  expect_identical(inv(cdf, grid, probs, 'grid'),
                   vapply(probs, function(q) {
                     loc <- which(cdf >= q)
                     if (length(loc)) grid[loc[[1L]]] else grid[length(grid)]
                   }, numeric(1)))
})
