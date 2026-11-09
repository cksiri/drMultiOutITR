make_dr_test_data <- function(n = 90L) {
  set.seed(126)
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  a <- ifelse(stats::runif(n) < stats::plogis(.3 + .2 * x), 'A', 'B')
  p <- stats::plogis(-.2 + .4 * x + .3 * (a == 'B'))
  ybin <- stats::rbinom(n, 1L, p)
  latent <- .2 * x + .1 * (a == 'B') + stats::rnorm(n, sd = 1.3)
  yord <- ordered(cut(latent, breaks = c(-Inf, -1, 0, 1, Inf),
                      labels = c('poor','fair','good','great')),
                  levels = c('poor','fair','good','great'))
  ycont <- .5 + x + .2 * z + (a == 'B') * .7 + stats::rnorm(n)
  event <- stats::rexp(n, rate = .2 * exp(-.2 * x))
  censor <- stats::rexp(n, rate = .15)
  data.frame(a = a, x = x, z = z, ybin = ybin, yord = yord,
             ycont = ycont, time = pmin(event, censor),
             status = as.integer(event <= censor))
}
