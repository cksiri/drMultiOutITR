# Synthetic two-treatment example with binary and continuous outcomes.
set.seed(2026)
n <- 160L
x <- rnorm(n)
z <- runif(n, -1, 1)
a <- ifelse(runif(n) < plogis(.2 + .4 * x), 'A', 'B')
ybin <- rbinom(n, 1, plogis(.2 + .3 * x + .5 * (a == 'B')))
ycont <- .4 * x - .3 * z + .6 * (a == 'A') + rnorm(n)
d <- data.frame(a, x, z, ybin, ycont)

spec <- list(
  nonuse = list(type = 'binary', column = 'ybin', favorable = 1),
  benefit = list(type = 'continuous', column = 'ycont',
                 direction = 'higher', probs = c(.25, .50, .75),
                 grid_length = 7L)
)

fit <- fit_dr_itr(
  data = d, treatment = 'a', covariates = c('x', 'z'),
  outcomes = spec, weights = c(nonuse = .4, benefit = .6),
  nuisance_folds = 2L, num.trees = 60L, propensity_trees = 60L,
  min.node.size = 5L, seed = 2026L
)
print(fit)

new_patients <- data.frame(x = c(-.5, 0, .5), z = c(.2, -.1, 0))
result <- predict(fit, new_patients, ids = c('P1', 'P2', 'P3'),
                  return_summaries = TRUE)
print(result$selection)
print(head(result$summaries, 8L))

# Reprioritize outcomes without refitting the conditional models or scaling.
prioritize_nonuse <- predict(fit, new_patients,
                             weights = c(nonuse = .7, benefit = .3))
print(prioritize_nonuse)
