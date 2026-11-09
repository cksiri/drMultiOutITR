# Synthetic example: binary and continuous DR summaries followed by ranking.
# This illustrative fit is not a held-out policy-value assessment.
library(drMultiOutITR)
set.seed(91)
n <- 160L
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
A <- ifelse(runif(n) < plogis(-.1 + .4 * X$x1), 'A', 'B')
dat <- within(X, {
  arm <- A
  binary <- rbinom(n, 1, plogis(.3 + .4 * x1 + .3 * (A == 'B')))
  benefit <- x1 + .5 * x2 + .7 * (A == 'B') + rnorm(n)
})
training <- dat[1:130, , drop = FALSE]
all_profiles <- dat[, c('x1','x2'), drop = FALSE]
# Three stratified nuisance folds: propensity and binary/CDF first stage
# are predicted outside each participant's nuisance fold.
ps <- fit_dr_propensity(training, 'arm', c('x1','x2'),
                        method = 'multinomial', folds = 3)
b <- estimate_dr_binary(training, 'arm', 'binary', c('x1','x2'), ps,
                        newdata = all_profiles,
                        outcome_learner = 'glm', second_stage = 'lm')
q <- estimate_dr_continuous(training, 'arm', 'benefit', c('x1','x2'), ps,
                            newdata = all_profiles,
                            grid = as.numeric(stats::quantile(training$benefit,
                                       probs = seq(.15,.85,length.out=15))),
                            outcome_learner = 'glm', second_stage = 'lm')
tab <- rbind(as_dr_summary_table(b, 'binary', outcome_weight = .5),
             as_dr_summary_table(q, 'benefit', outcome_weight = .5))
train_ids <- as.character(seq_len(130))
test_ids <- as.character(131:160)
rank_fit <- fit_rank_rule(tab[tab$new_id %in% train_ids, , drop = FALSE])
print(head(predict(rank_fit, tab[tab$new_id %in% test_ids, , drop = FALSE])))
# This simple example uses in-sample covariance estimates. Research analyses
# may require a separate held-out assessment.
