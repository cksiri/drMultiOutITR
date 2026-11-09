# Ranking from synthetic conditional outcome summaries.
# No causal nuisance models are fitted in this example.
library(drMultiOutITR)
ids <- sprintf('synthetic_%02d', 1:18)
trt <- c('A', 'B', 'C')
set.seed(42)
rows <- expand.grid(new_id = ids,
                    outcome_id = c('benefit', 'adverse_event_risk'),
                    treatment = trt, KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
xx <- match(rows$new_id, ids) / length(ids)
kk <- match(rows$treatment, trt)
rows$theta_hat <- ifelse(rows$outcome_id == 'benefit',
                         1 + kk * xx + 0.1 * sin(kk + 5 * xx),
                         0.1 + (4 - kk) * xx / 6 + 0.01 * cos(5 * xx + kk))
rows$summary_id <- 'primary'
rows$larger_is_better <- rows$outcome_id == 'benefit'
rows$outcome_weight <- 0.5
rows$rho <- 0.5
rows$summary_weight <- 1
rows$summary_scale <- 1
train <- rows[rows$new_id %in% ids[1:12], ]
test <- rows[rows$new_id %in% ids[13:18], ]
model <- fit_rank_rule(train)
print(predict(model, test))
print(predict(model, test,
              weights = c(benefit = 0.7, adverse_event_risk = 0.3)))
