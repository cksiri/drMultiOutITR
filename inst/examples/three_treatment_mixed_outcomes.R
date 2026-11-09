# Three-treatment example with continuous, ordinal, and survival outcomes.
# Simulated data and small forests are used for illustration.

library(drMultiOutITR)
set.seed(2026)
n <- 300L

# 1. Simulate observational data.
# Treatment assignment is related to baseline characteristics (not randomized).
age <- round(runif(n, 40, 80), 1)
risk <- rnorm(n)
group <- factor(sample(c("Group_1", "Group_2"), n, replace = TRUE))
age10 <- (age - 60) / 10

scores <- cbind(
  Care_A = 0,
  Care_B = 0.25 + 0.40 * risk - 0.15 * age10,
  Care_C = -0.15 - 0.30 * risk + 0.20 * age10
)
probabilities <- exp(scores - apply(scores, 1L, max))
probabilities <- probabilities / rowSums(probabilities)
treatment <- vapply(seq_len(n), function(i) {
  sample(colnames(probabilities), 1L, prob = probabilities[i, ])
}, character(1L))
treatment <- factor(treatment, levels = colnames(probabilities))

# Continuous outcome: fewer days are preferable.
los_days <- pmax(1, round(8 + 0.8 * age10 + 1.3 * risk +
  0.6 * (group == "Group_2") - 0.9 * (treatment == "Care_B") +
  0.5 * (treatment == "Care_C") + rnorm(n, sd = 2.8), 1))

# Ordinal outcome: higher categories are preferable.
latent <- -0.35 * risk - 0.10 * age10 + 0.35 * (treatment == "Care_C") +
  rnorm(n)
wellbeing <- ordered(cut(latent, breaks = c(-Inf, -0.5, 0.5, Inf),
  labels = c("Low", "Medium", "High")), levels = c("Low", "Medium", "High"))

# Right-censored survival: status = 1 means an observed event.
event_time <- rexp(n, rate = 0.014 * exp(0.35 * risk +
  0.15 * age10 - 0.18 * (treatment == "Care_B")))
censor_time <- rexp(n, rate = 0.005)
followup_days <- pmin(event_time, censor_time, 180)
status <- as.integer(event_time <= censor_time & event_time <= 180)

study <- data.frame(treatment, age, risk, group, los_days,
  wellbeing, followup_days, status)
print(table(study$treatment))

# 2. Specify outcomes and favorable directions.
responses <- list(
  stay = list(type = "continuous", column = "los_days",
              direction = "lower", probs = c(0.25, 0.50, 0.75),
              grid_length = 7L),
  wellbeing = list(type = "ordinal", column = "wellbeing",
                   direction = "higher"),
  survival = list(type = "survival", time = "followup_days",
                  status = "status", times = c(30, 90, 150),
                  censoring = "pooled_km", min_G = 0.05)
)

# 3. Fit the treatment-ranking model.
fit <- fit_dr_itr(
  data = study,
  treatment = "treatment",
  covariates = c("age", "risk", "group"),
  outcomes = responses,
  weights = c(stay = 1/3, wellbeing = 1/3, survival = 1/3),
  nuisance_folds = 2L,
  propensity_method = "ranger",
  outcome_learner = "ranger",
  num.trees = 60L,
  propensity_trees = 60L,
  min.node.size = 8L,
  seed = 2026L
)
print(fit)

# 4. Rank treatments for new baseline profiles.
new_profiles <- data.frame(
  age = c(45, 55, 65, 75),
  risk = c(-0.8, 0.2, 1.0, -0.2),
  group = factor(c("Group_1", "Group_2", "Group_1", "Group_2"),
                 levels = levels(study$group))
)
profile_ids <- paste0("Person_", seq_len(nrow(new_profiles)))
result <- predict(fit, new_profiles, ids = profile_ids,
                  return_summaries = TRUE)
print(result$selection[, c("new_id", "best_treatment", "optimal_ranking")])
print(head(result$summaries[, c("new_id", "outcome_id", "summary_id",
                               "treatment", "theta_hat")], 12L))
# optimal_ranking gives the rank of each treatment in fit$treatments order.
# For example, 2-3-1 means Care_A second, Care_B third, Care_C first.
print(fit$treatments)
readable_ranking <- vapply(
  strsplit(result$selection$optimal_ranking, "-", fixed = TRUE),
  function(x) paste(fit$treatments[order(as.integer(x))], collapse = " > "),
  character(1L)
)
print(data.frame(new_id = result$selection$new_id,
                 best_treatment = result$selection$best_treatment,
                 ordered_treatments = unname(readable_ranking)))

# 5. Change outcome weights without refitting the forests.
prioritize_survival <- predict(
  fit, new_profiles, ids = profile_ids,
  weights = c(stay = 0.20, wellbeing = 0.20, survival = 0.60)
)
print(prioritize_survival[, c("new_id", "best_treatment")])

# 6. Save and reload the fitted model.
model_file <- tempfile(fileext = ".rds")
saveRDS(fit, model_file)
restored <- readRDS(model_file)
restored_result <- predict(restored, new_profiles, ids = profile_ids,
                           return_summaries = TRUE)
stopifnot(identical(result$selection, restored_result$selection))
stopifnot(identical(result$summaries, restored_result$summaries))
unlink(model_file)
message("Three-treatment mixed-outcome example completed successfully.")
