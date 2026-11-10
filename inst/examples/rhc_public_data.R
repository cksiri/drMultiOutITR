# Public-data RHC example for drMultiOutITR.
# Source this file, then call run_rhc_public_example().
# No data files, fitted models, or individual predictions are written.

prepare_rhc_public_data <- function(csv_file) {
  if (!is.character(csv_file) || length(csv_file) != 1L ||
      is.na(csv_file) || !file.exists(csv_file)) {
    stop("Provide a path to the downloaded Vanderbilt rhc.csv file.",
         call. = FALSE)
  }

  d <- utils::read.csv(csv_file, stringsAsFactors = FALSE,
                       check.names = FALSE, na.strings = c("", "NA"))

  # The CSV has a leading, unnamed row-number column.
  if (ncol(d) && names(d)[1L] %in% c("", "X", "...1") &&
      is.numeric(d[[1L]]) &&
      identical(as.numeric(d[[1L]]), as.numeric(seq_len(nrow(d))))) {
    d[[1L]] <- NULL
  }

  covariates <- c(
    "cat1", "ca", "cardiohx", "chfhx", "dementhx", "psychhx",
    "chrpulhx", "renalhx", "liverhx", "gibledhx", "malighx", "immunhx",
    "transhx", "amihx", "age", "sex", "edu", "race", "income",
    "ninsclas", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1",
    "temp1", "pafi1", "alb1", "hema1", "bili1", "crea1", "sod1",
    "pot1", "paco21", "ph1", "wtkilo1", "dnr1", "resp", "card", "neuro",
    "gastr", "renal", "meta", "hema", "seps", "trauma", "ortho"
  )
  required <- unique(c("swang1", "sadmdte", "dschdte", "dthdte",
                       "lstctdte", covariates))
  missing <- setdiff(required, names(d))
  if (length(missing)) {
    stop("RHC CSV is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (nrow(d) != 5735L) {
    stop("Expected 5,735 rows in the public RHC CSV. Found ", nrow(d),
         "; verify the input file.", call. = FALSE)
  }

  d$A_label <- factor(d$swang1, levels = c("No RHC", "RHC"))
  if (anyNA(d$A_label)) stop("Unexpected RHC treatment labels.", call. = FALSE)

  d$los_days <- d$dschdte - d$sadmdte
  death_days <- d$dthdte - d$sadmdte
  contact_days <- d$lstctdte - d$sadmdte
  if (any(death_days < 0, na.rm = TRUE) ||
      any(contact_days < 0, na.rm = TRUE) ||
      any(d$los_days < 0, na.rm = TRUE)) {
    stop("Invalid admission, discharge, death, or contact dates.", call. = FALSE)
  }

  d$event180 <- as.integer(!is.na(d$dthdte) &
                             death_days >= 0 & death_days <= 180)
  d$surv_time_180 <- ifelse(
    d$event180 == 1L, death_days,
    ifelse(!is.na(d$dthdte) & death_days > 180, 180,
           pmin(contact_days, 180))
  )

  # One source record lacks a discharge date; LOS cannot be defined.
  d <- d[!is.na(d$los_days), , drop = FALSE]
  rownames(d) <- NULL

  categorical <- c(
    "cat1", "ca", "sex", "race", "income", "ninsclas", "dnr1",
    "resp", "card", "neuro", "gastr", "renal", "meta", "hema",
    "seps", "trauma", "ortho"
  )
  for (v in categorical) d[[v]] <- factor(d[[v]])

  needed <- c("A_label", "los_days", "surv_time_180", "event180", covariates)
  if (nrow(d) != 5734L || anyNA(d[, needed, drop = FALSE]) ||
      any(!is.finite(d$surv_time_180)) ||
      any(d$surv_time_180 <= 0 | d$surv_time_180 > 180)) {
    stop("Unexpected prepared RHC cohort; check the source data.",
         call. = FALSE)
  }
  if (!identical(unname(as.integer(table(d$A_label))), c(3551L, 2183L))) {
    stop("Prepared RHC treatment counts differ from the source analysis.",
         call. = FALSE)
  }

  # Do not retain source identifiers, dates, or other unused patient fields.
  d <- d[, needed, drop = FALSE]
  list(data = d, covariates = covariates)
}

.rhc_example_sample <- function(index, treatment, size, seed) {
  arms <- split(index, as.character(treatment[index]))
  if (length(arms) != 2L || size < 4L || size > length(index)) {
    stop("Invalid stratified-sampling request.", call. = FALSE)
  }
  counts <- vapply(arms, length, integer(1))
  n <- pmin(counts, pmax(2L, floor(size * counts / sum(counts))))
  while (sum(n) < size) {
    j <- which.max(counts - n)
    n[j] <- n[j] + 1L
  }
  while (sum(n) > size) {
    j <- which.max(n - 2L)
    n[j] <- n[j] - 1L
  }
  set.seed(seed)
  sort(unlist(Map(function(ii, nn) sample(ii, nn),
                  arms, as.list(n)), use.names = FALSE))
}

run_rhc_public_example <- function(csv_file = NULL,
                                   n_train = 1600L,
                                   n_holdout = 60L,
                                   num_trees = 60L,
                                   seed = 20261107L) {
  if (!requireNamespace("drMultiOutITR", quietly = TRUE) ||
      !requireNamespace("ranger", quietly = TRUE) ||
      !requireNamespace("randomForestSRC", quietly = TRUE) ||
      !requireNamespace("survival", quietly = TRUE)) {
    stop("Install drMultiOutITR, ranger, randomForestSRC, and survival.",
         call. = FALSE)
  }
  for (x in list(n_train, n_holdout, num_trees)) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != as.integer(x) || x < 4L) {
      stop("n_train, n_holdout, and num_trees must be integers >= 4.",
           call. = FALSE)
    }
  }

  if (is.null(csv_file)) {
    csv_file <- tempfile("rhc_public_", fileext = ".csv")
    on.exit(unlink(csv_file), add = TRUE)
    url <- "https://hbiostat.org/data/repo/rhc.csv"
    cat("Downloading public RHC CSV from Vanderbilt Biostatistics...\n")
    tryCatch(utils::download.file(url, csv_file, mode = "wb", quiet = TRUE),
             error = function(e) stop(
               "Download failed. Save rhc.csv from https://hbiostat.org/data/ ",
               "and pass csv_file='path/to/rhc.csv'.\n", conditionMessage(e),
               call. = FALSE))
  }
  prepared <- prepare_rhc_public_data(csv_file)
  d <- prepared$data
  baseline <- prepared$covariates

  # Form an arm-stratified held-out pool, then small disjoint subsets.
  set.seed(seed)
  by_arm <- split(seq_len(nrow(d)), d$A_label)
  heldout_pool <- sort(unlist(lapply(by_arm, function(ii)
    sample(ii, round(0.30 * length(ii)))), use.names = FALSE))
  train_pool <- setdiff(seq_len(nrow(d)), heldout_pool)
  train_ids <- .rhc_example_sample(train_pool, d$A_label, n_train, seed + 1L)
  test_ids <- .rhc_example_sample(heldout_pool, d$A_label,
                                  n_holdout, seed + 2L)
  stopifnot(!length(intersect(train_ids, test_ids)))
  training <- d[train_ids, , drop = FALSE]
  heldout <- d[test_ids, , drop = FALSE]

  unseen <- vapply(baseline, function(v) {
    if (!is.factor(training[[v]])) return(FALSE)
    any(!as.character(heldout[[v]]) %in%
          unique(as.character(training[[v]])))
  }, logical(1))
  if (any(unseen)) {
    stop("The training subset lacks held-out factor levels for: ",
         paste(baseline[unseen], collapse = ", "),
         ". Try a larger n_train.", call. = FALSE)
  }

  outcomes <- list(
    stay = list(type = "continuous", column = "los_days",
                direction = "lower", probs = c(0.25, 0.5, 0.75),
                grid_length = 9L),
    survival = list(type = "survival", time = "surv_time_180",
                    status = "event180", times = c(30, 90, 150),
                    censoring = "pooled_km", min_G = 0.05)
  )
  cat("Prepared", nrow(d), "patients; fitting", nrow(training),
      "training profiles and predicting for", nrow(heldout), "held-out profiles.\n")
  fit <- drMultiOutITR::fit_dr_itr(
    data = training, treatment = "A_label", covariates = baseline,
    outcomes = outcomes, weights = c(stay = 0.5, survival = 0.5),
    propensity_method = "ranger", nuisance_folds = 2L,
    propensity_eps = 0.05, propensity_truncation = "simplex",
    outcome_learner = "ranger", second_stage = "ranger",
    covariance_source = "oob", num.trees = as.integer(num_trees),
    propensity_trees = as.integer(num_trees), min.node.size = 10L,
    propensity_min_node_size = 50L, survival_outcome_learner = "rsf",
    survival_num.trees = as.integer(num_trees),
    survival_min_node_size = 10L,
    survival_second_stage_round_digits = NULL, seed = as.integer(seed)
  )

  profiles <- heldout[, baseline, drop = FALSE]
  ids <- sprintf("Profile_%03d", seq_len(nrow(profiles)))
  weight_schemes <- list(
    equal = c(stay = 0.5, survival = 0.5),
    LOS_heavy = c(stay = 0.7, survival = 0.3),
    survival_heavy = c(stay = 0.3, survival = 0.7)
  )
  predictions <- lapply(weight_schemes, function(w)
    stats::predict(fit, profiles, ids = ids, weights = w))

  counts <- do.call(rbind, lapply(names(predictions), function(nm) {
    observed <- table(factor(predictions[[nm]]$best_treatment,
                             levels = fit$treatments))
    data.frame(weight_scheme = nm, treatment = names(observed),
               n = as.integer(observed),
               percent = round(100 * as.numeric(observed) / length(ids), 1),
               row.names = NULL)
  }))
  rownames(counts) <- NULL
  cat("\nTreatment recommendations by outcome weighting:\n")
  print(counts, row.names = FALSE)
  cat("\nSource: Vanderbilt Biostatistics RHC data; Connors et al. (1996).\n")
  invisible(list(counts = counts, fit = fit))
}
