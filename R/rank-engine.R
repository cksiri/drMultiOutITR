# Covariance-scaled rank-discordance calculations.

smmr_all_permutations <- function(x) {

  x <- as.vector(x)

  if (length(x) == 1L) {
    return(matrix(x, nrow = 1L))
  }

  out <- lapply(
    seq_along(x),
    function(ii) {

      rest <- smmr_all_permutations(x[-ii])

      cbind(
        rep(x[ii], nrow(rest)),
        rest
      )
    }
  )

  out <- do.call(rbind, out)
  colnames(out) <- paste0("treatment_position_", seq_len(ncol(out)))

  out
}


# Validate matrix-scaling summary table

validate_matrix_summary_table <- function(summary_table) {

  required <- c(
    "new_id",
    "outcome_id",
    "summary_id",
    "treatment",
    "theta_hat",
    "larger_is_better",
    "outcome_weight",
    "rho"
  )

  miss <- setdiff(required, names(summary_table))

  if (length(miss) > 0L) {
    stop(
      "summary_table is missing required columns: ",
      paste(miss, collapse = ", ")
    )
  }

  if (!"summary_weight" %in% names(summary_table)) {
    summary_table$summary_weight <- 1
  }

  if (!"summary_scale" %in% names(summary_table)) {
    summary_table$summary_scale <- 1
  }

  if (nrow(summary_table) == 0L) {
    stop("summary_table has zero rows.")
  }

  if (anyNA(summary_table$new_id)) {
    stop("new_id contains missing values.")
  }

  if (anyNA(summary_table$outcome_id)) {
    stop("outcome_id contains missing values.")
  }

  if (anyNA(summary_table$summary_id)) {
    stop("summary_id contains missing values.")
  }

  if (anyNA(summary_table$treatment)) {
    stop("treatment contains missing values.")
  }

  if (any(!is.finite(summary_table$theta_hat))) {
    stop("Non-finite theta_hat values found.")
  }

  if (any(!is.finite(summary_table$rho)) ||
      any(summary_table$rho <= 0)) {
    stop("rho must contain positive finite values.")
  }

  if (any(!is.finite(summary_table$outcome_weight)) ||
      any(summary_table$outcome_weight < 0)) {
    stop("outcome_weight must contain non-negative finite values.")
  }

  if (any(!is.finite(summary_table$summary_weight)) ||
      any(summary_table$summary_weight < 0)) {
    stop("summary_weight must contain non-negative finite values.")
  }

  if (any(!is.finite(summary_table$summary_scale)) ||
      any(summary_table$summary_scale <= 0)) {
    stop("summary_scale must contain positive finite values.")
  }

  summary_table
}


# Extract consistent treatment ordering

get_matrix_treatment_order <- function(summary_table) {

  trt <- summary_table$treatment

  if (is.factor(trt)) {

    lev <- levels(trt)
    lev <- lev[lev %in% as.character(trt)]

    if (length(lev) > 0L) {
      return(as.character(lev))
    }
  }

  sort(unique(as.character(trt)))
}


# Extract outcome ordering and weights

get_matrix_outcome_info <- function(summary_table) {

  outcomes <- unique(as.character(summary_table$outcome_id))

  weights <- vapply(
    outcomes,
    function(jj) {

      ww <- unique(
        summary_table$outcome_weight[
          as.character(summary_table$outcome_id) == jj
        ]
      )

      ww <- ww[is.finite(ww)]

      if (length(ww) != 1L) {
        stop(
          "Outcome ", jj,
          " must have exactly one outcome_weight."
        )
      }

      ww
    },
    numeric(1)
  )

  if (any(weights < 0)) {
    stop("Outcome weights must be non-negative.")
  }

  if (all(weights == 0)) {
    stop("At least one outcome weight must be positive.")
  }

  list(
    outcomes = outcomes,
    weights = weights
  )
}


# Internal structural checks

check_matrix_summary_structure <- function(summary_table,
                                           treatment_order = NULL,
                                           outcome_order = NULL) {

  summary_table <- validate_matrix_summary_table(summary_table)

  if (is.null(treatment_order)) {
    treatment_order <- get_matrix_treatment_order(summary_table)
  }

  if (is.null(outcome_order)) {
    outcome_order <- get_matrix_outcome_info(summary_table)$outcomes
  }

  ids <- unique(summary_table$new_id)

  for (id_i in ids) {

    dat_i <- summary_table[
      summary_table$new_id == id_i,
      ,
      drop = FALSE
    ]

    for (outcome_j in outcome_order) {

      dat_j <- dat_i[
        as.character(dat_i$outcome_id) == outcome_j,
        ,
        drop = FALSE
      ]

      if (nrow(dat_j) == 0L) {
        stop(
          "No rows found for new_id=", id_i,
          ", outcome_id=", outcome_j, "."
        )
      }

      summary_ids <- unique(as.character(dat_j$summary_id))

      for (summary_l in summary_ids) {

        tmp <- dat_j[
          as.character(dat_j$summary_id) == summary_l,
          ,
          drop = FALSE
        ]

        present <- as.character(tmp$treatment)

        if (!setequal(present, treatment_order)) {
          stop(
            "Treatment structure is incomplete for new_id=", id_i,
            ", outcome_id=", outcome_j,
            ", summary_id=", summary_l, "."
          )
        }

        if (anyDuplicated(present)) {
          stop(
            "Duplicate treatment rows found for new_id=", id_i,
            ", outcome_id=", outcome_j,
            ", summary_id=", summary_l, "."
          )
        }
      }
    }
  }

  invisible(TRUE)
}


# Compute Pi(delta, x) for a single patient

compute_pi_one_patient <- function(summary_table_one,
                                   delta,
                                   treatment_order,
                                   outcome_order,
                                   magnitude_scale = c(
                                     "none",
                                     "summary_scale"
                                   ),
                                   ties.method = "average") {

  magnitude_scale <- match.arg(magnitude_scale)

  K <- length(treatment_order)

  if (length(delta) != K) {
    stop("Length of delta must equal number of treatments.")
  }

  if (!setequal(delta, seq_len(K))) {
    stop("delta must be a permutation of 1,...,K.")
  }

  pi_vec <- numeric(length(outcome_order))
  names(pi_vec) <- outcome_order

  for (jj in seq_along(outcome_order)) {

    outcome_j <- outcome_order[jj]

    dat_j <- summary_table_one[
      as.character(summary_table_one$outcome_id) == outcome_j,
      ,
      drop = FALSE
    ]

    if (nrow(dat_j) == 0L) {
      stop("Missing outcome ", outcome_j, " for one evaluation profile.")
    }

    summary_ids <- unique(as.character(dat_j$summary_id))

    pi_j <- 0

    for (ll in summary_ids) {

      tmp <- dat_j[
        as.character(dat_j$summary_id) == ll,
        ,
        drop = FALSE
      ]

      idx <- match(
        treatment_order,
        as.character(tmp$treatment)
      )

      if (anyNA(idx)) {
        stop(
          "Not every treatment is represented for outcome ",
          outcome_j, ", summary ", ll, "."
        )
      }

      tmp <- tmp[idx, , drop = FALSE]

      theta <- as.numeric(tmp$theta_hat)

      direction <- unique(tmp$larger_is_better)

      if (length(direction) != 1L || is.na(direction)) {
        stop(
          "larger_is_better must be one non-missing value ",
          "within each summary."
        )
      }

      rho_val <- unique(tmp$rho)

      if (length(rho_val) != 1L ||
          !is.finite(rho_val) ||
          rho_val <= 0) {
        stop("rho must be one positive finite value within each summary.")
      }

      sw <- unique(tmp$summary_weight)

      if (length(sw) != 1L ||
          !is.finite(sw) ||
          sw < 0) {
        stop(
          "summary_weight must be one non-negative finite value ",
          "within each summary."
        )
      }

      ss <- unique(tmp$summary_scale)

      if (length(ss) != 1L ||
          !is.finite(ss) ||
          ss <= 0) {
        stop(
          "summary_scale must be one positive finite value ",
          "within each summary."
        )
      }

      # Rank treatments according to the current outcome summary.

      if (isTRUE(direction)) {
        r_obs <- rank(
          -theta,
          ties.method = ties.method
        )
      } else {
        r_obs <- rank(
          theta,
          ties.method = ties.method
        )
      }

      # delta[k] is the rank assigned to treatment k, not a treatment index.
      # Compare theta[k] with the outcome value at that rank position.
      # Using theta[delta] instead would give a degenerate identity ranking.

      if (isTRUE(direction)) {
        theta_ranked <- sort(theta, decreasing = TRUE)
      } else {
        theta_ranked <- sort(theta, decreasing = FALSE)
      }

      theta_candidate <- theta_ranked[delta]

      abs_diff <- abs(theta - theta_candidate)

      if (magnitude_scale == "summary_scale") {
        abs_diff <- abs_diff / ss
      }

      Phi <- abs_diff^rho_val

      # Rank discrepancy.

      gamma <- abs(r_obs - delta)

      # Component contribution:

      component_distance <- sum(gamma * Phi)

      pi_j <- pi_j + sw * component_distance
    }

    pi_vec[jj] <- pi_j
  }

  pi_vec
}


# Compute Pi(delta, X_i) for every patient

compute_pi_population <- function(summary_table,
                                  delta,
                                  treatment_order = NULL,
                                  outcome_order = NULL,
                                  magnitude_scale = c(
                                    "none",
                                    "summary_scale"
                                  ),
                                  ties.method = "average") {

  magnitude_scale <- match.arg(magnitude_scale)

  summary_table <- validate_matrix_summary_table(summary_table)

  if (is.null(treatment_order)) {
    treatment_order <- get_matrix_treatment_order(summary_table)
  }

  if (is.null(outcome_order)) {
    outcome_order <- get_matrix_outcome_info(
      summary_table
    )$outcomes
  }

  ids <- unique(summary_table$new_id)

  Pi <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(outcome_order)
  )

  rownames(Pi) <- as.character(ids)
  colnames(Pi) <- outcome_order

  for (ii in seq_along(ids)) {

    tmp <- summary_table[
      summary_table$new_id == ids[ii],
      ,
      drop = FALSE
    ]

    Pi[ii, ] <- compute_pi_one_patient(
      summary_table_one = tmp,
      delta = delta,
      treatment_order = treatment_order,
      outcome_order = outcome_order,
      magnitude_scale = magnitude_scale,
      ties.method = ties.method
    )
  }

  Pi
}


# Regularize covariance matrix

regularize_smmr_covariance <- function(Sigma,
                                       ridge = 1e-6,
                                       eigen_floor = 1e-8) {

  Sigma <- as.matrix(Sigma)

  if (nrow(Sigma) != ncol(Sigma)) {
    stop("Sigma must be square.")
  }

  J <- nrow(Sigma)

  if (J < 1L) {
    stop("Sigma must have positive dimension.")
  }

  if (any(!is.finite(Sigma))) {
    stop("Sigma contains non-finite values.")
  }

  if (!is.finite(ridge) || ridge < 0) {
    stop("ridge must be a non-negative finite value.")
  }

  if (!is.finite(eigen_floor) || eigen_floor <= 0) {
    stop("eigen_floor must be a positive finite value.")
  }

  ## Force exact symmetry.
  Sigma <- (Sigma + t(Sigma)) / 2

  ## Preserve the truly unregularized eigenvalues for diagnostics.
  eig_raw <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  diag_mean <- mean(diag(Sigma))

  if (!is.finite(diag_mean) || diag_mean <= 0) {

    finite_vals <- abs(Sigma[is.finite(Sigma)])

    if (length(finite_vals) == 0L ||
        max(finite_vals) <= 0) {
      diag_mean <- 1
    } else {
      diag_mean <- max(finite_vals)
    }
  }

  # Ridge regularization.

  Sigma_ridge <- Sigma +
    diag(ridge * diag_mean, J)

  # Minimum covariance eigenvalue.

  eg <- eigen(
    Sigma_ridge,
    symmetric = TRUE
  )

  max_eig <- max(abs(eg$values))

  floor_val <- max(
    eigen_floor * max(max_eig, 1),
    .Machine$double.eps
  )

  eig_regularized <- pmax(
    eg$values,
    floor_val
  )

  Sigma_reg <- eg$vectors %*%
    diag(eig_regularized, nrow = J, ncol = J) %*%
    t(eg$vectors)

  Sigma_reg <- (Sigma_reg + t(Sigma_reg)) / 2

  condition_number <- max(eig_regularized) /
    min(eig_regularized)

  list(
    raw = Sigma,
    regularized = Sigma_reg,
    eigenvalues_raw = eig_raw,
    eigenvalues_after_ridge = eg$values,
    eigenvalues_regularized = eig_regularized,
    condition_number = condition_number,
    ridge = ridge,
    eigen_floor = eigen_floor
  )
}


# Estimate Sigma(delta) for every candidate ranking

estimate_smmr_scaling_matrices <- function(
    scaling_summary_table,
    magnitude_scale = c(
      "none",
      "summary_scale"
    ),
    ridge = 1e-6,
    eigen_floor = 1e-8,
    ties.method = "average",
    check_structure = TRUE) {

  magnitude_scale <- match.arg(magnitude_scale)

  scaling_summary_table <- validate_matrix_summary_table(
    scaling_summary_table
  )

  treatment_order <- get_matrix_treatment_order(
    scaling_summary_table
  )

  outcome_info <- get_matrix_outcome_info(
    scaling_summary_table
  )

  outcome_order <- outcome_info$outcomes

  if (isTRUE(check_structure)) {
    check_matrix_summary_structure(
      summary_table = scaling_summary_table,
      treatment_order = treatment_order,
      outcome_order = outcome_order
    )
  }

  K <- length(treatment_order)

  candidates <- smmr_all_permutations(
    seq_len(K)
  )

  out <- vector(
    "list",
    nrow(candidates)
  )

  candidate_labels <- apply(
    candidates,
    1,
    paste,
    collapse = "-"
  )

  names(out) <- candidate_labels

  for (dd in seq_len(nrow(candidates))) {

    delta <- as.integer(candidates[dd, ])

    Pi_train <- compute_pi_population(
      summary_table = scaling_summary_table,
      delta = delta,
      treatment_order = treatment_order,
      outcome_order = outcome_order,
      magnitude_scale = magnitude_scale,
      ties.method = ties.method
    )

    if (nrow(Pi_train) < 2L) {
      stop(
        "At least two training profiles are required ",
        "to estimate Sigma(delta)."
      )
    }

    Sigma_raw <- stats::cov(Pi_train)

    ## stats::cov() can return a scalar for one outcome.
    if (length(outcome_order) == 1L) {
      Sigma_raw <- matrix(
        Sigma_raw,
        nrow = 1L,
        ncol = 1L,
        dimnames = list(outcome_order, outcome_order)
      )
    }

    reg <- regularize_smmr_covariance(
      Sigma = Sigma_raw,
      ridge = ridge,
      eigen_floor = eigen_floor
    )

    out[[dd]] <- list(
      delta = delta,
      candidate_label = candidate_labels[dd],
      Pi_train = Pi_train,
      mean_Pi = colMeans(Pi_train),
      Sigma_raw = reg$raw,
      Sigma = reg$regularized,
      eigenvalues_raw = reg$eigenvalues_raw,
      eigenvalues_regularized = reg$eigenvalues_regularized,
      condition_number = reg$condition_number
    )
  }

  list(
    candidates = candidates,
    treatment_order = treatment_order,
    outcome_order = outcome_order,
    outcome_weights = outcome_info$weights,
    scaling = out,
    magnitude_scale = magnitude_scale,
    ridge = ridge,
    eigen_floor = eigen_floor
  )
}


# Matrix-scaled rank decision

select_rank_based_treatment_matrix <- function(
    summary_table,
    scaling_summary_table,
    magnitude_scale = c(
      "none",
      "summary_scale"
    ),
    ridge = 1e-6,
    eigen_floor = 1e-8,
    ties.method = "average",
    check_structure = TRUE,
    return_details = TRUE) {

  magnitude_scale <- match.arg(magnitude_scale)

  summary_table <- validate_matrix_summary_table(
    summary_table
  )

  scaling_summary_table <- validate_matrix_summary_table(
    scaling_summary_table
  )

  # Candidate-specific covariance matrices.

  scale_fit <- estimate_smmr_scaling_matrices(
    scaling_summary_table = scaling_summary_table,
    magnitude_scale = magnitude_scale,
    ridge = ridge,
    eigen_floor = eigen_floor,
    ties.method = ties.method,
    check_structure = check_structure
  )

  candidates <- scale_fit$candidates
  treatment_order <- scale_fit$treatment_order
  outcome_order <- scale_fit$outcome_order

  # Validate test-table structure against training.

  test_treatment_order <- get_matrix_treatment_order(
    summary_table
  )

  if (!setequal(
    treatment_order,
    test_treatment_order
  )) {
    stop(
      "Treatment sets differ between summary_table ",
      "and scaling_summary_table."
    )
  }

  test_outcome_info <- get_matrix_outcome_info(
    summary_table
  )

  if (!setequal(
    outcome_order,
    test_outcome_info$outcomes
  )) {
    stop(
      "Outcome sets differ between summary_table ",
      "and scaling_summary_table."
    )
  }

  if (isTRUE(check_structure)) {
    check_matrix_summary_structure(
      summary_table = summary_table,
      treatment_order = treatment_order,
      outcome_order = outcome_order
    )
  }

  # Use requested outcome weights.

  lambda <- vapply(
    outcome_order,
    function(jj) {
      test_outcome_info$weights[
        match(jj, test_outcome_info$outcomes)
      ]
    },
    numeric(1)
  )

  Lambda_half <- diag(
    sqrt(lambda),
    nrow = length(lambda),
    ncol = length(lambda)
  )

  rownames(Lambda_half) <- outcome_order
  colnames(Lambda_half) <- outcome_order

  # Candidate rankings for new profiles.

  ids <- unique(summary_table$new_id)

  detail_list <- vector(
    "list",
    length(ids) * nrow(candidates)
  )

  row_counter <- 1L

  selection_list <- vector(
    "list",
    length(ids)
  )

  for (ii in seq_along(ids)) {

    id_i <- ids[ii]

    tmp_i <- summary_table[
      summary_table$new_id == id_i,
      ,
      drop = FALSE
    ]

    W_vec <- numeric(
      nrow(candidates)
    )

    for (dd in seq_len(nrow(candidates))) {

      delta <- as.integer(
        candidates[dd, ]
      )

      Pi_test <- compute_pi_one_patient(
        summary_table_one = tmp_i,
        delta = delta,
        treatment_order = treatment_order,
        outcome_order = outcome_order,
        magnitude_scale = magnitude_scale,
        ties.method = ties.method
      )

      Sigma_delta <- scale_fit$scaling[[dd]]$Sigma

      weighted_Pi <- as.vector(
        Lambda_half %*% Pi_test
      )

      invSigma_times_Pi <- tryCatch(
        solve(
          Sigma_delta,
          weighted_Pi
        ),
        error = function(e) {
          qr.solve(
            Sigma_delta,
            weighted_Pi
          )
        }
      )

      W_val <- sum(
        weighted_Pi * invSigma_times_Pi
      )

      if (!is.finite(W_val)) {
        W_val <- Inf
      }

      W_vec[dd] <- W_val

      detail_row <- data.frame(
        new_id = id_i,
        candidate_id = dd,
        candidate_ranking = paste(delta, collapse = "-"),
        W = W_val,
        condition_number = scale_fit$scaling[[dd]]$condition_number,
        stringsAsFactors = FALSE
      )

      for (jj in seq_along(outcome_order)) {

        detail_row[[
          paste0(
            "pi_",
            outcome_order[jj]
          )
        ]] <- Pi_test[jj]
      }

      detail_list[[row_counter]] <- detail_row

      row_counter <- row_counter + 1L
    }

    # Select the lowest-discrepancy candidate ranking.

    if (all(!is.finite(W_vec))) {
      stop(
        "All candidate ranking distances are non-finite ",
        "for new_id=", id_i, "."
      )
    }

    best_dd <- which.min(W_vec)

    best_delta <- as.integer(
      candidates[best_dd, ]
    )

    best_treatment <- treatment_order[
      which.min(best_delta)
    ]

    selection_list[[ii]] <- data.frame(
      new_id = id_i,
      best_treatment = best_treatment,
      optimal_ranking = paste(best_delta, collapse = "-"),
      min_W = W_vec[best_dd],
      candidate_id = best_dd,
      stringsAsFactors = FALSE
    )
  }

  selection <- do.call(
    rbind,
    selection_list
  )

  details <- do.call(
    rbind,
    detail_list
  )

  rownames(selection) <- NULL
  rownames(details) <- NULL

  out <- list(
    selection = selection,
    scaling_fit = scale_fit,
    settings = list(
      magnitude_scale = magnitude_scale,
      ridge = ridge,
      eigen_floor = eigen_floor,
      ties.method = ties.method
    )
  )

  if (isTRUE(return_details)) {
    out$details <- details
  }

  out
}


