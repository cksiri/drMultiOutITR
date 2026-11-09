#' Convert conditional DR estimates into rank-aggregation summary tables
#'
#' Returns the common long-format format accepted by [fit_rank_rule()] and
#' [rank_treatments()]. If an ordinal response was explicitly oriented so that
#' all its component probabilities represent favorable events, choose
#' `larger_is_better=TRUE`. For length of stay and similar original-scale
#' quantiles, choose `FALSE`.
#'
#' @param estimate A `dr_summary_estimate` object, or a list with `theta_long`.
#' @param outcome_id Unique outcome name, e.g. `"survival"`.
#' @param larger_is_better Single logical specifying clinical direction.
#' @param outcome_weight Nonnegative clinical weight on this outcome.
#' @param rho Positive magnitude-discrepancy exponent (paper default 0.5).
#' @param summary_weight Nonnegative within-outcome component weight.
#' @param summary_scale Positive within-component magnitude scale.
#' @return Long-format table with rank-engine columns.
#' @export
as_dr_summary_table <- function(estimate, outcome_id, larger_is_better = TRUE,
                                outcome_weight = 1, rho = 0.5,
                                summary_weight = 1, summary_scale = 1) {
  tl <- estimate$theta_long
  if (!is.data.frame(tl) || !all(c('new_id', 'summary_id', 'treatment', 'theta_hat') %in% names(tl))) {
    stop('estimate must contain $theta_long with standard DR summary columns.')
  }
  if (!is.character(outcome_id) || length(outcome_id) != 1L ||
      is.na(outcome_id) || !nzchar(outcome_id)) stop('outcome_id must be nonempty.')
  if (!is.logical(larger_is_better) || length(larger_is_better) != 1L ||
      is.na(larger_is_better)) stop('larger_is_better must be TRUE or FALSE.')
  if (any(lengths(list(outcome_weight, rho, summary_weight, summary_scale)) != 1L)) {
    stop('Outcome weight, rho, summary_weight and summary_scale must be scalar.')
  }
  out <- data.frame(new_id = as.character(tl$new_id),
                    outcome_id = outcome_id,
                    summary_id = as.character(tl$summary_id),
                    treatment = as.character(tl$treatment),
                    theta_hat = as.numeric(tl$theta_hat),
                    larger_is_better = larger_is_better,
                    outcome_weight = outcome_weight, rho = rho,
                    summary_weight = summary_weight, summary_scale = summary_scale,
                    stringsAsFactors = FALSE)
  validate_matrix_summary_table(out)
}
