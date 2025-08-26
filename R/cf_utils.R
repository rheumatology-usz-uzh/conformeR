#' Compute normalized weights and conformal threshold eta(x)
#'
#' Internal helper function used in weighted split-CQR to normalize calibration
#' and test weights and compute the conformal threshold \eqn{\eta(x)}.
#'
#' @param scores Numeric vector of conformal scores from the calibration set (output of `compute_cqr_scores`).
#' @param weights_cal Numeric vector of calibration weights.
#' @param weights_test Numeric scalar, weight associated with the test point.
#' @param alpha Numeric scalar, miscoverage level in (0,1).
#'
#' @return A numeric scalar giving the conformal threshold \eqn{\eta(x)}.
eta_conformeR <-function(scores,weights_cal,weights_test,alpha=alpha){

  denom <- sum(weights_cal)+weights_test
  p_hat <- weights_cal/denom
  p_hat_infty <- weights_test/denom

  ord <- order(scores)
  scores_s <- scores[ord]
  p_s <- p_hat[ord]
  cum_p <- cumsum(p_s)
  target <- 1 - alpha

  if (cum_p[length(cum_p)] >= target) {
    idx <- which(cum_p >= target)[1]
    eta <- scores_s[idx]
  } else {
    eta <- p_hat_infty
  }
  return(eta)
}

#' Compute conformal CQR scores
#'
#' Internal helper function that computes calibration scores for prediction
#' intervals in counterfactual CQR models. Scores are defined as the maximum
#' deviation of the observed outcome from the predicted lower and upper
#' quantiles.
#'
#' @param qr_model A quantile regression model fitted with conformalCQR.
#' @param data_cal A data frame containing calibration covariates and outcomes.
#' @param gene_idx Column index of the outcome gene.
#' @param alpha Numeric scalar, miscoverage level in (0,1).
#'
#' @return A numeric vector of conformal scores based on some calibration set.
compute_cqr_scores <- function(qr_model, data_cal, gene_idx, alpha) {
  # Predict lower and upper quantiles
  q_lo  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = alpha/2)$predictions
  q_hi  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = 1 - alpha/2)$predictions

  # Compute conformal scores
  scores <- pmax(q_lo - data_cal[[gene_idx]], data_cal[[gene_idx]] - q_hi)

  return(as.vector(scores))
}

#' Build conformal prediction intervals for difference between observed and counterfactual logcounts
#'
#' Internal helper function that constructs prediction intervals for the
#' difference between observed and counterfactual logcounts using weighted
#' split-CQR.
#'
#' @param test_row A single-row data frame containing covariates (all-but-gene-of-interest logcounts) for the test cell.
#' @param qr_model A quantile regression model fitted with conformalCQR.
#' @param other_scores Numeric vector of calibration scores.
#' @param other_weights Numeric vector of calibration weights.
#' @param w_x Numeric scalar, weight associated with the test sample.
#' @param alpha Numeric scalar, miscoverage level in (0,1).
#' @param gene_idx Column index of the outcome gene.
#' @param gene_names Character vector of all gene names.
#' @param tr_flag Integer flag (0 or 1) indicating treatment/control setting.
#'
#' @return A numeric vector of length two, giving the lower and upper bounds of the conformal prediction interval.
build_intervals <- function(test_row, qr_model, other_scores, other_weights, w_x, alpha, gene_idx, gene_names, tr_flag) {
  if (tr_flag == 0) {
    eta_x <- eta_conformeR(other_scores, other_weights, w_x, alpha)
    q_lo <- predict(qr_model, data = test_row, type = "quantiles", quantiles = alpha/2)$predictions
    q_hi <- predict(qr_model, data = test_row, type = "quantiles", quantiles = 1 - alpha/2)$predictions
    return(c(q_lo - eta_x - test_row[, gene_idx], q_hi + eta_x - test_row[, gene_idx]))
  } else {
    eta_x <- eta_conformeR(other_scores, other_weights, w_x, alpha)
    q_lo <- predict(qr_model, data = test_row, type = "quantiles", quantiles = alpha/2)$predictions
    q_hi <- predict(qr_model, data = test_row, type = "quantiles", quantiles = 1 - alpha/2)$predictions
    return(c(test_row[, gene_idx] - q_hi - eta_x, test_row[, gene_idx] - q_lo + eta_x))
  }
}
