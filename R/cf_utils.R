#' Compute normalized weights and conformal threshold eta(x)
#'
#' Internal helper function used in weighted split-CQR to normalize calibration
#' and test weights and compute the conformal threshold \eqn{\eta(x)}.
#'
#' @param scores Numeric matrix of conformal scores from the calibration set (output of `compute_cqr_scores`) for each level in `alphas`.
#' @param weights_cal Numeric vector of calibration weights.
#' @param weights_test Numeric vector, weight associated with the test point.
#' @param alphas Numeric vector of miscoverage levels in (0,1).
#'
#' @return A numeric matrix giving the conformal threshold \eqn{\eta(x)} with dimension `nrow(data_cal) x length(alphas)`.
eta_conformeR <-function(scores,weights_cal,weights_test,alphas){ # version 3.0 vectorize over alpha

  # Weight's normalization
  denom <- sum(weights_cal)+weights_test # \perp alpha
  # Normalized calibration weights
  p_hat <- weights_cal%*%t(1/denom)# nrow(cal) x nrow(test) # \perp alpha

  # Weight at infinity
  p_hat_infty <- weights_test/denom # \perp alpha

  # Work col-wise across alpha
  ord <- apply(scores, 2, order) # nrow(cal) x length(alphas)
  score_s <- apply(scores, 2, sort)

  # Reorder weights (p_hat) according to score order for each alpha
  p_hat_s <- lapply(seq_along(alphas), function(a) p_hat[ord[, a], , drop = FALSE])  # list of (n_cal x n_test) per alpha

  # Cumulative sum of weights
  cum_p_list <- lapply(p_hat_s, function(p) apply(p, 2, cumsum))  # list of (n_cal x n_test) matrices

  # Compute quantiles per alpha
  eta_list <- lapply(seq_along(alphas), function(a) {
    idx_first <- apply(cum_p_list[[a]], 2, function(col) which(col >= 1 - alphas[a])[1])

    # Fallback to weight at infinity if threshold not reached
    eta <- ifelse(!is.na(idx_first), score_s[idx_first,a], p_hat_infty)
    return(eta)
  })

  # Return as matrix: (n_test x n_alpha)
  do.call(cbind, eta_list)
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
#' @param gene Character string giving the name of the target gene.
#' @param alphas Numeric vector of miscoverage levels in (0,1).
#'
#' @return A numeric matrix of conformal scores with dimension `nrow(data_cal) x length(alphas)`.
compute_cqr_scores <- function(qr_model, data_cal, gene, alphas) { # vectorized version

  # Predict lower and upper quantiles
  q_lo  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = alphas/2)$predictions
  q_hi  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = 1 - alphas/2)$predictions # matrix of dim nrow(cal_set) x length(alphas)

  # Compute conformal scores
  scores <- pmax(q_lo - data_cal[[gene]], data_cal[[gene]] - q_hi) # returns a vector of dimension nrow(cal_set)*length(alphas)

  return(matrix(scores, nrow=nrow(data_cal), byrow=F))
}


#' Build conformal prediction intervals for difference between observed and counterfactual logcounts
#'
#' Constructs prediction intervals for the difference between observed and
#' counterfactual logcounts using weighted split-CQR. Vectorized: accepts
#' multiple test rows at once.
#'
#' @param test_data A data frame of test cells (rows = cells, cols = covariates).
#' @param qr_model A quantile regression model fitted with conformalCQR.
#' @param scores Numeric matrix of conformal scores with dimension `nrow(data_cal) x length(alphas)`.
#' @param weights_cal Numeric vector of calibration weights.
#' @param w_test Numeric vector of weights associated with the test samples.
#' @param alphas Numeric vector of miscoverage levels in (0,1).
#' @param gene Character string of the outcome gene.
#' @param gene_names Character vector of all gene names.
#' @param tr_flag Integer-vector flag (0 or 1) indicating treatment/control setting for each row of `test_data`
#'
#' @return A matrix with `nrow(test_data)*length(alphas)` with three columns ("lower","upper","cell_id","alpha").
build_intervals <- function(test_data, qr_model, scores, weights_cal,
                            w_test, alphas, gene, gene_names, tr_flag) {

  tr_flag <- rep(tr_flag, nrow(test_data)) # tr_flag = 0 or 1
  eta_mat <- eta_conformeR(scores, weights_cal, w_test, alphas) # returns mat of dim nrow(test) x length(alphas)
  # predict lower/upper quantiles for all test rows
  q_lo <- predict(qr_model, data = test_data, type = "quantiles",
                  quantiles = alphas/2)$predictions
  q_hi <- predict(qr_model, data = test_data, type = "quantiles",
                  quantiles = 1 - alphas/2)$predictions

  y_obs <- test_data[[gene]]
  lower <- tr_flag*(y_obs - q_hi - eta_mat) + (1-tr_flag)*(q_lo - eta_mat - y_obs) # mat of dim nrow(test) x length(alphas)
  upper <- tr_flag*(y_obs - q_lo + eta_mat) + (1-tr_flag)*(q_hi + eta_mat - y_obs) # mat of dim nrow(test) x length(alphas)
  int <- cbind.data.frame(lower = c(lower), upper =  c(upper), cell_id = rep(paste(tr_flag,seq_len(nrow(test_data)),sep="."),length(alphas)), alpha=rep(alphas, each=nrow(test_data)))
  rownames(int) <- NULL
  int
}
