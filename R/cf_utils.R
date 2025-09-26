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
#' @param gene Character string giving the name of the target gene.
#' @param alpha Numeric scalar, miscoverage level in (0,1).
#'
#' @return A numeric vector of conformal scores based on some calibration set.
compute_cqr_scores <- function(qr_model, data_cal, gene, alpha) {
  # Predict lower and upper quantiles
  q_lo  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = alpha/2)$predictions
  q_hi  <- predict(qr_model, data = data_cal, type = "quantiles", quantiles = 1 - alpha/2)$predictions

  # Compute conformal scores
  scores <- pmax(q_lo - data_cal[[gene]], data_cal[[gene]] - q_hi)

  return(as.vector(scores))
}


#' Build conformal prediction intervals for difference between observed and counterfactual logcounts
#'
#' Constructs prediction intervals for the difference between observed and
#' counterfactual logcounts using weighted split-CQR. Vectorized: accepts
#' multiple test rows at once.
#'
#' @param test_data A data frame of test cells (rows = cells, cols = covariates).
#' @param qr_model A quantile regression model fitted with conformalCQR.
#' @param other_scores Numeric vector of calibration scores.
#' @param other_weights Numeric vector of calibration weights.
#' @param w_x Numeric vector of weights associated with the test samples.
#' @param alpha Numeric scalar, miscoverage level in (0,1).
#' @param gene_idx Column index of the outcome gene.
#' @param gene_names Character vector of all gene names.
#' @param tr_flag Integer-vector flag (0 or 1) indicating treatment/control setting for each row of `test_data`
#'
#' @return A numeric matrix with two columns ("lower","upper"), one row per test cell.
build_intervals <- function(test_data, qr_model, other_scores, other_weights,
                            w_x, alpha, gene, gene_names, tr_flag) {
  n <- nrow(test_data)

  # conformal correction per test row
  eta_x <- vapply(w_x, function(w) eta_conformeR(other_scores, other_weights, w, alpha),
                  numeric(1L))

  # predict lower/upper quantiles for all test rows
  q_lo <- predict(qr_model, data = test_data, type = "quantiles",
                  quantiles = alpha/2)$predictions
  q_hi <- predict(qr_model, data = test_data, type = "quantiles",
                  quantiles = 1 - alpha/2)$predictions

  y_obs <- test_data[[gene]]

  lower <- ifelse(tr_flag == 0,
                  q_lo - eta_x - y_obs,
                  y_obs - q_hi - eta_x)

  upper <- ifelse(tr_flag == 0,
                  q_hi + eta_x - y_obs,
                  y_obs - q_lo + eta_x)
  out <- cbind(lower = lower, upper = upper)
  rownames(out) <- NULL
  out
}

#' Process a single group to compute conformal prediction intervals
#'
#' This internal helper function computes conformal prediction intervals for a specific
#' conformal group. It splits further the data into training,
#' calibration, and test subsets for the group, fits propensity score and quantile
#' regression models for each gene, and constructs conformal prediction intervals
#' across multiple significance levels (alphas).
#' @importFrom dplyr filter
#'
#' @param g Character or factor. The identifier of the conformal group to process.
#' @param dt_train A \code{data.frame} containing the training set, including covariates,
#'   outcome genes, and the \code{conf_group} column.
#' @param dt_cal A \code{data.frame} containing the calibration set.
#' @param dt_test A \code{data.frame} containing the test set.
#' @param gene_names Character vector of all gene names in the dataset.
#' @param alphas Numeric vector of miscoverage levels for which conformal intervals
#'   are computed.
#' @param obs_condition Character. Name of the column representing the observed condition
#'   (treatment vs. control).
#' @param gene_chunk Integer. Number of genes to process per chunk in order to manage memory usage.
#' @param param3 A \code{BiocParallelParam} object controlling parallelization for per-gene computations.
#'
#' @return A \code{data.frame} containing conformal prediction intervals for all genes in the
#'   group. Columns include:
#'   \itemize{
#'     \item \code{lower} - Lower bound of the conformal interval
#'     \item \code{upper} - Upper bound of the conformal interval
#'     \item \code{alpha} - Miscoverage level corresponding to the interval
#'     \item \code{gene} - Gene name
#'     \item \code{covered} - Factor indicating whether 0 falls inside the interval ("inside"/"outside")
#'     \item \code{conf_group} - Identifier of the processed group
#'   }
#'
#' @details
#' This function manages the group-level workflow by:
#' \enumerate{
#'   \item Subsetting the training, calibration, and test sets by conformal group.
#'   \item Preparing combined training/calibration data frames for propensity score modeling.
#'   \item Splitting the list of genes into chunks and delegating computation to
#'         \code{process_gene_chunk()}.
#'   \item Combining per-chunk results into a single data.frame for the group.
#' }
#'
#' @seealso \code{\link{process_gene_chunk}}

process_group <- function(g, dt_train, dt_cal, dt_test, gene_names, alphas, obs_condition, gene_chunk,param3) {
  # Subset by group & condition
  train_g <- dt_train |> filter(conf_group == g)
  cal_g   <- dt_cal |> filter(conf_group == g)
  test_g  <- dt_test |> filter(conf_group == g)


  # combined training set (converted once to data.frame)
  train_combined_df <- rbind(train_g[train_g[[obs_condition]] == 0,],
                             train_g[train_g[[obs_condition]] == 1,])
  # calibration
  cal_combined_df <- rbind(cal_g[cal_g[[obs_condition]] == 0,],
                           cal_g[cal_g[[obs_condition]] == 1,])

  all_gene_indices <- seq_along(gene_names)
  chunk_starts <- seq(1, length(all_gene_indices), by = gene_chunk)
  group_gene_outs <- lapply(chunk_starts, function(s) {
    chunk_idxs <- all_gene_indices[s:min(s + gene_chunk - 1, length(all_gene_indices))]
    process_gene_chunk(chunk_idxs,train_combined_df, cal_combined_df, test_g, gene_names, alphas, obs_condition, param3)
  })

  group_df <- do.call(rbind, group_gene_outs)
  group_df$conf_group <- rep(g,nrow(group_df))
  group_df
}


#' Process a chunk of genes to compute conformal prediction intervals
#'
#' This internal helper function computes conformal prediction intervals for a subset of
#' genes (a "chunk") within a given conformal group. For each gene, it estimates
#' propensity scores, fits quantile regression models for treatment and control groups,
#' computes calibration scores, and constructs conformal prediction intervals across
#' multiple alphas.
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr mutate filter pull
#'
#' @param some_idxs Integer vector. Indices of the genes to process in this chunk.
#' @param train_combined_df A \code{data.frame}. Combined training set (treatment and control) for
#'   propensity score modeling.
#' @param cal_combined_df A \code{data.frame}. Combined calibration set (treatment and control).
#' @param test_g A \code{data.frame} containing the test set for this group.
#' @param gene_names Character vector of all gene names in the dataset.
#' @param alphas Numeric vector of miscoverage levels for which conformal intervals
#'   are computed.
#' @param obs_condition Character. Name of the column representing the observed condition
#'   (treatment vs. control).
#' @param param3 A \code{BiocParallelParam} object controlling parallelization across genes.
#'
#' @return A \code{data.frame} with conformal prediction intervals for all genes in the chunk.
#'   The returned data.frame contains the following columns:
#'   \itemize{
#'     \item \code{lower} - Lower bound of the conformal interval
#'     \item \code{upper} - Upper bound of the conformal interval
#'     \item \code{alpha} - Miscoverage level corresponding to the interval
#'     \item \code{gene} - Gene name
#'     \item \code{covered} - Factor indicating whether 0 falls inside the interval ("inside"/"outside")
#'   }
#'
#' @details
#' For each gene in the chunk, the function:
#' \enumerate{
#'   \item Fits a propensity score model using the combined training set.
#'   \item Computes calibration weights for treatment and control samples.
#'   \item Fits quantile regression models separately for treatment and control samples.
#'   \item Computes calibration scores for all specified alphas.
#'   \item Constructs conformal prediction intervals for all test cells, taking into account
#'         their observed treatment/control condition.
#' }
#'
#' Parallelization over genes is handled via \code{bplapply}. Computation within test sets
#' is vectorized to reduce memory usage.
#'
#' @seealso \code{\link{process_group}}

process_gene_chunk <- function( some_idxs, train_combined_df, cal_combined_df, test_g, gene_names, alphas, obs_condition, param3) {
  gene_results <- bplapply(some_idxs, function(gene_idx) {
    gene <- gene_names[gene_idx]
    # propensity model
    ps_model <- prop_score(train_combined_df, gene, gene_names, obs_condition)
    ps_cal_full <- predict(ps_model, cal_combined_df, type = "prob") |> dplyr::pull(.pred_1)
    w_cal_full <- (1 - ps_cal_full) / ps_cal_full
    wC0 <- w_cal_full[1:nrow(cal_combined_df[cal_combined_df[[obs_condition]] == 0,])]
    wC1 <- w_cal_full[(nrow(cal_combined_df[cal_combined_df[[obs_condition]] == 0,]) + 1):length(w_cal_full)]
    # train QR
    qrT0 <- train_qr(train_combined_df[train_combined_df[[obs_condition]] == 0,], gene, gene_names)
    qrT1 <- train_qr(train_combined_df[train_combined_df[[obs_condition]] == 1,], gene, gene_names)

    # test props
    ps_test <- predict(ps_model, test_g, type = "prob") |> dplyr::pull(.pred_1)
    w_test  <- (1 - ps_test) / ps_test

    # calibration scores
    scoresT0_list <- lapply(alphas, function(a) compute_cqr_scores(qrT0, cal_combined_df[cal_combined_df[[obs_condition]] == 0,], gene, a))
    scoresT1_list <- lapply(alphas, function(a) compute_cqr_scores(qrT1, cal_combined_df[cal_combined_df[[obs_condition]] == 1,], gene, a))
    print(scoresT0_list)
    # build intervals (vectorized version you just updated)
    test_idx0 <- which(test_g[[obs_condition]] == 0)
    test_idx1 <- which(test_g[[obs_condition]] == 1)
    out_mat <- matrix(NA_real_, nrow = length(alphas) * nrow(test_g), ncol = 2)

    rowpos <- 1
    for (ai in seq_along(alphas)) {
      a <- alphas[ai]
      intervals_full <- matrix(NA_real_, nrow = nrow(test_g), ncol = 2)
      if (length(test_idx0) > 0) {
        intervals_full[test_idx0, ] <- build_intervals(
          test_g[test_idx0, , drop = FALSE],
          qrT1, scoresT1_list[[ai]], wC1, w_test[test_idx0],
          a, gene, gene_names, 0)
      }
      if (length(test_idx1) > 0) {
        intervals_full[test_idx1, ] <- build_intervals(
          test_g[test_idx1, , drop = FALSE],
          qrT0, scoresT0_list[[ai]], wC0, w_test[test_idx1],
          a, gene, gene_names, 1
        )
      }
      out_mat[rowpos:(rowpos + nrow(test_g) - 1), ] <- intervals_full
      rowpos <- rowpos + nrow(test_g)
    }

    df_out <- as.data.frame(out_mat)
    colnames(df_out) <- c("lower","upper")
    df_out$alpha   <- rep(alphas, each = nrow(test_g))
    df_out$gene    <- gene
    df_out$cell_id <- rep(seq_len(nrow(test_g)),length(alphas))
    df_out$covered <- factor(ifelse(df_out$lower < 0 & 0 < df_out$upper, "inside", "outside"))
    df_out
  }, BPPARAM = param3)
  do.call(rbind, gene_results)
}

