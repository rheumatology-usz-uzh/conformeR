#' Prediction interval builder for difference in conditions, based on counterfactual scenario prediction.
#' @importFrom BiocParallel MulticoreParam
#' @import dplyr
#'
#' @param sce SingleCellExperiment with replicate_id, obs_condition, cell_type in colData
#' @param replicate_id column name for biological replicate (string)
#' @param obs_condition column name for observed condition (string)
#' @param cell_type column name for cell types (string)
#' @param size_train proportion for training set.
#' @param size_cal proportion of training set for calibration set.
#' @param spacing level of precision for alpha's grid.
#' @param cores number of workers for BiocParallel job.
#'
#' @return a data frame of prediction intervals for difference in observed vs. counterfactual logcount for each observed gene x cell x confidence level
#' @export

conformeR_cf <- function(
    sce, obs_condition, replicate_id, cell_type,
    spacing = 0.01, size_train = 0.5, size_cal = 0.25, cores = 32) {

  set.seed(123)
  param <- MulticoreParam(workers = cores, RNGseed = 123)

  # Split data
  splits <- data_processing(sce, replicate_id, obs_condition, cell_type, size_train, size_cal)
  properT0 <- subset_by_condition(splits$train_set, obs_condition, 0)
  properT1 <- subset_by_condition(splits$train_set, obs_condition, 1)
  calT0    <- subset_by_condition(splits$cal_set, obs_condition, 0)
  calT1    <- subset_by_condition(splits$cal_set, obs_condition, 1)
  test     <- splits$test_set

  groups     <- levels(splits$train_set$conf_group)
  alphas     <- seq(spacing, 1 - spacing, spacing)
  gene_names <- rownames(sce)

  # Iterate over groups
  results <- lapply(groups, function(g) {
    # Subset by group
    gsets <- list(
      T0 = properT0[properT0$conf_group == g, ],
      T1 = properT1[properT1$conf_group == g, ],
      C0 = calT0[calT0$conf_group == g, ],
      C1 = calT1[calT1$conf_group == g, ],
      Te = test[test$conf_group == g,]
    )
    stopifnot(ncol(gsets$T0) > 2, ncol(gsets$T1) > 2)

    idx0 <- which(gsets$Te[[obs_condition]] == 0)
    idx1 <- which(gsets$Te[[obs_condition]] == 1)

    # Loop over genes in parallel
    gene_pvalues <- bplapply(gene_names, function(gene) {

      # Train quantile regressions
      qrT0 <- train_qr(gsets$T0, gene, gene_names)
      qrT1 <- train_qr(gsets$T1, gene, gene_names)

      # Propensity score
      ps_model <- prop_score(rbind(gsets$T0, gsets$T1), gene, gene_names, obs_condition)
      ps_cal   <- predict(ps_model, rbind(gsets$C0, gsets$C1), type = "prob") |> pull(.pred_1)
      ps_test <- predict(ps_model, gsets$Te, type = "prob") |> pull(.pred_1)
      w_cal    <- (1 - ps_cal) / ps_cal
      w_test <- (1 - ps_test) / ps_test
      wC0      <- w_cal[1:nrow(gsets$C0)]
      wC1      <- w_cal[(nrow(gsets$C0) + 1):length(w_cal)]

      # calibration scores
      scoresT0 <- compute_cqr_scores(qrT0, gsets$C0, gene, alphas)
      scoresT1 <- compute_cqr_scores(qrT1, gsets$C1, gene, alphas)

      int0 <- build_intervals(
        gsets$Te[idx0, , drop = FALSE],
        qrT1, scoresT1, wC1, w_test[idx0],
        alphas, gene, gene_names, 0)

      int1 <- build_intervals(
        gsets$Te[idx1, , drop = FALSE],
        qrT0, scoresT0, wC0, w_test[idx1],
        alphas, gene, gene_names, 0)

      int <- rbind(int0,int1) |>
        mutate(gene=gene)

      int
    }, BPPARAM = param)
    # Combine all genes for the group
    data.table::rbindlist(gene_pvalues)[, conf_group := g]
  })

  # Combine all groups
  data.table::rbindlist(results) |>
    mutate(covered = sign(lower)==sign(upper))
}
