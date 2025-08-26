#' Prediction interval builder for difference in conditions, based on counterfactual scenario prediction.
#' @import BiocParallel
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

  suppressMessages(suppressWarnings({
    set.seed(123)
    param2 <- MulticoreParam(workers = cores)
    param3 <- MulticoreParam(workers = cores)

    # Split data
    splits <- data_processing(sce, replicate_id, obs_condition, cell_type, size_train, size_cal)
    properT0 <- subset_by_condition(splits$train_set, obs_condition, 0)
    properT1 <- subset_by_condition(splits$train_set, obs_condition, 1)
    calT0    <- subset_by_condition(splits$cal_set, obs_condition, 0)
    calT1    <- subset_by_condition(splits$cal_set, obs_condition, 1)
    test     <- splits$test_set

    groups     <- levels(colData(splits$original)$conf_group)
    alphas     <- seq(spacing, 1 - spacing, spacing)
    gene_names <- rownames(sce)

    # Iterate over groups
    results <- lapply(groups, function(g) {
      # Subset by group
      gsets <- list(
        T0 = properT0[, properT0$conf_group == g],
        T1 = properT1[, properT1$conf_group == g],
        C0 = calT0[, calT0$conf_group == g],
        C1 = calT1[, calT1$conf_group == g],
        Te = test[, test$conf_group == g]
      )
      stopifnot(ncol(gsets$T0) > 2, ncol(gsets$T1) > 2)

      # Convert to tibbles
      dataT0 <- sce_to_tibble(gsets$T0, obs_condition, gene_names)
      dataT1 <- sce_to_tibble(gsets$T1, obs_condition, gene_names)
      dataC0 <- sce_to_tibble(gsets$C0, obs_condition, gene_names)
      dataC1 <- sce_to_tibble(gsets$C1, obs_condition, gene_names)
      dataTe <- sce_to_tibble(gsets$Te, obs_condition, gene_names)

      # Loop over genes
      gene_pvalues <- bplapply(seq_along(gene_names), function(gene_idx) {
        gene <- gene_names[gene_idx]

        # Propensity score
        ps_model <- prop_score(rbind(dataT0, dataT1), gene_idx, gene_names)
        ps_cal   <- predict(ps_model, rbind(dataC0, dataC1), type = "prob") |> pull(.pred_1)
        w_cal    <- (1 - ps_cal) / ps_cal
        wC0      <- w_cal[1:nrow(dataC0)]
        wC1      <- w_cal[(nrow(dataC0) + 1):length(w_cal)]

        # Quantile regressions
        qrT0 <- train_qr(dataT0, gene, gene_names)
        qrT1 <- train_qr(dataT1, gene, gene_names)

        # Iterate over alphas
        df_grid <- bplapply(alphas, function(alpha) {
          # Compute scores
          scoresT0 <- compute_cqr_scores(qrT0, dataC0, gene_idx, alpha)
          scoresT1 <- compute_cqr_scores(qrT1, dataC1, gene_idx, alpha)

          intervals <- lapply(seq_len(nrow(dataTe)), function(i) {
            prop_x <- predict(ps_model, dataTe[i, ], type = "prob") |> pull(.pred_1)
            w_x    <- (1 - prop_x) / prop_x
            if (dataTe$Tr[i] == 0) {
              build_intervals(dataTe[i, ], qrT1, scoresT1, wC1, w_x, alpha, gene_idx, gene_names, 0)
            } else {
              build_intervals(dataTe[i, ], qrT0, scoresT0, wC0, w_x, alpha, gene_idx, gene_names, 1)
            }
          })
          do.call(rbind, intervals) |> as.data.frame() |> setNames(c("lower", "upper"))
        }, BPPARAM = param2)

        do.call(rbind, df_grid) |>
          as.data.frame() |>
          mutate(cell_id = rep(seq_len(nrow(dataTe)), length(alphas)),
                 alpha   = rep(alphas, each = nrow(dataTe)),
                 gene    = gene,
                 covered = factor(ifelse(lower < 0 & 0 < upper, "inside", "outside")))
      }, BPPARAM = param3)

      do.call(rbind, gene_pvalues) |> as.data.frame() |> mutate(conf_group = g)
    })

    do.call(rbind, results) |> as.data.frame()
  }))
}
