#' Prediction interval builder for difference in conditions, based on counterfactual scenario prediction.
#' @import BiocParallel
#' @import dplyr
#' @importFrom cfcausal conformalIte
#'
#' @param sce SingleCellExperiment with replicate_id, obs_condition, cell_type in colData
#' @param replicate_id column name for biological replicate (string)
#' @param obs_condition column name for observed condition (string)
#' @param cell_type column name for cell types (string)
#' @param size_train proportion for training set.
#' @param spacing level of precision for alpha's grid.
#' @param cores number of workers for BiocParallel job.
#'
#' @return a data frame of prediction intervals for difference in observed vs. counterfactual logcount for each observed gene x cell x confidence level
#' @export

conformeR_cfcausal<- function(sce,
                              obs_condition,
                              replicate_id,
                              cell_type,
                              spacing = 0.01,
                              size_train = 0.5,
                              cores = 32) {
  suppressMessages(suppressWarnings({
    set.seed(123)
    params <- BiocParallel::BatchtoolsParam(workers = cores)

    # Split data
    splits <- data_processing(sce, replicate_id, obs_condition, cell_type, size_train, size_cal = 0)
    train_sce <- splits$train_set
    test_sce  <- splits$test_set
    test_sce$cell_id <- seq_len(ncol(test_sce))

    groups <- levels(colData(train_sce)$conf_group)
    alphas <- seq(spacing, 1 - spacing, spacing)
    gene_names <- rownames(sce)

    # Iterate over conformal groups
    results <- lapply(groups, function(g) {
      train_group <- train_sce[, colData(train_sce)$conf_group == g]
      test_group  <- test_sce[, colData(test_sce)$conf_group == g]

      data_train <- sce_to_tibble(train_group, obs_condition, gene_names)
      data_test <- sce_to_tibble(test_group, obs_condition, gene_names) |>
        mutate(
          cell_id = test_group$cell_id
        )

      names(data_train) <- c(gene_names, "Tr")
      names(data_test)  <- c(gene_names, "Tr", "cell_id")

      # Per-gene processing
      gene_intervals <- BiocParallel::bplapply(seq_along(gene_names), function(i) {
        gene <- gene_names[i]

        fun_cfcausal <- function(alpha) {
          X_train <- data_train |> dplyr::select(-dplyr::all_of(c(gene, "Tr"))) |> as.matrix()
          Y_train <- data_train[[gene]]
          T_train <- data_train$Tr

          X_test <- data_test |> dplyr::select(-dplyr::all_of(c(gene, "Tr", "cell_id"))) |> as.matrix()
          Y_test <- data_test[[gene]]
          T_test <- data_test$Tr

          mod <- conformalIte(
            X = X_train, Y = Y_train, T = T_train,
            alpha = alpha, algo = "counterfactual",
            type = "CQR", side = "two",
            quantiles = c(alpha / 2, 1 - alpha / 2),
            outfun = "quantRF", useCV = FALSE
          )

          int <- mod(X = X_test, Y = Y_test, T = T_test)

          as_tibble(cbind("cell_id"=data_test$cell_id,int)) |>
            mutate(alpha = alpha)
        }

        intervals <- lapply(alphas, fun_cfcausal) |>
          bind_rows() |>
          mutate(gene = gene,
                 covered = factor(ifelse(lower < 0 & 0 < upper, "inside", "outside")))
        return(intervals)
      }, BPPARAM = params)

      # Combine gene-wise pvalues into single tibble
      bind_rows(gene_intervals) |>
        mutate(conf_group = g)
    })
    bind_rows(results)
  }))
}
