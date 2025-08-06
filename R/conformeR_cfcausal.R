#' Imputation of counterfactual scenario and computation of conformal p-values for each gene x patient id.
#' @import BiocParallel
#' @import dplyr
#' @importFrom cfcausal conformalIte
#'
#' @param sce SingleCellExperiment with replicate_id, obs_condition, cell_type in colData
#' @param replicate_id column name for biological replicate (string)
#' @param obs_condition column name for observed condition (string)
#' @param cell_type column name for cell types (string)
#' @param size_train proportion for training set.
#' @param cutoff confidence level for FDR
#' @param size_cal proportion for calibration set.
#' @param spacing level of precision for alpha's grid.
#' @param cores number of workers for BiocParallel job.
#'
#' @return a data frame of dimension number of `conf_group` x number of genes containing FDR(level=`cutoff`) for each gene in a given  `conf_group`
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

    # Iterate over conformal groups
    results <- lapply(groups, function(g) {
      train_group <- train_sce[, colData(train_sce)$conf_group == g]
      test_group  <- test_sce[, colData(test_sce)$conf_group == g]

      data_train <- as_tibble(t(assay(train_group, "logcounts"))) |>
        mutate(Tr = colData(train_group)[[obs_condition]])
      data_test <- as_tibble(t(assay(test_group, "logcounts"))) |>
        mutate(
          Tr = colData(test_group)[[obs_condition]],
          cell_id = test_group$cell_id
        )

      gene_names <- rownames(train_group)
      names(data_train) <- c(gene_names, "Tr")
      names(data_test)  <- c(gene_names, "Tr", "cell_id")

      # Per-gene processing
      gene_pvalues <- BiocParallel::bplapply(seq_along(gene_names), function(i) {
        gene <- gene_names[i]

        fun_cfcausal <- function(alpha) {
          X_train <- data_train |> select(-all_of(c(gene, "Tr"))) |> as.matrix()
          Y_train <- data_train[[gene]]
          T_train <- data_train$Tr

          X_test <- data_test |> select(-all_of(c(gene, "Tr", "cell_id"))) |> as.matrix()
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

          tibble::as_tibble(cbind(as.data.frame(colData(test_group)), int)) |>
            mutate(alpha = alpha)
        }

        intervals <- lapply(alphas, fun_cfcausal) |> bind_rows()

        pvals <- intervals |>
          group_by(cell_id) |>
          summarise(
            pvalue = (1 + sum((lower < 0) & (0 < upper), na.rm = TRUE)) / (1 + n()),
            .groups = "drop"
          ) |>
          mutate(gene = gene)

        return(pvals)
      }, BPPARAM = params)

      # Combine gene-wise pvalues into single tibble
      bind_rows(gene_pvalues) |>
        mutate(conf_group = g)
    })

    # Combine all groups
    final_results <- bind_rows(results)

    return(list(final_results, test_sce))
  }))
}
