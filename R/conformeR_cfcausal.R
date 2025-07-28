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

conformeR_cfcausal <- function(sce,
                    obs_condition,
                    replicate_id,
                    cell_type,
                    spacing = 0.1,
                    size_train = 0.5,
                    cutoff = .05,
                    cores = 32) {
  suppressMessages(suppressWarnings({
    set.seed(123)
    params <- BatchtoolsParam(workers = cores)

    # Split data using data_processing()
    splits <- data_processing(sce, replicate_id, obs_condition, cell_type, size_train, size_cal=0)
    train_sce <- splits$train_set
    test_sce <- splits$test_set

    # Get list of conformal groups
    groups <- levels(colData(train_sce)$conf_group)
    alphas <- seq(spacing, 1 - spacing, spacing)

    # Process each group
    results <- lapply(groups, function(g) {
      # Subset by group
      train_group <- train_sce[, colData(train_sce)$conf_group == g]
      test_group  <- test_sce[, colData(test_sce)$conf_group == g]

      # Extract logcounts and prepare data
      data_train <- as_tibble(t(assay(train_group, "logcounts"))) |>
       mutate(Tr = colData(train_group)[[obs_condition]])
      data_test <- as_tibble(t(assay(test_group, "logcounts"))) |>
        mutate(Tr = colData(test_group)[[obs_condition]])

      names(data_train) <- c(rownames(train_group), "Tr")
      names(data_test) <- c(rownames(test_group), "Tr")

      # Loop over genes
      gene_pvalues <- bplapply(seq_len(nrow(train_group)), function(i) {
        # Standard function performing for gene i, alpha the cfcausal WF
        fun_cfcausal <- function(alpha) {
          # Get gene name corresponding to index i
          gene_name <- rownames(train_group)[i]

          # Prepare X, Y, T for training
          X_train <- data_train %>%
            dplyr::select(-all_of(c(gene_name, "Tr"))) %>%
            as.matrix()

          Y_train <- data_train[[gene_name]]
          T_train <- data_train$Tr

          # Prepare X, Y, T for testing
          X_test <- data_test %>%
            dplyr::select(-all_of(c(gene_name, "Tr"))) %>%
            as.matrix()

          Y_test <- data_test[[gene_name]]
          T_test <- data_test$Tr

          # Fit conformal model
          mod <- conformalIte(
            X = X_train, Y = Y_train, T = T_train,
            alpha = alpha, algo = "counterfactual",
            type = "CQR", side = "two", quantiles = c(alpha / 2, 1 - alpha / 2),
            outfun = "quantRF", useCV = FALSE
          )

          # Predict on test set
          int <- mod(X = X_test, Y = Y_test, T = T_test)

          test <- cbind(as.data.frame(colData(test_group)), int) %>%
            dplyr::mutate(alpha = alpha, id = seq_len(nrow(.)))
          return(test)
        }
        intervals <- lapply(alphas, fun_cfcausal) |> bind_rows()
        print(intervals)
        pvalues <- intervals |>
          group_by(id) |>
          summarise(
            pvalue = (1 + sum((lower < 0) & (0 < upper), na.rm = TRUE)) / (1 + n())
          )
        return(pvalues$pvalue)
      }, BPPARAM = params)
      # Collapse into a row (1 row per group, length = number of genes)
      res <- matrix(unlist(gene_pvalues),ncol=length(gene_pvalues),byrow=T)
      res <- as.data.frame(res)
      colnames(res) <- rownames(sce)
      res$conf_group <- rep(g,nrow(res))
      return(res)
    })

    # Combine result across conf_group
    res <- do.call(rbind, results)
    res$conf_group <- as.factor(res$conf_group)
    fdr_res <- fdr(res,cutoff,cores)
    return(fdr_res)
  }))
}
