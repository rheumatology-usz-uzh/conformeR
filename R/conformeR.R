# Outputs a matrix of (adjusted) p-values encoding the difference between the observed and predicted expression levels for each pair gene x cell in a test set given training set and calibration set, all three randomly built from an original sc-RNA-seq dataset of class SingleCellExperiment with replicate_id, obs_condition, cell_type annotated in colData.
#
# @param sce an original sc-RNA-seq dataset of class SingleCellExperiment with replicate_id, obs_condition, cell_type annotated in colData.
# @param replicate_id name of the column in colData(sce) for biological replicate id specified as a string.
# @param obs_condition name of the column in colData(sce) for the cell's observed condition specified as a string.
# @param cell_type name of the column in colData(sce) for cell-types specified as a string.
# @param size_train size in proportion of original data set for training set.
# @param size_cal size in proportion of original data set for calibration set.
# @param alpha level of confidence between 0 and 1.
# @param cores number of workers for BiocParallel job.
#
# @return a test_set with two new assays named "pvalues" resp. "adj_pvalues" containing a matrix of p_values for each pair gene x cell resp. BH-adjusted p_values for each pair gene x cell..
#

conformeR <- function(sce,
                      replicate_id,
                      obs_condition,
                      cell_type,
                      size_train = .25,
                      size_cal = .25,
                      alpha = .1,
                      cores = 32,
                      file_name) {
  param <- BatchtoolsParam(workers = cores)
  data_processed <- data_processing(
    sce,
    replicate_id,
    obs_condition,
    cell_type,
    size_train = .25,
    size_cal = .25
  )
  train_set <- data_processed$train_set
  cal_set <- data_processed$cal_set
  test_set <- data_processed$test_set
  lemur_imputed <- lemur_imputation(data_processed, replicate_id, obs_condition)
  train_set_imp <- lemur_imputed$pred_trainset
  cal_set_imp <- lemur_
  imputed$pred_calset
  pval <- bptry({
    bplapply(seq_len(ncol(test_set)), function(j) {
      sapply(seq_len(nrow(test_set)), function(i) {
        conf_pval(
          assay(train_set, "logcounts")[i, colData(train_set)$conf_group == colData(test_set)$conf_group[j]],
          train_set_imp[i, colData(train_set)$conf_group ==
                          colData(test_set)$conf_group[j]],
          assay(cal_set, "logcounts")[i, colData(cal_set)$conf_group ==
                                        colData(test_set)$conf_group[j]],
          cal_set_imp[i, colData(cal_set)$conf_group ==
                        colData(test_set)$conf_group[j]],
          assay(test_set, "logcounts")[i, j],
          alpha
        )
      })
    }, BPPARAM = param)
  })
  clean_pval <- lapply(pval, function(sublist) {
    lapply(sublist, function(x) {
      if (is.numeric(x)) {
        x
      } else {
        NA
      }
    })
  })
  pval_mat <- matrix(
    data = unlist(clean_pval),
    ncol = ncol(test_set),
    nrow = nrow(test_set)
  )
  adjpval_mat <- matrix(ncol = ncol(test_set), nrow = nrow(test_set))
  for (i in seq_len(nrow(test_set))) {
    df <- data.frame(
      pval = pval_mat[i, ],
      group = colData(test_set)$conf_group,
      original_order = seq_along(pval_mat[i, ])
    )
    df <- df %>%
      group_by(group) %>%
      mutate(adj_pval = p.adjust(pval, method = "BH")) %>%
      ungroup() %>%
      arrange(original_order)
    adjpval_mat[i, ] <- df$adj_pval
  }
  assay(test_set, "pvalues") <- pval_mat
  assay(test_set, "adj_pvalues") <- adjpval_mat
  return(test_set)
}
