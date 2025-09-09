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
    spacing = 0.01, size_train = 0.5, size_cal = 0.25,
    cores = 32, gene_chunk = 200) {

  suppressMessages(suppressWarnings({
    set.seed(123)
    param3 <- MulticoreParam(workers = cores)

    splits <- data_processing(sce, replicate_id, obs_condition, cell_type,
                              size_train, size_cal)
    dt_train <- splits$train_set |> as.data.frame()
    dt_cal   <- splits$cal_set |> as.data.frame()
    dt_test  <- splits$test_set |> as.data.frame()

    groups     <- dt_train |> select(conf_group) |> unique()
    alphas     <- seq(spacing, 1 - spacing, spacing)
    gene_names <- sce |> rownames()

    res_list <- lapply(groups, function(g) process_group(g, dt_train, dt_cal, dt_test, gene_names, alphas, obs_condition, gene_chunk, param3))
    res_df <- do.call(rbind, res_list)
    rownames(res_df) <- NULL
    res_df
    }))
}
