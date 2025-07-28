#' Lemur imputation of training and calibration set using lemur fit trained on training set.
#' @importFrom lemur lemur
#' @importFrom lemur align_by_grouping
#'
#' @param data_processed an output of the function data_processing. If provided, it replaces train_set and cal_set.
#' @param replicate_id name of the column in colData() for biological replicate id specified as a string.
#' @param obs_condition name of the column in colData() for the cell's observed condition specified as a string.
#' @param cell_type name of the column in colData(sce) for cell-types specified as a string.
#'
#' @return a list of length 2 containing the lemur-imputed training matrix and the lemur-imputed calibration matrix (no observed conditions but imputed values alone).
#' @export

lemur_imputation <- function(data_processed,
                             replicate_id,
                             obs_condition,
                             cell_type) {
  train_set <- data_processed$train_set
  cal_set <- data_processed$cal_set

  fitlemur_train <- lemur(
    train_set,
    design = reformulate(c(replicate_id, obs_condition)),
    n_embedding = 15,
    verbose = FALSE )

  fitlemur_train <- align_by_grouping(
    fitlemur_train,
    grouping = colData(train_set)[[cell_type]],
    verbose = FALSE)

    fitlemur_cal <- lemur(
      cal_set,
      design = reformulate(c(replicate_id, obs_condition)),
      n_embedding = 15,
      verbose = FALSE )

    fitlemur_cal <- align_by_grouping(
      fitlemur_cal,
      grouping = colData(cal_set)[[cell_type]],
      verbose = FALSE)

    dm_cal_set <- cbind(
      model.matrix(reformulate(replicate_id), data = as.data.frame(colData(cal_set))),
      abs(as.numeric(colData(cal_set)[[obs_condition]]) - 2)
    )
    dm_train_set <- cbind(
      model.matrix(reformulate(replicate_id), data = as.data.frame(colData(train_set))),
      abs(as.numeric(colData(train_set)[[obs_condition]]) - 2)
    )
    pred_trainset <- predict(fitlemur_train, newdesign = dm_train_set)
    pred_calset <- predict(fitlemur_cal, newdesign = dm_cal_set)
    return(
      list(pred_trainset = pred_trainset, pred_calset = pred_calset)
    )
}
