#' Two-step data processing of the original dataset ;
#' 1) creation of conformal groups.
#' 2) splitting of the original dataset
#'
#' @param sce an original sc-RNA-seq dataset of class SingleCellExperiment with replicate_id, obs_condition, cell_type annotated in colData.
#' @param replicate_id name of the column in colData(sce) for biological replicate id specified as a string.
#' @param obs_condition name of the column in colData(sce) for the cell's observed condition specified as a string.
#' @param cell_type name of the column in colData(sce) for cell-types specified as a string.
#' @param size_train size in proportion of original data set for training set.
#' @param size_cal size in proportion of original data set for calibration set.
#'
#' @return a list of SingleCellExperiment objects with original dataset, training set, calibration set, test set, including conformal groups id as a new column in colData()
#' @export

data_processing <- function(sce,
                            replicate_id,
                            obs_condition,
                            cell_type,
                            size_train = .25,
                            size_cal = .25) {
  set.seed(1234)
  stopifnot(is(sce, "SingleCellExperiment"))
  colData(sce)[[replicate_id]] <- as.factor(colData(sce)[[replicate_id]])
  colData(sce)[[obs_condition]] <- as.factor(colData(sce)[[obs_condition]])
  colData(sce)[[cell_type]] <- as.factor(colData(sce)[[cell_type]])
  # Build the conformal groups.
  conf_group <- (as.data.frame(colData(sce)) %>%
                   mutate(conf_group = group_indices(., !!!syms(
                     c(replicate_id, cell_type)
                   ))))$conf_group
  colData(sce)$conf_group <- as.factor(conf_group)
  # Split data.
  id_train <- sample(1:ncol(sce), floor(size_train * ncol(sce)),
                     replace =
                       FALSE
  )
  id_cal <- sample(subset(1:ncol(sce), !1:ncol(sce) %in% id_train),
                   floor(size_cal * ncol(sce)),
                   replace = FALSE
  )
  id_test <- subset(1:ncol(sce), (!1:ncol(sce) %in% id_train) &
                      (!1:ncol(sce) %in% id_cal))
  # Prepare each training / calibration / test set.
  train_set <- sce[, id_train]
  cal_set <- sce[, id_cal]
  test_set <- sce[, id_test]
  return(list(
    original = sce,
    train_set = train_set,
    cal_set = cal_set,
    test_set = test_set
  ))
}
