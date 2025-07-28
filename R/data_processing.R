#' 1) creation of conformal groups `conf_group`.
#' 2) splitting of the original dataset using `initial_split` and balancing on `conf_group`.
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom dplyr group_indices mutate
#' @importFrom rlang syms
#' @importFrom rsample initial_split training testing
#'
#' @param sce SingleCellExperiment with replicate_id, obs_condition, cell_type in colData
#' @param replicate_id column name for biological replicate (string)
#' @param obs_condition column name for observed condition (string)
#' @param cell_type column name for cell types (string)
#' @param size_train proportion for training set
#' @param size_cal proportion for calibration set of training set.
#'
#' @return list with original sce, proper training, calibration, and test sets
#' @export


data_processing <- function(sce,
                            replicate_id,
                            obs_condition,
                            cell_type,
                            size_train = 0.5,
                            size_cal = 0.25) {
  stopifnot(is(sce, "SingleCellExperiment"))
  set.seed(1234)

  # Ensure factors
  colData(sce)[[replicate_id]] <- as.factor(colData(sce)[[replicate_id]])
  colData(sce)[[obs_condition]] <- as.factor(colData(sce)[[obs_condition]])
  colData(sce)[[cell_type]] <- as.factor(colData(sce)[[cell_type]])

  # Convert to data.frame and build conformal groups
  coldata_df <- as.data.frame(colData(sce)) %>%
    mutate(conf_group = factor(paste0(.data[[replicate_id]], " x ", .data[[cell_type]])))


  # Assign conformal groups back to colData
  colData(sce)$conf_group <- coldata_df$conf_group

  # Row index for traceability
  coldata_df$row <- seq_len(nrow(coldata_df))

  # Split into training + test, stratified
  split1 <- initial_split(coldata_df, prop = size_train, strata = conf_group)
  remaining <- training(split1)
  test_idx <- testing(split1)$row

  # Further split remaining into proper + calibration
  if (size_cal == 0) {
    train_idx <- training(split1)$row
    train_set <- sce[, train_idx]
    test_set <- sce[, test_idx]

    return(list(
    original = sce,
    train_set = train_set,
    test_set = test_set
  ))
  }

  else {
  split2 <- initial_split(remaining, prop = 1 - size_cal, strata = conf_group)
  proper_idx <- training(split2)$row
  cal_idx <- testing(split2)$row

  # Subset SCE
  train_set <- sce[, proper_idx]
  cal_set <- sce[, cal_idx]
  test_set <- sce[, test_idx]

  return(list(
    original = sce,
    train_set = train_set,
    cal_set = cal_set,
    test_set = test_set
  )) }
}

