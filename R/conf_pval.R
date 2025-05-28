#' Outputs a p-value encoding the difference between the observed and predicted expression levels for a new observed count.
#' All data given as input has to refer to a same gene and belonging to cells within a same conformal group.
#'
#' @importFrom quantregForest quantregForest
#'
#' @param obs_train a vector of observed logcounts, on which the quantile regression forest is trained.
#' @param pred_train a vector of predicted lemur logcounts of the same length as obs_train.
#' @param obs_cal a vector of observed logcounts, to be passed as argument of the score function.
#' @param pred_cal a vector of predicted lemur logcounts of the same length as obs_cal.
#' @param obs_test a new observed logcount.
#' @param alpha level of confidence between 0 and 1.
#'
#' @return a single p-value encoding the difference between the observed and predicted expression levels for a new observed count.
#' @export

conf_pval <- function(obs_train,
                      pred_train,
                      obs_cal,
                      pred_cal,
                      obs_test,
                      alpha = 0.1) {
  n <- length(obs_cal)
  qrmodel <- quantregForest(
    x = data.frame(obs_counts = obs_train), y =
      pred_train
  )
  quantpred <- predict(qrmodel,
                       data.frame(obs_counts = obs_test),
                       what = c(alpha / 2, 1 - alpha / 2)
  )
  qscore <- floor((n + 1) * (1 - alpha)) / n
  if (qscore <= 0 | qscore >= 1) {
    return(1)
  }
  t_qscore <- quantile(score(obs_train, pred_train, obs_cal, pred_cal, alpha),
                       probs = qscore
  )
  # Build the bounds so that the observed condition lies within the interval.
  bounds <- c(
    min(obs_test, quantpred[, 1] - t_qscore),
    max(obs_test, quantpred[, 2] + t_qscore)
  )
  while (obs_test %in% bounds == FALSE) {
    alpha <- alpha + 0.01
    quantpred <- predict(qrmodel,
                         data.frame(obs_counts = obs_test),
                         what = c(alpha / 2, 1 - alpha / 2)
    )
    # estimate quantiles for the score function.
    qscore <- floor((n + 1) * (1 - alpha)) / n
    if (qscore <= 0 | qscore >= 1) {
      return(1)
    }
    t_qscore <- quantile(score(obs_train, pred_train, obs_cal, pred_cal, alpha),
                         probs = qscore
    )
    # Build the bounds so that the observed condition lies within the interval.
    bounds <- c(
      min(obs_test, quantpred[, 1] - t_qscore),
      max(obs_test, quantpred[, 2] + t_qscore)
    )
  }
  return(alpha)
}
