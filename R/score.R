# Score function underlying the conformal prediction step.
#
# @param obs_train a vector of observed logcounts for the same gene, and on which the quantile regression forest is trained.
# @param pred_train a vector of predicted lemur logcounts for the same gene and of the same length as obs_train.
# @param obs_cal a vector of observed logcounts for the same gene, to be passed as argument of the score function.
# @param pred_cal a vector of predicted lemur logcounts for the same gene and of the same length as obs_cal.
# @param alpha level of confidence between 0 and 1.
#
# @return a vector of score values of the same length as obs_cal.
#

score <- function(obs_train,
                  pred_train,
                  obs_cal,
                  pred_cal,
                  alpha) {
  stopifnot((length(obs_train) == length(pred_train)) &
              (length(obs_cal) == length(pred_cal)))
  qrmodel <- quantregForest(
    x = data.frame(obs_counts = obs_train), y =
      pred_train
  )
  quant <- predict(qrmodel,
                   data.frame(obs_counts = obs_cal),
                   what = c(alpha / 2, 1 - alpha / 2)
  )
  res <- apply(cbind(quant[, 1] - pred_cal, pred_cal - quant[, 2]), 1, max)
  return(res)
}
