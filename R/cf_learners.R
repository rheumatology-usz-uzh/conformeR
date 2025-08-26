#' Fit a propensity score model
#'
#' Internal helper function that trains a propensity score model using all
#' observed logcounts except the gene of interest. The model is estimated
#' with a random forest classifier from the \pkg{ranger} engine.
#'
#' @param proper_set A data frame containing predictors (logcounts) and the
#'   treatment indicator \code{Tr}.
#' @param gene_index Integer index of the gene of interest to be excluded
#'   from the predictors.
#' @param gene_names Character vector of all gene names.
#'
#' @return A fitted workflow object containing the trained propensity score model.
prop_score <- function(proper_set, gene_index, gene_names){
  model_prop_score <-
    rand_forest() |>
    set_engine("ranger") |>
    set_mode("classification")

  prop_scores <- workflow() |>
    add_formula(
      as.formula(
        paste("as.factor(Tr) ~", paste(gene_names[-gene_index], collapse = "+"))
      )) |>
    add_model(model_prop_score) |>
    fit(proper_set)

  return(prop_scores)
}

#' Train a quantile regression model
#'
#' Internal helper function that fits a quantile regression forest for a given
#' gene using the expression of all other genes as predictors. The model is
#' estimated using the \pkg{ranger} implementation of quantile regression.
#'
#' @param data A data frame containing predictors (logcounts) and outcomes.
#' @param gene Character string giving the name of the target gene.
#' @param gene_names Character vector of all gene names.
#'
#' @return A fitted \code{ranger} object trained in quantile regression mode.
train_qr <- function(data, gene, gene_names) {
  ranger(
    as.formula(paste(gene, "~", paste(gene_names[gene_names != gene], collapse = "+"))),
    data = data,
    quantreg = TRUE
  )
}
