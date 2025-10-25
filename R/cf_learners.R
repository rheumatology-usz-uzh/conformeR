#' Fit a propensity score model
#'
#' Internal helper function that trains a propensity score model using all
#' observed logcounts except the gene of interest. The model is estimated
#' with a random forest classifier from the \pkg{ranger} engine.
#'
#' @import workflows
#' @import parsnip
#'
#' @param proper_set A data frame containing predictors (logcounts) and the
#'   treatment indicator \code{Tr}.
#' @param gene Character string giving the name of the target gene.
#' @param gene_names Character vector of all gene names.
#'
#' @return A fitted workflow object containing the trained propensity score model.
prop_score <- function(proper_set, gene, gene_names, obs_condition){
  model_prop_score <-
    rand_forest() |>
    set_engine("ranger") |>
    set_mode("classification")


weight_function<-function(fitted_model, data, gene_names){
  cal <- rbind(
    data$C0,
    data$C1
  )
  test<-data$Te
  
  #for justine to edit
  columns_to_remove=list("replicate_id", "cell_type") #only Genes and obs_condition columns left
  proper_all <- proper_all[, !names(proper_all) %in% c(columns_to_remove)]
  
  weight_cal_dict <- list()
  weight_test_dict <- list()
  
  for (gene_of_interest in gene_names){
    #cal set
    prop_score_cal<-predict_without_gene(model=fitted_model,gene_name=gene_of_interest,data=cal)
    #weights
    weights_cal <- (1 - prop_score_cal) / prop_score_cal
    weight_cal_dict[[gene_of_interest]] <- weights_cal
    
    # test set
    prop_score_test<-predict_without_gene(model=fitted_model, gene_name=gene_of_interest,data=test)
    #weights
    weights_test <- (1 - prop_score_test) / prop_score_test
    weight_test_dict[[gene_of_interest]] <- weights_test
    
  }
  
  return(list(
    weight_cal = weight_cal_dict,
    weight_test = weight_test_dict
  ))
}

  # Build formula dynamically with the actual treatment column
  response <- paste0("as.factor(", obs_condition, ")")
  predictors <- paste(gene_names[gene_names != gene], collapse = " + ")
  f <- as.formula(paste(response, "~", predictors))

  prop_scores <- workflow() |>
    add_formula(f) |>
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
#' @importFrom ranger ranger
#'
#' @param data A data frame containing predictors (logcounts) and outcomes.
#' @param gene Character string giving the name of the target gene.
#' @param gene_names Character vector of all gene names.
#'
#' @return A fitted \code{ranger} object trained in quantile regression mode.
train_qr <- function(data, gene, gene_names) {
  print(data)
  predictors <- setdiff(gene_names, gene)
  formula <- as.formula(paste(make.names(gene), "~", paste(make.names(predictors), collapse = "+")))
  ranger(formula, data = data, quantreg = TRUE)
}


