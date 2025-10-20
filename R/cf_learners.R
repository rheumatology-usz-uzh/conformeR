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

 # light propensity score estimator
  #colnames is the list of the column of sce (obs is not included)
  colnames<-list("Batch", "Group","sim.Lib.Size","cell","sizeFactor","replicate_id","cell_type","conf_group") #<-example

propensity_score<-function(model, sce, colnames){
  
  useless_col<-colnames
  splits<-sce
  
  #LITTLE BIT OF PREPROCESSING
  splits[[1]]$obs_condition <- factor(ifelse(splits[[1]]$obs_condition == 2, 1, 0), levels = c(0, 1))
  splits[[2]]$obs_condition <- factor(ifelse(splits[[2]]$obs_condition == 2, 1, 0), levels = c(0, 1))
  splits[[3]]$obs_condition <- factor(ifelse(splits[[3]]$obs_condition == 2, 1, 0), levels = c(0, 1))
  
  
  train <- as.data.frame(splits[[1]])
  train <- train[, !names(train) %in% useless_col]
  cal <- as.data.frame(splits[[2]])
  cal <- cal[, !names(cal) %in% c(useless_col,"obs_condition")]
  test <- as.data.frame(splits[[3]])
  test <- test[, !names(test) %in% c(useless_col,"obs_condition")]
  
  weight_cal_dict <- list()
  weight_test_dict <- list()

  all_genes <- colnames(cal)

    #training model only once
    model_fit <- training_model(model=model,proper_set=train, all_genes=all_genes,             obs_condition="obs_condition")
    
    for (gene_of_interest in all_genes){
      #cal set
      prop_score_cal<-predict_without_gene(model=model_fit,gene_name=gene_of_interest,        data=cal)
      #weights
      weights_cal <- (1 - prop_score_cal) / prop_score_cal
      weight_cal_dict[[gene_of_interest]] <- weights_cal
      
      # test set
      prop_score_test<-predict_without_gene(model=model_fit,gene_name=gene_of_interest,       data=test)
      #weights
      weights_test <- (1 - prop_score_test) / prop_score_test
      weight_test_dict[[gene_of_interest]] <- weights_test

    }
    
    return(list(
      weight_cal = weight_cal_dict,
      weight_test = weight_test_dict
    ))
}

res<-propensity_score(model="ridge", sce=sce1, colnames=colnames)
  res$weight_cal 
  res$weight_test
  

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


