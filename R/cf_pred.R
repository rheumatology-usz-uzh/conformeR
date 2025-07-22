#' Imputation of counterfactual scenario and computation of conformal p-values for each gene x patient id.
#' @import BiocParallel
#' @import dplyr
#' @import probably
#'
#' @param sce an original sc-RNA-seq dataset of class SingleCellExperiment with obs_condition.
#' @param obs_condition name of the column in colData(sce) for the cell's observed condition specified as a string.
#' @param prop_train size in proportion of original data set for training set.
#' @param prop_proper size in proportion of training data set for proper set.
#' @param spacing level of precision for alpha's grid.
#' @param cores number of workers for BiocParallel job.
#'
#' @return a tibble of size n_cells x n_genes containing conformal p-values.
#' @export


cf_pred <- function(sce, obs_condition, spacing = 0.1, prop_train = 0.5, prop_proper = 0.75, cores = 32) {
  suppressMessages(
    suppressWarnings({
      params <- BatchtoolsParam(workers = cores)
      set.seed(123)

      # Prepare data
      data <- as_tibble(t(assay(sce, "logcounts"))) |>
        mutate(Tr = factor(colData(sce)[[obs_condition]]))
      names(data) <- c(rownames(sce), "Tr")

      # Split data
      split1 <- initial_split(data, prop = prop_train, strata = Tr)
      train <- training(split1)
      test <- testing(split1)

      split2 <- initial_split(train, prop = prop_proper, strata = Tr)
      proper <- training(split2)
      cal <- testing(split2)

      # Pre-split subsets by treatment
      properT0 <- filter(proper, Tr == levels(data$Tr)[1])
      properT1 <- filter(proper, Tr == levels(data$Tr)[2])
      calT0 <- filter(cal, Tr == levels(data$Tr)[1])
      calT1 <- filter(cal, Tr == levels(data$Tr)[2])
      testT0 <- filter(test, Tr == levels(data$Tr)[1])
      testT1 <- filter(test, Tr == levels(data$Tr)[2])
      alphas <- seq(spacing, 1 - spacing, spacing)

      # Launch parallel processing
      res <- bptry({
        bplapply(seq_len(nrow(sce)), function(i) {
          model_spec <- rand_forest() |> set_engine("ranger") |> set_mode("regression")

          response_var <- names(data)[i]
          predictors <- setdiff(names(data), c(response_var, "Tr"))
          formula_obj <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))

          # Fit two models
          fitT0 <- workflow() |> add_formula(formula_obj) |> add_model(model_spec) |> fit(properT0)
          fitT1 <- workflow() |> add_formula(formula_obj) |> add_model(model_spec) |> fit(properT1)

          # Conformal intervals
          intT0 <- int_conformal_split(fitT0, calT0)
          intT1 <- int_conformal_split(fitT1, calT1)

          int_builder <- function(alpha) {
            t0 <- cbind(testT0, predict(intT0, testT0, level = alpha))
            t1 <- cbind(testT1, predict(intT1, testT1, level = alpha))
            out <- rbind(t0, t1)
            out$alpha <- alpha
            out$id <- 1:nrow(out)
            return(out)
          }

          intervals <- lapply(alphas, int_builder) |> bind_rows()
          gene_name <- names(intervals)[i]
          pvalues <- intervals |>
            group_by(id) |>
            summarise(
              pvalue = (1 + sum((.pred_lower < .data[[gene_name]]) & (.data[[gene_name]] < .pred_upper), na.rm=T))/(1 + n()),
              .groups = "drop"
            )

          return(pvalues)
        }, BPPARAM = params)
      })
      res <- res |>
        bind_cols() |>
        select(-matches("id"))
      colnames(res) <- rownames(sce)
      return(res)
    })
  )
}
