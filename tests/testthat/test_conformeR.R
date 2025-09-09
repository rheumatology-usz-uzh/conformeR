library(testthat)
library(SingleCellExperiment)

# Create a small mock dataset
mock_sce <- SingleCellExperiment(
  assays = list(logcounts = matrix(rnorm(200), nrow = 5, ncol = 40)),
  colData = DataFrame(
    Tr = rep(c(0,1),20),
    replicate = rep(c("A","B"),each=20),
    cell_type = rep("T-cells",40)
))

rownames(mock_sce) <- paste0("gene", 1:5)

test_that("eta_conformeR returns numeric scalar", {
  scores <- runif(4)
  weights_cal <- runif(4)
  weights_test <- 0.5
  expect_type(eta_conformeR(scores, weights_cal, weights_test,alpha=.05), "double")
})

test_that("interval_to_pval returns correct columns", {
  intervals_df <- data.frame(
    cell_id = 1:4,
    conf_group = factor(c("A","A","B","B")),
    gene = rep("gene1",4),
    lower = runif(4),
    upper = runif(4),
    covered = sample(c("inside","outside"), 4, replace = TRUE)
  )
  res <- interval_to_pval(intervals_df)
  expect_true(all(c("conf_group","gene","pvalue") %in% colnames(res)))
})

test_that("fdr returns list with expected components", {
  intervals_df <- data.frame(
    cell_id = 1:4,
    conf_group = factor(c("A","A","B","B")),
    gene = rep("gene1",4),
    lower = runif(4),
    upper = runif(4),
    covered = sample(c("inside","outside"), 4, replace = TRUE)
  )
  res <- fdr(intervals_df, cores = 1)
  expect_true(all(c("gene_fdr_fun","comb_fdr","pvaldf") %in% names(res)))
})

test_that("conformeR returns list with expected components", {
  res <- conformeR(
    sce = mock_sce,
    obs_condition = "Tr",
    replicate_id = "replicate",
    cell_type = "cell_type",
    cf_method = "conformeR",
    cores = 32
  )
  expect_true(all(c("pred_intervals","conformeR_pvalues","gene_fdr_fun","comb_fdr") %in% names(res)))
})

