#' Conformal inference with counterfactual methods for multi-condition single-cell data DE detection
#'
#' This function computes conformal prediction intervals and FDR-adjusted
#' results for gene expression data stored in a
#' \linkS4class{SingleCellExperiment} object. Results are aggregated at `cell_type` level.
#' Two counterfactual-based methods
#' are available: `"conformeR"` (default) and `"cfcausal"`.
#'
#' @param sce A \linkS4class{SingleCellExperiment} object containing
#'   gene expression data ("logcounts" assay).
#' @param obs_condition Character scalar. The name of the column in
#'   `colData(sce)` indicating the observed condition (e.g., treatment vs control).
#' @param replicate_id Character scalar. The name of the column in
#'   `colData(sce)` identifying biological replicates (e.g. patient ID).
#' @param cell_type Character scalar. The name of the column in
#'   `colData(sce)` specifying cell types.
#' @param cf_method Character scalar, either `"conformeR"` (default) or
#'   `"cfcausal"`. Determines which counterfactual method is used for interval estimation.
#' @param spacing Numeric scalar. Grid spacing for prediction intervals computation.
#'   Default is `0.01`.
#' @param size_train Numeric scalar between 0 and 1. Proportion of data used
#'   for training. Default is `0.5`.
#' @param size_cal Numeric scalar between 0 and 1. Proportion of data used
#'   for calibration. Default is `0.25`.
#' @param cores Integer. Number of parallel workers for computation.
#'   Default is `32`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{pred_intervals}{A data.frame of conformal prediction intervals
#'     for each gene x cell x level of confidence.}
#'   \item{conformeR_pvalues}{A data.frame containing per-gene x cell p-values.}
#'   \item{gene_fdr_fun}{A list of FDR estimation functions,
#'     one per gene × conformal group.}
#'   \item{comb_fdr}{A tibble summarizing combined FDR per gene.}
#' }
#'
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' sce <- mock_sce_object()  # example input
#' res <- conformeR(sce, obs_condition="treatment",
#'                  replicate_id="patient_id", cell_type="celltype",
#'                  cf_method="conformeR", cores=4)
#' }
#'
#' @export

conformeR <- function(sce,
                       obs_condition,
                       replicate_id,
                       cell_type, cf_method=c("conformeR","cfcausal"),
                       spacing = 0.01,
                       size_train = 0.5,
                       size_cal = .25,
                       cores = 32) {

  if (cf_method == "conformeR"){
    df_intervals <- conformeR_cf(sce,
                                 obs_condition,
                                 replicate_id,
                                 cell_type,
                                 spacing = 0.01,
                                 size_train = 0.5,
                                 size_cal = .25,
                                 cores = 32)
  }
  else {
    df_intervals <- conformeR_cfcausal(sce,
                                 obs_condition,
                                 replicate_id,
                                 cell_type,
                                 spacing = 0.01,
                                 size_train = 0.5,
                                 cores = 32)
  }

  fdr_res <- fdr(df_intervals,cores=cores)

  return(list(pred_intervals=df_intervals,
         conformeR_pvalues=fdr_res$pvaldf,
         gene_fdr_fun = fdr_res$gene_fdr_fun,
         comb_fdr = fdr_res$comb_fdr))
}
