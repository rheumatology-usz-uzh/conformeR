#' Computation of FDR at level `cutoff` for each `conf_group`
#' @import BiocParallel
#' @import dplyr
#' @importFrom qvalue qvalue
#'
#' @param p_val_df data frame of conformal p-values as output by function `conformeR_cfcausal`
#' @param cutoff confidence level for FDR
#' @param cores number of workers for BiocParallel job.
#'
#' @return a data frame of dimension number of `conf_group` x number of genes containing FDR(level=`cutoff`) for each gene in a given  `conf_group`
#' @export

fdr <- function(p_val_df, cutoff = 0.05, cores = 32) {
  suppressMessages(suppressWarnings({
    set.seed(123)
    params <- BatchtoolsParam(workers = cores)
    results <- lapply(unique(p_val_df$conf_group), function(g) {
      gene_fdr <- bplapply(seq_len(ncol(p_val_df) - 1), function(i) {
        gene_name <- colnames(p_val_df)[i]
        df_group <- p_val_df %>%
          filter(conf_group == g) %>%
          pull(!!sym(gene_name)) %>%
          as.numeric()
        if (any(is.na(df_group))) {
          warning(paste0("NA values found in ", gene_name, " for group ", g))
          return(NA)
        }
        qobj <- qvalue(p = df_group, lambda = seq(min(df_group), min(max(df_group),0.95), 0.05))
        FDR <- ifelse(max(qobj$qvalues[qobj$pvalues <= cutoff], na.rm = TRUE)==-Inf, 1,max(qobj$qvalues[qobj$pvalues <= cutoff], na.rm = TRUE))
        return(FDR)
      }, BPPARAM = params)
      return(unlist(gene_fdr))
    })
    res_mat <- matrix(unlist(results), ncol=ncol(p_val_df)-1,byrow=T)
    colnames(res_mat) <- colnames(p_val_df)[-ncol(p_val_df)]
    rownames(res_mat) <- unique(p_val_df$conf_group)
    return(res_mat)
  }))
}
