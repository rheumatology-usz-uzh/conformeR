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

    conf_groups <- unique(p_val_df$conf_group)
    genes <- unique(p_val_df$gene)

    results <- lapply(conf_groups, function(g) {
      gene_fdr <- bplapply(genes, function(gene_name) {
        tryCatch({
          df_group <- p_val_df %>%
            filter(conf_group == g, gene == gene_name)
          pvals <- df_group$pvalue
          qobj <- qvalue(p = pvals, lambda = seq(min(pvals), min(max(pvals), 0.95), 0.05))
          fdr_val <- ifelse(
            all(qobj$pvalues > cutoff),
            1,
            max(qobj$qvalues[qobj$pvalues <= cutoff], na.rm = TRUE)
          )
          return(fdr_val)

        })
      }, BPPARAM = params)

      return(unlist(gene_fdr))
    })

    res_mat <- matrix(unlist(results), ncol = length(genes), byrow = TRUE)
    colnames(res_mat) <- genes
    rownames(res_mat) <- conf_groups

    return(res_mat)
  }))
}
