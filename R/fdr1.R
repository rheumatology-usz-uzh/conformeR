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

fdr1 <- function(p_val_df, cutoff = 0.05, cores = 32) {
  suppressMessages(suppressWarnings({
    set.seed(123)
    params <- BatchtoolsParam(workers = cores)

    genes <- unique(p_val_df$gene)
    conf_groups <- unique(p_val_df$conf_group)

    # Loop over genes
    gene_results <- bplapply(genes, function(gene_name) {
      df_gene <- p_val_df %>% filter(gene == gene_name)

      # Parallel FDR computation for each conf_group
      group_fdrs <- lapply(conf_groups, function(g) {
        tryCatch({
          df_group <- df_gene %>% filter(conf_group == g)
          pvals <- df_group$pvalue
          pvals <- pvals[!is.na(pvals) & pvals >= 0 & pvals <= 1]

          if (length(pvals) < 2) return(NA_real_)

          qobj <- qvalue(p = pvals, lambda = seq(min(pvals), min(max(pvals), 0.95), 0.05))
          fdr_val <- ifelse(
            all(qobj$pvalues > cutoff),
            1,
            max(qobj$qvalues[qobj$pvalues <= cutoff], na.rm = TRUE)
          )
          return(fdr_val)
        }, error = function(e) {
          return(NA_real_)
        })
      })

      # Convert to named vector
      group_fdrs <- setNames(unlist(group_fdrs), conf_groups)

      # Compute combined FDR for this gene
      group_terms <- sapply(conf_groups, function(g) {
        df_group <- df_gene %>% filter(conf_group == g)
        Rg <- sum(df_group$pvalue <= cutoff, na.rm = TRUE)
        fdr_g <- group_fdrs[g]
        # Give value 1 everywhere if there are no discoveries. (pi0 \sim 1)
        if (is.na(fdr_g)) fdr_g <- 1
        if (is.na(Rg)) Rg <- 1
        return(Rg * fdr_g)
      })

      total_R <- sum(df_gene$pvalue <= cutoff, na.rm = TRUE)
      comb_fdr <- if (total_R == 0) 1 else sum(group_terms) / total_R

      return(c(group_fdrs, combined = comb_fdr))
    }, BPPARAM = params)

    # Build result matrix
    res_mat <- do.call(cbind, gene_results)
    rownames(res_mat) <- c(conf_groups, "combined")
    colnames(res_mat) <- genes

    return(res_mat)
  }))
}
