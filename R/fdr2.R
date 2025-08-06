

#' Computation of FDR at level `cutoff` for each `conf_group`
#' @import BiocParallel
#' @import dplyr
#' @importFrom qvalue qvalue
#' @importFrom stats approxfun
#'
#' @param p_val_df data frame of conformal p-values as output by function `conformeR_cfcausal`
#' @param cutoff confidence level for FDR
#' @param cores number of workers for BiocParallel job.
#'
#' @return a list of length `unique(p_val_df$gene)` each containing a list of length `unique(conf_group)`+1 of functions for computing tail-area FDR, and combined FDR across `unique(conf_group)` at level `cutoff`
#' @export


fdr2 <- function(p_val_df, cutoff=.05, cores=32) {
  genes <- unique(p_val_df$gene)
  set.seed(123)
  params <- BatchtoolsParam(workers = cores)

  # Loop over genes
  gene_fdr_lists <- bplapply(genes, function(gene_name) {
    df_gene <- p_val_df %>% filter(gene == gene_name)
    conf_groups <- unique(df_gene$conf_group)

    # Estimate FDR functions per group
    Fdr_list <- lapply(conf_groups, function(g) {
      df_group <- df_gene %>%
        filter(conf_group == g)
      pvals <- df_group$pvalue

      # First attempt: adaptive lambda
      qobj <- tryCatch({
        lambda_grid <- seq(min(pvals), min(max(pvals), 0.95), 0.05)
        qvalue(pvals, pfdr = TRUE, lambda = lambda_grid)
      }, error = function(e) {
        # On error, try with manually specified pi0
        tryCatch({
          pi0_manual <- 1 / (1 + length(pvals))
          qvalue(pvals, pi0 = pi0_manual)
        }, error = function(e2) {
          # On second failure, return NULL
          return(NULL)
        })
      })

      df_group$lfdr <- qobj$lfdr
      df <- na.omit(data.frame(pvalue = df_group$pvalue, lfdr = df_group$lfdr))
      df <- df[order(df$pvalue), ]

      thresholds <- sort(unique(df$pvalue))

      cond_expect <- sapply(thresholds, function(p) {
        if (sum(df$pvalue < p) == 0) return(NA_real_)
        mean(df$lfdr[df$pvalue < p])
      })

      valid_idx <- which(!is.na(cond_expect))
      if (length(valid_idx) >= 2) {
        return(approxfun(x = thresholds[valid_idx], y = cond_expect[valid_idx],
                         method = "linear", rule = 2, ties = "ordered"))
      } else {
        constant_fdr <- 1 / (1 + length(qobj$qvalues))
        return(approxfun(x = 0:1, y = rep(constant_fdr, 2), rule = 2)) # Fallback function for very small p-values
      }

    })

    names(Fdr_list) <- conf_groups

    # Compute combined FDR
    group_terms <- sapply(conf_groups, function(patient) {
      group_df <- df_gene %>% filter(conf_group == patient)
      Rg <- sum(group_df$pvalue <= cutoff)
      fdr_g <- Fdr_list[[patient]](cutoff)
      if (is.na(fdr_g)) fdr_g <- 1
      return(Rg * fdr_g)
    })
    total_R <- sum(df_gene$pvalue <= cutoff)

    Fdr_list[["comb_fdr"]] <- if (total_R == 0) 1 else sum(group_terms) / total_R

    return(Fdr_list)
  }, BPPARAM = params)

  names(gene_fdr_lists) <- genes
  return(gene_fdr_lists)
}
