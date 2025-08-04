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
        if (is.na(fdr_g)) fdr_g <- 1
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

      qobj <- tryCatch({
        qvalue(df_group$pvalue, pfdr = TRUE,
               lambda = seq(min(df_group$pvalue), min(max(df_group$pvalue), 0.95), 0.05))
      }, error = function(e) return(NULL))

      if (is.null(qobj)) return(approxfun(x = 0:1, y = rep(1, 2), rule = 2))

      df_group$lfdr <- qobj$lfdr
      df <- na.omit(data.frame(pvalue = df_group$pvalue, lfdr = df_group$lfdr))
      df <- df[order(df$pvalue), ]

      thresholds <- sort(unique(df$pvalue))

      cond_expect <- sapply(thresholds, function(p) {
        if (sum(df$pvalue < p) == 0) return(NA_real_)
        mean(df$lfdr[df$pvalue < p])
      })

      return(approxfun(x = thresholds, y = cond_expect,
                       method = "linear", rule = 2, ties = "ordered"))
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
