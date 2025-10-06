#' Computation of FDR at level `cutoff` for each `conf_group`
#' @import BiocParallel
#' @import dplyr
#' @importFrom qvalue qvalue
#'
#' @param interval_df data frame of prediction intervals as output by
#'   function `conformeR_cfcausal` or `conformeR_cf`
#' @param cutoff confidence level for aggregated FDR
#' @param cores number of workers for BiocParallel job.
#'
#' @return a list with:
#'   - `comb_fdr`: tibble with gene-level combined FDR at cutoff
#'   - `group_fdr`: tibble with gene × conf_group FDR values at cutoff
#'   - `pvaldf`: the cell × gene p-values
#' @export
fdr <- function(interval_df, cutoff = 0.05, cores = 32) {
  set.seed(123)
  params <- BiocParallel::BatchtoolsParam(workers = cores)
  pvaldf <- interval_df |> interval_to_pval()
  genes <- interval_df |> pull(gene) |> unique()
  conf_group <- interval_df |> pull(conf_group) |> unique()

  gene_results <- BiocParallel::bplapply(genes, function(gene_name) { # Start looping over genes
    df_gene <- pvaldf |> filter(gene == gene_name)

    group_fdr_tbl <- dplyr::bind_rows(lapply(conf_group, function(g) { # Start looping over conformal groups
      df_group <- df_gene |> filter(conf_group == g)
      pvals <- df_group |> select(pvalue) |> unlist()

      # Fit qvalue with fallback
      qobj <- tryCatch({
        qvalue::qvalue(pvals, pfdr = TRUE, lambda = seq(0.05, 0.95, 0.05))
      }, error = function(e) {
        tryCatch({
          pi0_manual <- 1 / (1 + length(pvals))
          qvalue::qvalue(pvals, pi0 = pi0_manual)
        }, error = function(e2) NULL)
      })

      if (is.null(qobj) || is.null(qobj$lfdr)) {
        fdr_val <- 0
      } else {
        df_group$lfdr <- qobj$lfdr
        df <- stats::na.omit(data.frame(pvalue = df_group$pvalue, lfdr = df_group$lfdr))
        df <- df[order(df$pvalue), ]
        if (sum(df$pvalue < cutoff) == 0) {
          fdr_val <- 1
        } else {
          fdr_val <- mean(df$lfdr[df$pvalue < cutoff])
        }
      }

      Rg <- sum(df_group$pvalue <= cutoff)
      tibble::tibble(gene = gene_name,
                     conf_group = g,
                     fdr_at_cutoff = fdr_val,
                     Rg = Rg)
    }))

    # combined FDR
    total_R <- sum(group_fdr_tbl$Rg)
    comb_fdr_val <- if (total_R == 0) {
      mean(group_fdr_tbl$Rg * group_fdr_tbl$fdr_at_cutoff)
    } else {
      sum(group_fdr_tbl$Rg * group_fdr_tbl$fdr_at_cutoff) / total_R
    }

    list(group_fdr = group_fdr_tbl,
         comb_fdr = tibble::tibble(gene = gene_name,
                                   comb_fdr = comb_fdr_val))
  }, BPPARAM = params)

  group_fdr <- dplyr::bind_rows(lapply(gene_results, `[[`, "group_fdr"))
  comb_fdr  <- dplyr::bind_rows(lapply(gene_results, `[[`, "comb_fdr"))

  return(list(
    group_fdr = group_fdr,
    comb_fdr  = comb_fdr,
    pvaldf    = pvaldf
  ))
}



