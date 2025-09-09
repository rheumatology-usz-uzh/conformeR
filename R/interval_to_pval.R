#' Convert conformal intervals to p-values
#'
#' Internal helper function that aggregates conformal prediction intervals
#' computed over a grid of confidence levels alpha for each cell x gene into one cell x gene
#' p-value computed as \deqn{(1 + \#\{ \text{intervals covering 0} \}) / (1 + length(grid))}
#'
#' @param interval_df A data frame of prediction intervals as output by function `conformeR_cfcausal` or `conformeR_cf`
#'
#' @return A data frame with one row per \code{cell_id}-\code{gene}
#'   combination, and a column \code{pvalue} giving the aggregated p-value.

interval_to_pval <- function(interval_df){
  pval_df <-interval_df |>
    group_by(cell_id,conf_group,gene) |>
    summarise(
      pvalue = (1 + sum(covered=="inside")) / (1 + n()),
      .groups = "drop"
    ) |>
    mutate(gene = gene,
           conf_group= conf_group) |>
    select(-cell_id)
  return(pval_df)
}
