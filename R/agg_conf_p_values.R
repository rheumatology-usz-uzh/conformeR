# Outputs the sce returned by conformeR with gene-level aggregated p-values within each conformal groups. Averaging is performed.
#
# @param sce output of conformeR with the assay "adj_pvalues".
#
# @return the sce returned by conformeR with gene-level aggregated p-values within each conformal groups.
#
agg_conf_p_values <- function(sce) {
  sce <- aggregateAcrossCells(
    sce,
    ids = sce$conf_group,
    statistics = "mean",
    use.assay.type = "adj_pvalues"
  )
  return(sce)
}
