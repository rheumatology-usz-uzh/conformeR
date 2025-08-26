# Subset SCE object according to a condition
subset_by_condition <- function(sce, condition, value) {
  sce[, colData(sce)[[condition]] == value]
}

# Convert SCE object to tibble containing logcounts assay + condition
sce_to_tibble <- function(sce, condition, gene_names) {
  as_tibble(t(assay(sce, "logcounts"))) |>
    mutate(Tr = colData(sce)[[condition]]) |>
    setNames(c(gene_names, "Tr"))
}

