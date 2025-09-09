# Subset SCE object according to a condition
subset_by_condition <- function(sce, condition, value) {
  sce[, colData(sce)[[condition]] == value]
}


# helper: convert SCE to wide tibble with conf_group
sce_to_wide_tibble <- function(sce,obs_condition) {
  meta_df <- as.data.frame(colData(sce)) %>%
    select(conf_group, all_of(obs_condition))

  expr_df <- as.data.frame(t(assay(sce)))  # transpose: rows=cells, cols=genes
  colnames(expr_df) <- rownames(sce)

  bind_cols(meta_df, expr_df) %>% as_tibble()
}

