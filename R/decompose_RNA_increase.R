decompoose_RNA_increase <- function(
  a_idx = '24h_1_2_sn',
  b_idx = '2h_unstim',
  source_data = rna_data_s,
  tf = function(x) x) {

  p_dat <- tibble(
    gene = rownames(source_data),
    diff = tf(source_data[, a_idx]) - tf(source_data[, b_idx]),
    a = tf(source_data[, a_idx]),
    b = tf(source_data[, b_idx])
  ) %>%
  left_join(selection_crit_table_u, by = 'gene')

  rel_size <- p_dat %>% summarize(across(c(a, b), sum))
  if (rel_size$a / rel_size$b < 1) return(NULL)

  p_dat %>%
    group_by(geneset) %>%
    summarize(diff = sum(diff)) %>%
    dplyr::mutate(rel_diff = diff / sum(diff)) %>%
    dplyr::mutate(a = a_idx) %>%
    dplyr::mutate(b = b_idx) %>%
    { . }
}

plot_diff_comp_mat <- function(
  source_data = rna_data_s,
  ref_sample = '2h_unstim') {

  p_dat <- tidyr::expand_grid(
    a = colnames(source_data),
    b = ref_sample) %>%
    dplyr::filter(a != b) %>%
    purrr::pmap_dfr(decompoose_RNA_increase)

  ggplot(p_dat, aes(x = geneset, y = rel_diff, fill = geneset)) +
    geom_col(alpha = .8) +
    facet_wrap(~a) +
    coord_flip()
}
