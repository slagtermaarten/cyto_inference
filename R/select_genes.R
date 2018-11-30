l21 <- function(x) log2(x + 1)
## Ratio with the largest of two in the denominator
mr <- function(x, y)  l21(x) - l21(y)

## Stim name look up
stim_name_lu <-
  c('tnf' = '10 ng/ml TNFa',
    'ifn' = '100 ng/ml IFNy',
    'combo' = '100 ng/ml IFNy 10 ng/ml TNFa')
stim_name_lui <- rlang::set_names(names(stim_name_lu), stim_name_lu)


render_table <- function(dtf)  {
  dtf %>%
    # dplyr::select(-message_cm, -message_cp) %>%
    knitr::kable(digits = 2) %>%
    kableExtra::kable_styling(latex_options = 'scale_down') %>%
    # row_spec(1:, color = 'white', background = 'black')
    { . }
}


prepare_max_diff_stats <- function(
  gene = 'APOBEC3G',
  min_lfc = 1,
  min_exp_diff = 10,
  epsilon = .1375035,
  lookup_data = targets::tar_read('kallisto_5029')) {

  p_dat <-
    lookup_data[which(rownames(lookup_data) == gene), ] %>%
    unlist() %>%
    tibble::enframe('sample_name', 'gexp') %>%
    # dplyr::mutate(gexp = log2(gexp + 1)) %>%
    merge_sample_annotation()

  stats <- p_dat %>%
    # dplyr::filter(sn_rank == 1) %>%
    dplyr::filter(tnf_rank %in% c(1, 4) & ifn_rank %in% c(1, 4)) %>%
    dplyr::group_by(duration) %>%
    dplyr::summarize(
      stim_group, tnf_rank, ifn_rank,
      gexp,
      abs_diff_from_us =
        gexp - gexp[stim_group == 'Unstimulated in vitro'],
      abs_diff_from_combo =
        gexp - gexp[stim_group == '100 ng/ml IFNy 10 ng/ml TNFa'],
      abs_diff_from_tnf =
        gexp - gexp[stim_group == '10 ng/ml TNFa'],
      abs_diff_from_ifn =
        gexp - gexp[stim_group == '100 ng/ml IFNy'],
      abs_diff_from_sn =
        gexp - gexp[stim_group == '1/2 SN'],
      lfc_from_us =
        mr(gexp, gexp[stim_group == 'Unstimulated in vitro']),
      lfc_from_tnf =
        mr(gexp, gexp[stim_group == '10 ng/ml TNFa']),
      lfc_from_ifn =
        mr(gexp, gexp[stim_group == '100 ng/ml IFNy']),
      lfc_from_combo =
        mr(gexp, gexp[stim_group == '100 ng/ml IFNy 10 ng/ml TNFa']),
      lfc_from_sn = mr(gexp, gexp[stim_group == '1/2 SN'])
    ) %>%
    dplyr::filter(stim_group != 'Unstimulated in vitro') %>%
    # dplyr::select(-stim_group) %>%
    # debug_pipe() %>%
    dplyr::ungroup() %>%
    # dplyr::filter(abs(diff_from_us) == max(diff_from_us)) %>%
    # pull(duration) %>%
    { . }

  # browser()
  # colSums(in_vitro_bulk_kallisto_cpm)
  # stats$lfc_from_us
  # stats$abs_diff_from_us
  # winning_stats <- stats[which.max(stats$abs_diff_from_us), ]
  # p_dat[which.max(p_dat$gexp), ]

  ## Determine winning stimulus; the stimulus that effectuates the
  ## largest fold change in comparison to unstimulated across all
  ## available time points
  sel <- tibble::tibble(ifn_rank = c(4, 4, 1), tnf_rank = c(4, 1, 4))
  winning_stats <- stats %>%
    dplyr::right_join(sel, by = c('ifn_rank', 'tnf_rank')) %>%
    dplyr::slice_max(lfc_from_us, n = 1)
  if (nrow(winning_stats) > 1) winning_stats <- winning_stats[1, ]
  if (maartenutils::eps(winning_stats$lfc_from_us, 0, 1e-3))
    return(NULL)

  if (winning_stats$stim_group == stim_name_lu[['combo']]) {
    ## Does the combo truly exceed both stimuli for all timepoints?
    ## Else, make the largest single stimulus the 'winner'
    diffs <- stats %>%
      dplyr::filter(stim_group == stim_name_lu[['combo']])

    if (all(abs(diffs$lfc_from_ifn) <= epsilon)) {
      cc <- stim_name_lu[['ifn']]
    } else if (all(abs(diffs$lfc_from_tnf) <= epsilon)) {
      cc <- stim_name_lu[['tnf']]
    } else {
      cc <- stim_name_lu[['combo']]
    }

    winning_stats <- stats %>%
      dplyr::filter(stim_group == cc) %>%
      dplyr::arrange(desc(lfc_from_us)) %>%
      dplyr::slice_max(1)
  }

  return(list(
      'stats' = stats,
      'winning_stats' = winning_stats,
      'p_dat' = p_dat)
  )
}


#' Determine for each of the stimulating signals whether the gene
#' response is monotonic over time. Return TRUE if all are monotonic,
#' FALSE otherwise
#'
extract_monotonicity <- function(
  gene = 'APOBEC3G',
  min_lfc = 1,
  min_exp_diff = 10,
  epsilon = .1375035,
  lookup_data = targets::tar_read('kallisto_5029')) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  e_stats <- prepare_max_diff_stats(gene = gene,
    min_lfc = min_lfc,
    min_exp_diff = min_exp_diff,
    epsilon = epsilon,
    lookup_data = lookup_data)

  if (is.null(e_stats)) return(NULL)

  responsive_stimuli <-
    e_stats$stats %>%
    dplyr::filter(lfc_from_us >= 1) %>%
    dplyr::pull(stim_group) %>%
    unique()

  if (length(responsive_stimuli) == 0)
    return(NULL)

  out <- responsive_stimuli %>%
    purrr::map(function(rs) {
      # qd <- e_stats$winning_stats[, 'stim_group']
      qd <- tibble('stim_group' = rs)
      v <- e_stats[['p_dat']] %>%
        dplyr::right_join(qd, by = 'stim_group') %>%
        dplyr::arrange(duration) %>%
        pull(gexp)
      length(unique(sign(diff(v)))) == 1
    }) %>%
    { list('gene' = gene, 'monotonic_response' = all(.)) }

  return(out)
}


extract_max_diff_stats <- function(
  gene = 'APOBEC3G',
  min_lfc = 1,
  min_exp_diff = 10,
  epsilon = .1375035,
  lookup_data = targets::tar_read('kallisto_5029')) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  e_stats <- prepare_max_diff_stats(gene = gene,
    min_lfc = min_lfc,
    min_exp_diff = min_exp_diff,
    epsilon = epsilon,
    lookup_data = lookup_data)

  stats <- e_stats$stats
  p_dat <- e_stats$p_dat
  winning_stats <- e_stats$winning_stats

  ## Extract stats for both stimuli at the timepoint of max diff
  tnf_stats <- stats %>% dplyr::right_join(
    tibble::tibble(
      stim_group = '10 ng/ml TNFa',
      duration = winning_stats$duration
    ), by = c('stim_group', 'duration')) %>%
    dplyr::select(-duration, -stim_group, -tnf_rank, -ifn_rank,
      -lfc_from_tnf, -abs_diff_from_tnf) %>%
    dplyr::rename_with(~paste0('tnf_', .x))

  ifn_stats <- stats %>% dplyr::right_join(
    tibble::tibble(
      stim_group = '100 ng/ml IFNy',
      duration = winning_stats$duration
    ), by = c('stim_group', 'duration')) %>%
    dplyr::select(-duration, -stim_group, -tnf_rank, -ifn_rank,
      -lfc_from_ifn, -abs_diff_from_ifn) %>%
    dplyr::rename_with(~paste0('ifn_', .x))

  ## Extract raw gene expression at max_diff timepoint
  sel <- tibble::tibble(duration = winning_stats$duration,
    ifn_rank = c(4, 4, 1), tnf_rank = c(4, 1, 4))
  exp_stats <- p_dat %>%
    dplyr::right_join(sel,
      by = c('duration', 'tnf_rank', 'ifn_rank')) %>%
    dplyr::arrange(ifn_rank, tnf_rank) %>%
    dplyr::pull(gexp) %>%
    { tibble::tibble_row('combo_exp' = .[[3]], 'ifn_exp' = .[[2]],
      'tnf_exp' = .[[1]]) } %>%
    { . }

  ## Is the expression and response strong enough?
  minimum_fc_from_us <- stats %>%
    dplyr::filter(stim_group == winning_stats$stim_group) %>%
    dplyr::pull(lfc_from_us) %>%
    { any(abs(.) >= min_lfc) }
  minimum_expression_diff <- stats %>%
    dplyr::filter(stim_group == winning_stats$stim_group) %>%
    dplyr::pull(abs_diff_from_us) %>%
    { max(.) >= min_exp_diff }

  if (winning_stats$stim_group != stim_name_lu[['combo']]) {
    ## We now determine whether this gene could be a mono reporter,
    ## i.e. a gene that responds solely to one stimulus whereas the
    ## other signal doesn't at all

    ## Does the stimulus NOT differ from combo for all timepoints?
    mono_no_diff_from_combo <- stats %>%
      dplyr::filter(stim_group == winning_stats$stim_group) %>%
      # dplyr::pull(lfc_from_combo) %>% { all(abs(.) <= epsilon) } %>%
      dplyr::summarize(
        norm_lfc_from_combo = weighted.mean(lfc_from_combo, gexp)) %>%
      { unlist(.) <= epsilon } %>%
      { . }

    ## Does the other stimulus not differ from unstimulated for all
    ## timepoints?
    oc <- setdiff(stim_name_lu[1:2], winning_stats$stim_group)
    other_ck_no_diff_from_us <- stats %>%
      dplyr::filter(stim_group == oc) %>%
      # dplyr::pull(lfc_from_us) %>%
      # { all(abs(.) <= epsilon) }
      dplyr::summarize(
        norm_lfc_from_us = weighted.mean(lfc_from_us, gexp)) %>%
      { unlist(.) <= epsilon }

    test_exp <- stats %>%
      dplyr::arrange(duration) %>%
      dplyr::filter(stim_group == oc) %>%
      pull(gexp)

    control_exp <- stats %>%
      dplyr::arrange(duration) %>%
      dplyr::filter(stim_group == winning_stats$stim_group) %>%
      pull(gexp)

    other_ck_response_factor <-
      weighted.mean(test_exp / control_exp, control_exp)

    mono <- mono_no_diff_from_combo &&
      (other_ck_no_diff_from_us || other_ck_response_factor <=
        2^epsilon)
    all_stims_not_different_from_us <- NA
    combo_different_from_us <- NA
  } else {
    ## Could this be a synergistic mono gene? For such genes, mono
    ## stimuli do not differ from unstimulated, but the combo is
    ## detectable larger than unstim for at least one timepoint
    sel <- tibble::tibble(ifn_rank = c(4, 1), tnf_rank = c(1, 4))
    all_stims_not_different_from_us <-
      stats %>% dplyr::right_join(sel) %>%
      # pull(lfc_from_us) %>%
      # { all(. <= epsilon) } %>%
      dplyr::summarize(weighted.mean(lfc_from_us, gexp)) %>%
      { unlist(.) <= epsilon } %>%
      { . }

    sel <- tibble::tibble(ifn_rank = c(4, 4), tnf_rank = c(4, 4))
    combo_different_from_us <-
      stats %>% dplyr::right_join(sel) %>%
      dplyr::pull(lfc_from_us) %>%
      { any(. > min_lfc) }
      # dplyr::summarize(weighted.mean(lfc_from_us, gexp)) %>%
      # { any(unlist(.) > min_lfc } %>%

    mono_no_diff_from_combo <- NA
    other_ck_no_diff_from_us <- NA
    other_ck_response_factor <- NA
    mono <- all_stims_not_different_from_us && combo_different_from_us
  }

  sel <- tibble::tibble(ifn_rank = 4, tnf_rank = 1)
  ifn_max_diff_from_us <- stats %>%
    dplyr::right_join(sel, by = c('ifn_rank', 'tnf_rank')) %>%
    dplyr::slice_max(lfc_from_us, n = 1) %>%
    dplyr::pull(lfc_from_us)
  ifn_max_diff_exp <- stats %>%
    dplyr::right_join(sel, by = c('ifn_rank', 'tnf_rank')) %>%
    dplyr::slice_max(lfc_from_us, n = 1) %>%
    dplyr::pull(gexp)

  sel <- tibble::tibble(ifn_rank = 1, tnf_rank = 4)
  tnf_max_diff_from_us <- stats %>%
    dplyr::right_join(sel, by = c('ifn_rank', 'tnf_rank')) %>%
    dplyr::slice_max(lfc_from_us, n = 1) %>%
    dplyr::pull(lfc_from_us)
  tnf_max_diff_exp <- stats %>%
    dplyr::right_join(sel, by = c('ifn_rank', 'tnf_rank')) %>%
    dplyr::slice_max(lfc_from_us, n = 1) %>%
    dplyr::pull(gexp)

  sel <- tibble(ifn_rank = 4, tnf_rank = 4)
  combo_max_diff_from_us <- stats %>%
    dplyr::right_join(sel, by = c('ifn_rank', 'tnf_rank')) %>%
    dplyr::slice_max(lfc_from_us, n = 1) %>%
    dplyr::pull(lfc_from_us)
  combo_max_diff_exp <- stats %>%
    dplyr::right_join(sel, by = c('ifn_rank', 'tnf_rank')) %>%
    dplyr::slice_max(lfc_from_us, n = 1) %>%
    dplyr::pull(gexp)

  max_diff_bias_A <- ifn_max_diff_exp + tnf_max_diff_exp
  max_diff_bias_M <-
    abs(ifn_max_diff_from_us) /
    (abs(ifn_max_diff_from_us) + abs(tnf_max_diff_from_us))

  w_cyto <-
    stim_name_lui[[as.character(winning_stats$stim_group)]] %>%
    str_replace('\\d ng/ml ', '')

  gene_class <- glue::glue('{w_cyto}{if_else(mono, \' mono\', \'\')}')
  ## A row of stats
  out <- tibble::tibble(
    winning_stats, tnf_stats, ifn_stats, exp_stats,
    class = as.character(gene_class),
    mono_no_diff_from_combo = mono_no_diff_from_combo,
    other_ck_no_diff_from_us = other_ck_no_diff_from_us,
    all_stims_not_different_from_us =
      all_stims_not_different_from_us,
    combo_different_from_us = combo_different_from_us,
    minimum_fc_from_us = minimum_fc_from_us,
    minimum_expression_diff = minimum_expression_diff,
    other_ck_response_factor = other_ck_response_factor,
    bias_A = max_diff_bias_A,
    bias_M = max_diff_bias_M,
    combo_from_us = combo_max_diff_from_us,
    ifn_from_us = ifn_max_diff_from_us,
    tnf_from_us = tnf_max_diff_from_us,
  )
  out <- unique(out)
  if (nrow(out) > 1) browser()

  out <- out %>%
    dplyr::mutate(tnfa_bias = tnf_exp / (tnf_exp + ifn_exp)) %>%
    dplyr::rename_with(~paste0('max_diff_', .x)) %>%
    dplyr::mutate(gene = gene)

  return(out)
}


extract_max_diff_tnfa_bias <- function(
  # gene = 'CD274',
  gene = 'EPSTI1',
  lookup_data = targets::tar_read('kallisto_5029')) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  p_dat <- lookup_data[which(rownames(lookup_data) == gene), ] %>%
    unlist() %>%
    tibble::enframe('sample_name', 'gexp') %>%
    # dplyr::mutate(gexp = log2(gexp + 1)) %>%
    merge_sample_annotation()

  out <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::filter(tnf_rank %in% c(1, 4) & ifn_rank %in% c(1, 4)) %>%
    dplyr::group_by(duration) %>%
    dplyr::summarize(
      stim_group, tnf_rank, ifn_rank,
      diff_from_us = gexp / gexp[stim_group == 'Unstimulated']) %>%
    dplyr::filter(stim_group != 'Unstimulated') %>%
    dplyr::select(-stim_group) %>%
    # debug_pipe() %>%
    reshape2::dcast(duration ~ tnf_rank + ifn_rank,
      value.var = 'diff_from_us') %>%
    # dplyr::mutate(TNFa_bias = 2^`4_1` / (2^`4_1` + 2^`1_4`)) %>%
    dplyr::mutate(TNFa_bias = `4_1` / (`4_1` + `1_4`)) %>%
    dplyr::select(duration, TNFa_bias) %>%
    dplyr::arrange(duration) %>%
    # reshape2::dcast(~ duration, value.var = 'TNFa_bias') %>%
    # { rlang::set_names(.$TNFa_bias, .$duration) } %>%
    # { tibble(.$TNFa_bias, .$duration) } %>%
    # as_tibble()
    dplyr::mutate(gene = gene) %>%
    { . }
}


make_heatmap <- function(M, value_name = 'beta') {
  library(ComplexHeatmap)
  max_val <- max(abs(as.numeric(M)))
  col_fun <- circlize::colorRamp2(c(-max_val, 0, max_val),
    c('blue', 'white', 'red'))
  Heatmap(M, col = col_fun, name = glue('Stimulus\n{value_name}'),
    show_row_dend = F, show_row_names = F)
}


inspect_gene <- function(g, presentation=T) {
  missing_genes <- setdiff(g, selection_crit_table$gene)
  if (length(missing_genes) > 0) {
    warning('Missing: ', paste(missing_genes, collapse = ', '))
  }
  g <- unique(g)
  idx <- match(g, selection_crit_table$gene, )
  idx <- setdiff(idx, NA)
  d <- selection_crit_table[idx, ] %>%
    dplyr::arrange(TNFa_bias, -synergy_score_n) %>%
    dplyr::select(-synergy_score_n) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), signif, 2)) %>%
    { . }

  if (test_rendering() && presentation) {
    d %>%
      dplyr::select(-message_cm, -message_cp) %>%
      knitr::kable(digits = 2) %>%
      kableExtra::kable_styling(latex_options = 'scale_down') %>%
      # row_spec(1:, color = 'white', background = 'black')
      { . }
  } else {
    # print(d, width = 120L, n = 2000L)
    return(d)
  }
}


#' Determine timepoint of max difference between stimulic and check
#' whether concentrations are informative with respect to gene
#' expression for this timepoint
#'
#'
determine_concentration_correctness_b1 <- function(
  gene = 'CD274',
  lookup_data = in_vitro_bulk) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  p_dat <- lookup_data[which(rownames(lookup_data) == gene), ] %>%
    unlist %>%
    tibble::enframe('sample_name', 'gexp') %>%
    merge_sample_annotation(exp = '5092')

  max_diff_from_us <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::group_by(duration) %>%
    dplyr::summarize(stim_group, duration, tnf_rank, ifn_rank,
                     diff_from_us = gexp -
                       gexp[stim_group == 'Unstimulated']) %>%
    dplyr::ungroup() %>%
    dplyr::filter(abs(diff_from_us) == max(abs(diff_from_us))) %>%
    { . }

  mrs <- c( '100 ng/ml IFNy', '10 ng/ml TNFa')

  dominant_stim <- p_dat %>%
    dplyr::filter(duration == max_diff_from_us$duration) %>%
    dplyr::filter(stim_group %in% mrs) %>%
    dplyr::filter(abs(gexp) == max(gexp)) %>%
    # dplyr::pull(ifn_rank) %>%
    dplyr::pull(stim_group) %>%
    as.character

  if (grepl('IFNy', dominant_stim)) {
    order_var <- 'ifn_rank'
  } else {
    order_var <- 'tnf_rank'
  }
  other_cytokine <- setdiff(c('ifn_rank', 'tnf_rank'), order_var)

  sample_ordering <- p_dat %>%
    dplyr::filter(duration == max_diff_from_us$duration) %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::filter(.data[[other_cytokine]] == 1) %>%
    dplyr::arrange(gexp) %>%
    dplyr::pull(.data[[order_var]])

  final_test <-
    all(sample_ordering == 1:4) || all(sample_ordering == 4:1)
  list('gene' = gene, consistent_concentration_ordering = final_test)
}


tpi_to_tp <- function(i) {
  levels(sample_annotation$duration)[i] %>% as.numeric()
}


#' Determine timepoint of max difference between stimulic and check
#' whether gene expression is informative with respect to
#' concentration differences at this timepoint
#'
#'
determine_concentration_correctness <- function(
  gene = 'CD274',
  lookup_data = targets::tar_read('kallisto_5029')) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  p_dat <- lookup_data[which(rownames(lookup_data) == gene), ] %>%
    unlist() %>%
    tibble::enframe('sample_name', 'gexp') %>%
    merge_sample_annotation()

  ## Identify the biggest change from unstim to later identify the
  ## dominant stimulus
  max_diff_from_us <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::filter(!(tnf_rank == 4 & ifn_rank == 4)) %>%
    dplyr::group_by(duration) %>%
    dplyr::summarize(stim_group, duration, tnf_rank, ifn_rank,
                     diff_from_us = gexp -
                       gexp[stim_group == 'Unstimulated'],
                     .groups = 'drop') %>%
    dplyr::ungroup() %>%
    dplyr::filter(abs(diff_from_us) == max(abs(diff_from_us))) %>%
    { . }

  if (nrow(max_diff_from_us) > 1) {
    return(list(
        'gene' = gene,
        message = 'multiple_samples_max_diff',
        n_congruent_timepoints = NA_integer_)
    )
  }

  mrs <- c(
    '100 ng/ml IFNy',
    '10 ng/ml TNFa'
  )

  if (max_diff_from_us$tnf_rank < 4 &&
      max_diff_from_us$ifn_rank < 4) {
    return(list(
        'gene' = gene,
        message = 'max_GE_at_lower_concentration',
        n_congruent_timepoints = NA_integer_)
    )
  }

  dominant_stim <- p_dat %>%
    dplyr::filter(duration == max_diff_from_us$duration) %>%
    dplyr::filter(stim_group %in% mrs) %>%
    dplyr::filter(abs(gexp) == max(gexp)) %>%
    # { . }
    # dplyr::pull(ifn_rank) %>%
    dplyr::pull(stim_group) %>%
    as.character

  if (length(dominant_stim) > 1) browser()

  if (grepl('IFNy', dominant_stim)) {
    order_var <- 'ifn_rank'
  } else {
    order_var <- 'tnf_rank'
  }
  other_cytokine <- setdiff(c('ifn_rank', 'tnf_rank'), order_var)

  ## Exclude the lowest non-zero concentration from sample orderings,
  ## as this one turned out to be too noisy for our dataset
  sample_ordering <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::filter(.data[[other_cytokine]] == 1) %>%
    dplyr::filter(.data[[order_var]] != 1) %>%
    dplyr::group_by(duration) %>%
    dplyr::arrange(gexp, .by_group=T) %>%
    dplyr::summarise(ordering =
      paste(.data[[order_var]], collapse = ','), .groups = 'drop') %>%
    # dplyr::ungroup() %>%
    # dplyr::pull(.data[[order_var]])
    # dplyr::select(duration, .data[[order_var]]) %>%
    { . }

  ## The maximum
  n_congruent_timepoints <-
    table(sample_ordering$ordering) %>%
    as.numeric %>%
    max()

  list('gene' = gene,
    message = 'OK',
    n_congruent_timepoints = n_congruent_timepoints)
}


#' Determine timepoint of max difference between stimulic and check
#' whether gene expression is informative with respect to
#' concentration differences at this timepoint
#'
#'
nrow(targets::tar_read('kallisto_5029'))
determine_stimulus_correlation <- function(
  gene = 'CD274', lookup_data = targets::tar_read('kallisto_5029')) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  ## Exclude the lowest non-zero concentration from sample orderings,
  ## as this one turned out to be too noisy for our dataset
  p_dat <- lookup_data[which(rownames(lookup_data) == gene), ] %>%
    unlist %>%
    tibble::enframe('sample_name', 'gexp') %>%
    merge_sample_annotation(exp = '5092') %>%
    dplyr::filter(tnf_rank != 2 & ifn_rank != 2 & sn_rank != 2) %>%
    dplyr::arrange(duration, tnf_rank, ifn_rank, sn_rank) %>%
    { . }

  time_levs <- setdiff(levels(p_dat$duration), c('Unknown', NA))
  inter_timepoint_corrs <-
    map(auto_name(1:(length(time_levs)-1)), function(i) {
    a1 <- p_dat %>% dplyr::filter(duration == time_levs[i]) %>%
      pull(gexp)
    a2 <- p_dat %>% dplyr::filter(duration == time_levs[i+1]) %>%
      pull(gexp)
    if (all(a1 == 0) | all(a2 == 0)) return(0)
    suppressWarnings(cor(a1, a2, method = 'pearson'))
  }) %>% { . }

  extra_stats <-
    auto_name(c('mean', 'median', 'min', 'max')) %>%
    map(function(f) get(f)(unlist(inter_timepoint_corrs)))

  out <- c(inter_timepoint_corrs, extra_stats) %>%
    { rlang::set_names(., paste0('timecor_', names(.))) } %>%
    { c(list('gene' = gene), .) }

  return(out)
}


#' Determine stimulus associated with highest (lowest) gene
#' expression for each timepoint individually and score how consistent
#' that is between timepoints.
#'
#'
determine_combo_maximality <- function(
  gene = 'NFIX',
  lookup_data = targets::tar_read('kallisto_5029')) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  p_dat <- lookup_data[which(rownames(lookup_data) == gene), ] %>%
    unlist %>%
    tibble::enframe('sample_name', 'gexp') %>%
    merge_sample_annotation()

  ## Minimum gexp
  min_diff_from_us <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::group_by(duration) %>%
    dplyr::arrange(gexp, .by_group=T) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(stim_group) %>%
    as.character()

  ## Maximum gexp
  max_diff_from_us <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::group_by(duration) %>%
    dplyr::arrange(gexp, .by_group=T) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::pull(stim_group) %>%
    as.character()

  min_max_congruence <- max(table(min_diff_from_us))
  max_max_congruence <- max(table(max_diff_from_us))

  list('gene' = gene,
    message = 'OK',
    min_max_congruence = min_max_congruence,
    max_max_congruence = max_max_congruence)
}


#' Determine stimulus associated with highest gene
#' expression diff for a particular gene and find out whether that
#' maximal diff is acquired relatively early or late
#'
#'
determine_combo_maximality <- function(
  gene = 'NFIX',
  lookup_data = in_vitro_bulk) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  p_dat <- lookup_data[which(rownames(lookup_data) == gene), ] %>%
    unlist() %>%
    tibble::enframe('sample_name', 'gexp') %>%
    merge_sample_annotation(exp = '5092')

  ## Minimum gexp
  min_diff_from_us <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::group_by(duration) %>%
    dplyr::arrange(gexp, .by_group=T) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(stim_group) %>%
    as.character()

  ## Maximum gexp
  max_diff_from_us <- p_dat %>%
    dplyr::filter(sn_rank == 1) %>%
    dplyr::group_by(duration) %>%
    dplyr::arrange(gexp, .by_group=T) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::pull(stim_group) %>%
    as.character()

  min_max_congruence <- max(table(min_diff_from_us))
  max_max_congruence <- max(table(max_diff_from_us))

  list('gene' = gene,
    message = 'OK',
    min_max_congruence = min_max_congruence,
    max_max_congruence = max_max_congruence)
}


fit_forest <- function(all_dat,
  exp_vars = c('max_t', 'max_beta', 'max_diff_from_us',
    'n_congruent_timepoints', 'median_diff_from_us',
    'min_diff_from_us', 'Amax', 'min_max_congruence',
    'max_max_congruence', 'TNFa_bias', 'timecor_median',
    'max_diff_TNFa_bias', 'timecor_max', 'timecor_min')) {


  library(tidymodels)
  all_dat$geneset[all_dat$geneset == 'dubious'] <- 'blacklist'
  missing_vars <- setdiff(exp_vars, colnames(all_dat))
  if (length(missing_vars) > 0) {
    message('Missing these vars from data:',
      paste0(missing_vars, collapse = ', '))
  }
  exp_vars <- intersect(exp_vars, colnames(all_dat))

  all_dat <- all_dat[complete.cases(all_dat[, exp_vars]), ]

  t_dat <- all_dat %>% dplyr::filter(geneset != 'none')
  t_dat$geneset <- droplevels(t_dat$geneset)
  t_dat_split <- initial_split(t_dat, prop = 3/4)
  t_dat_cv <- vfold_cv(training(t_dat_split), v = 3)

  c_formula <-
    sprintf('geneset ~ %s', paste(exp_vars, collapse = ' + ')) %>%
    as.formula()

  rf_model <-
    rand_forest() %>%
    set_args(mtry = tune()) %>%
    set_engine('ranger', importance = 'impurity') %>%
    set_mode('classification')

  t_recipe <- recipe(c_formula, data = t_dat)

  rf_workflow <- workflow() %>%
    add_recipe(t_recipe) %>%
    add_model(rf_model)

  param_final <- rf_workflow %>%
    tune_grid(
      resamples = t_dat_cv,
      grid = tidyr::expand_grid(mtry = 3:length(exp_vars))
    ) %>%
    # select_best(metric = bal_accuracy())
    select_best()
    # select_best(metric = yardstick::bal_accuracy())

  rf_workflow <- rf_workflow %>%
    finalize_workflow(param_final)

  rf_fit <- rf_workflow %>%
    last_fit(t_dat_split)

  ACC <- rf_fit %>%
    collect_metrics() %>%
    # filter(.metric == 'accuracy') %>%
    # pull(.estimate)
    { . }

  final_model <- fit(rf_workflow, t_dat)
  ranger_obj <- pull_workflow_fit(final_model)$fit
  var_imp <-
    ranger_obj$variable.importance %>%
    tibble::enframe('term', 'imp') %>%
    dplyr::mutate(norm_imp = imp / sum(imp, na.rm = T)) %>%
    { . }

  # model_preds <-
  #   tibble(gene = all_dat$gene) %>%
  #   bind_cols(
  #     predict(final_model, new_data = all_dat)
  #   ) %>%
  #   dplyr::rename(geneset_predicted = `.pred_class`)

  pred_classes <- setdiff(levels(all_dat$geneset), 'none')
  probs <- predict(ranger_obj, data = all_dat)$predictions %>%
    set_colnames(pred_classes)
  max_classes <- pred_classes[apply(probs, 1, which.max)]
  # stopifnot(all(model_preds$geneset_predicted == max_classes))
  model_preds <- tibble(
    'gene' = all_dat$gene,
    'geneset_predicted' = max_classes,
    'confidence' = apply(probs, 1, max))

  return(list(
      mod = final_model, accuracy = ACC,
      var_imp = var_imp, preds = model_preds))

  return(model_preds)
}

recode_subs <- function(v) {
  subs_table <- c(
    'MS' = 'mono synergy',
    'SM' = 'mono synergy',
    'M' = 'mono reporter',
    'MR' = 'mono reporter',
    'B' = 'blacklist',
    'W' = 'whitelist',
    'I' = 'whitelist',
    'L' = 'lowly expressed',
    'S' = 'synergy',
    'AS' = 'anti-synergy')
  if (length(v) == 0) return(v)
  tryCatch(return(plyr::revalue(v, subs_table, warn_missing = F)),
    error = function(e) { v })
}


#' Update the selection_crit_table geneset column with ML-informed
#' (curated) geneset info
#'
#'
update_selection_crit_table_gs <- function(
  selection_crit_table, version_idx = 6) {

  selection_crit_table_u <- selection_crit_table
  selection_crit_table_u$geneset <-
    selection_crit_table_u$geneset %>%
    forcats::fct_recode(
      'mono reporter' = 'expected_monoreporter_curated',
      'blacklist' = 'dubious'
    ) %>%
    forcats::fct_expand('synergy',
      'mono synergy', 'anti-synergy', 'lowly expressed')

  expected_levs <- levels(selection_crit_table_u$geneset)

  for (c_version_idx in 0:(version_idx-1)) {
    fo <- maartenutils::gen_file_overview(
      file.path(data_raw_dir, 'gene_lists', 'corrections'),
      pat = glue('v{c_version_idx}'),
      include_full = T) %>%
      dplyr::arrange(short_fn)
    for (corr_fn in fo$full_fn) {
      if (!file.exists(corr_fn)) next()
      fh <- readr::read_tsv(corr_fn, show_col_types = F)
      if (ncol(fh) != 3) browser()
      if (nrow(fh) == 0) next()
      # if ('W' %in% fh[[3]]) browser()
      # print(table(fh[[3]]))
      cn <- grep('geneset_corrected', colnames(fh), value = T)
      if (length(cn) == 0) browser()
      fh[[cn]] <- recode_subs(fh[[cn]])
      ## Column name of corrected labels
      while (any(fh[[cn]] %nin% expected_levs)) {
        setdiff(fh$geneset_corrected,
          levels(selection_crit_table_u$geneset))
        # fh$geneset_corrected
        browser()
      }
      gene_idxs <- match(fh$gene, selection_crit_table_u$gene)
      stopifnot(fh$gene == selection_crit_table_u$gene[gene_idxs])
      selection_crit_table_u$geneset[gene_idxs] <- fh[[cn]]
      print(basename(corr_fn))
      print(table(selection_crit_table_u$geneset))
    }
  }

  return(selection_crit_table_u)
}


gs_cols <- rlang::set_names(
  c('grey80', 'green', 'forestgreen', 'forestgreen',
    'orange', 'indianred3',
    # 'purple', 'darkpurple', 'darkpurple', 'lightpurple',
    'purple', maartenutils::darken('purple'),
    maartenutils::darken('purple'), maartenutils::lighten('purple'),
    'blue'),
  c('none', 'whitelist', 'mono reporter',
    'expected_monoreporter_curated', 'dubious', 'blacklist',
    'synergy', 'synergy mono', 'mono synergy', 'anti-synergy',
    'lowly expressed')
)


genestat_labels <- rlang::set_names(
  c(
    'maximum FC between\nstim. and unstim. sample', 
    'maximum FC between\nIFNy stim. and unstim. sample', 
    'maximum FC between\nTNFa stim. and unstim. sample', 
    'median correlation\nbetween timepoints across stimuli', 
    'log maximal expression',
    'TNFa vs. IFNy responsiveness bias\n(1 := only responsive to TNFa,\n0 := only responsive to IFNy)',
    'Maximum limma t-statistic across stimuli',
    'Maximum limma beta-coefficient across stimuli'
  ),
  c('max_diff_from_us', 
    'max_diff_ifn_from_us', 
    'max_diff_tnf_from_us', 
    'timecor_median', 'Amax', 
  'TNFa_bias', 'max_t', 'max_beta')
)

plot_genestats <- function(
  x_var = 'max_t',
  y_var = 'max_diff_from_us',
  source_dat = selection_crit_table,
  hl_genes = NULL) {

  stopifnot(all(
      levels(selection_crit_table_u$geneset) %in% names(gs_cols)))
  gs_cols <- gs_cols[levels(source_dat$geneset)]

  hl_dat <- source_dat %>%
    dplyr::filter(gene %in% hl_genes)
  # hl_dat <- NULL
  pacman::p_load('ggrastr')

  p <- source_dat %>%
    dplyr::arrange(geneset) %>%
    ggplot(aes_string(x = x_var, y = y_var,
        label = 'gene', colour = 'geneset')) +
      scale_colour_manual(values = gs_cols) +
      ggrastr::rasterise(geom_point(alpha = .5, size = .2), dpi = 500) +
      xlab(genestat_labels[x_var]) +
      ylab(genestat_labels[y_var]) +
      guides(
        colour = guide_legend(
          override.aes = list(size = 2, alpha = 1)
        )
      )

  if (F) {
    p <- add_marginal_density(p, group_color = T)
  }

  if (nrow(hl_dat) > 0) {
    p <- p + geom_text(data = hl_dat, show.legend=F)
  }

  return(p)
}

plot_gene_stats_boxplot <- function(
  y_var = 'max_t',
  source_dat = selection_crit_table_u,
  hl_genes = NULL) {

  p <- source_dat %>%
    ggplot(aes_string(
        x = 'geneset', 
        # fill = 'geneset',
        y = y_var)) +
      # scale_fill_manual(values = gs_cols) +
      # geom_boxplot(alpha = .5, outlier.shape = NA, fill = 'goldenrod3') +
      geom_violin(alpha = .5, outlier.shape = NA, 
        draw_quantiles = c(0.25, .5, .75), fill = 'goldenrod3') +
      xlab('Gene class') +
      rotate_x_labels(45) +
      ylab(genestat_labels[y_var]) +
      guides(
        colour = guide_legend(
          override.aes = list(size = 2, alpha = 1)
        )
      )

  return(p)
}


fit_limma <- function(plot_panel = F) {
  pacman::p_load('limma')
  pacman::p_load('edgeR')
  pacman::p_load('Glimma')
  ## Use unnormalized counts here
  d0 <- DGEList(counts = tar_read('kallisto_5029'))
  d0 <- calcNormFactors(d0)
  sa <- tar_read(sample_annotation_exp5029)
  stopifnot(sa$sample_name ==
    colnames(tar_read('kallisto_5029')))
  cutoff <- 1
  gene_max_e <- apply(edgeR::cpm(d0), 1, max)
  gene_max_e_rc <- apply(tar_read('kallisto_5029'), 1, max)
  # cor(gene_max_e, gene_max_e_rc)
  drop <- which(gene_max_e < cutoff)
  mean(read_geneset('whitelist') %in% names(drop))
  mean(read_geneset('blacklist') %in% names(drop))
  d <- d0[-drop, ]
  group <- with(sa, interaction(stim_group, duration))
  sa$duration <- droplevels(sa$duration)
  # sa$stim_group
  mm <- model.matrix(~1 + stim_group + duration, data = sa)
  # mm <- model.matrix(~1 + stim_group * duration, data = sa)
  ## Leave out the unstimulated sample, such that the intercept of the
  ## model will mostly describe the unstimulated sample
  mm <- mm[, !grepl('Unstimulated', colnames(mm))]
  ## Infer the global relationship between gene expression mean and
  ## variance. Genes that fall far from this 'average' line, will have
  ## increasing statistical power and lowering false positive
  ## observations.
  v <- limma::voom(d0, mm, plot=plot_panel)
  vfit <- lmFit(v, mm)
  efit <- eBayes(vfit)

  if (F) {
    B <- efit[['coefficients']]
    B_R <- B['EPSTI1', ]
    B_R <- B['IFIT5', ]
    B_R <- B['CTSL', ]
    B_R <- B['HLA-DRB1', ]
    B_R <- B['LAP3', ]
    B_R <- B['BATF2', ]

    B_R['stim_group10 ng/ml TNFa'] /
    (B_R['stim_group100 ng/ml IFNy'] + B_R['stim_group10 ng/ml TNFa'])

    exp(B_R['stim_group10 ng/ml TNFa']) /
    (exp(B_R['stim_group100 ng/ml IFNy']) +
      exp(B_R['stim_group10 ng/ml TNFa']))

    exp(abs(B_R['stim_group10 ng/ml TNFa'])) /
    (exp(abs(B_R['stim_group100 ng/ml IFNy'])) +
      exp(abs(B_R['stim_group10 ng/ml TNFa'])))
  }

  return(efit)
}

#' Load/compile gene selection criteria data.frame
#'
#' @param unconvincing_lfc The minimum fold difference required for a
#' difference in gene expression to be called 'substantial'
load_gs_data <- function(
  version_idx = 14, redo = F,
  unconvincing_lfc = round(log2(1.25), 3)) {

  conflicted::conflict_prefer('exprs', 'rlang')

  rc <- maartenutils::gen_time_comparator(
    '2021-09-18 08:30', verbose = T)

  if (!exists('exp5310_sc')) {
    exp5310_sc <- readRDS(file.path(rds_dir, 'exp5310_sc-proc.rds'))
    exp5310_sc_split <- exp5310_sc %>%
      SplitObject(split.by = 'group_id') %>%
      purrr::map(~as.SingleCellExperiment(.x)) %>%
      purrr::map(~scater::logNormCounts(.x)) %>%
      { . }
  }

  if (F && !exists('df.list')) {
    df.list <- c(list(
      'bulk' = log(as.matrix(tar_read('kallisto_5029')) + 1)),
      exp5310_sc_split)
  }

  if (F && !exists('informative_features')) {
    dec <- map(df.list, scran::modelGeneVar)
    informative_features <- dec %>%
      map(~dplyr::filter(as.data.frame(.x), p.value <= .05)) %>%
      map(~rownames(.x))
  }

  ## 5310 gene detection rate
  compute_GDR <- function(so) {
    apply(so[['RNA']]@data, 1, function(x) mean(x > 0))
  }

  in_vitro_means <-
    subset(exp5310_sc, group_id == 'sc-in_vitro')[['RNA']] %>%
    rowMeans %>%
    tibble::enframe(name = 'gene', value = 'vitro_mean')
  in_vivo_means <-
    subset(exp5310_sc, group_id == 'sc-in_vivo')[['RNA']] %>%
    rowMeans %>%
    tibble::enframe(name = 'gene', value = 'vivo_mean')

  GDR <- compute_GDR(exp5310_sc)
  in_vitro_GDR <-
    compute_GDR(subset(exp5310_sc, group_id == 'sc-in_vitro')) %>%
    tibble::enframe(name = 'gene', value = 'vitro_GDR')
  in_vivo_GDR <-
    compute_GDR(subset(exp5310_sc, group_id == 'sc-in_vivo')) %>%
    tibble::enframe(name = 'gene', value = 'vivo_GDR')

  o_fn <- file.path(data_raw_dir, 'ge_combo_maximality.tsv')
  if (rc(o_fn) || redo) {
    features <- rownames(targets::tar_read('kallisto_5029'))
    combo_max_GE_test <-
      furrr::future_map_dfr(features, determine_combo_maximality,
        lookup_data = targets::tar_read('kallisto_5029'))
    readr::write_tsv(combo_max_GE_test, o_fn)
  } else {
    combo_max_GE_test <- readr::read_tsv(o_fn, show_col_types = F)
  }

  o_fn <- file.path(data_raw_dir,
                    'gene_behavior_concentration_consistency.tsv')
  if (rc(o_fn) || redo) {
    features <- rownames(targets::tar_read('kallisto_5029'))
    concentration_consistency_patterns <-
      furrr::future_map_dfr(features,
        determine_concentration_correctness)
    readr::write_tsv(concentration_consistency_patterns, o_fn)
  } else {
    concentration_consistency_patterns <- readr::read_tsv(o_fn,
      show_col_types = F)
  }

  if (F) {
    ## This file is generated in select_genes_new.Rmd
    o_fn <- file.path(data_raw_dir,
      'human_gene_behavior_in_response_to_cytokines.tsv')
    feature_patterns <- readr::read_tsv(o_fn, show_col_types = F)

    overrep_genes <- feature_patterns %>%
      group_by(gene, stim) %>%
      dplyr::summarize(n = n()) %>%
      dplyr::filter(n > 3)

    ## All of the 'weird' genes are just zero, let's discard them
    feature_patterns %>%
      right_join(dplyr::select(overrep_genes, gene, stim)) %>%
      pull(diff_from_us) %>%
      unique %>%
      { stopifnot(length(.) == 1) }

    correct_genes <- feature_patterns %>%
      group_by(gene, stim) %>%
      dplyr::summarize(n = n()) %>%
      dplyr::filter(n <= 1)

    feature_patterns <-
      feature_patterns %>%
      right_join(dplyr::select(correct_genes, gene, stim)) %>%
      { . }

    pmedian <- function(...) apply(cbind(...), 1, median, na.rm = T)
    feature_patterns_wide <- feature_patterns %>%
      dcast(gene ~ stim, value.var = 'diff_from_us') %>%
      dplyr::mutate(min_diff_from_us = pmin(comb, IFNy, TNFa)) %>%
      dplyr::mutate(median_diff_from_us = pmedian(comb, IFNy, TNFa)) %>%
      dplyr::mutate(max_diff_from_us = pmax(comb, IFNy, TNFa)) %>%
      dplyr::mutate(
        diff_max_mono_reporter =
          max_diff_from_us > 2 &
          min_diff_from_us < .5 &
          maartenutils::eps(max_diff_from_us, median_diff_from_us, .2),
        diff_max_mono_reporter =
          ifelse(diff_max_mono_reporter,
            ifelse(abs(IFNy) > abs(TNFa), 'IFNy', 'TNFa'), 'none'),
        diff_max_mono_reporter =
          ifelse(is.na(diff_max_mono_reporter), 'none',
            diff_max_mono_reporter),
        max_diff_TNFa_bias = abs(TNFa) / (abs(TNFa) + abs(IFNy))
      ) %>%
      dplyr::rename(
        comb_max_diff_from_us = comb,
        IFNy_max_diff_from_us = IFNy,
        TNFa_max_diff_from_us = TNFa
      ) %>%
      dplyr::arrange(max_diff_TNFa_bias) %>%
      { . }
  }

  o_fn <- file.path(rds_dir, 'limma_5029_efit.rds')
  if (rc(o_fn) || redo) {
    efit <- fit_limma()
    saveRDS(efit, o_fn)
  } else {
    efit <- readRDS(o_fn)
  }

  beta_M <- efit$coefficients
  beta_M <- beta_M[, grep('stim_group', colnames(beta_M))]
  colnames(beta_M) <- gsub('stim_group', '', colnames(beta_M))
  beta_M <- beta_M

  t_M <- efit$t
  t_M <- t_M[, grep('stim_group', colnames(t_M))]
  colnames(t_M) <- gsub('stim_group', '', colnames(t_M))
  t_M <- t_M
  sigma_table <- 
    rlang::set_names(efit$s2.post, rownames(efit$lods)) %>%
    tibble::enframe(name = 'gene', value = 'sigma') %>%
    { . }

  ## Inclusion not based on SN samples
  max_t <- apply(t_M[, !grepl('SN', colnames(t_M))], 1,
    function(x) max(abs(x)))
  max_betas <- apply(beta_M[, !grepl('SN', colnames(beta_M))], 1,
    function(x) max(abs(x)))

  o_fn <- file.path(data_raw_dir,
                    'determine_stimulus_correlation.tsv')
  if (rc(o_fn) || redo) {
    features <- rownames(targets::tar_read('kallisto_5029'))
    stimulus_correlation <-
      furrr::future_map_dfr(features, determine_stimulus_correlation)

    # NA_idxs <- which(is.na(stimulus_correlation), arr.ind = T)
    # for (i in seq(nrow(NA_idxs))) {
    #   stimulus_correlation[NA_idxs[i, 1], NA_idxs[i, 2]] <- 0
    # }
    readr::write_tsv(stimulus_correlation, o_fn)
  } else {
    stimulus_correlation <- readr::read_tsv(o_fn, show_col_types = F)
  }

  o_fn <- file.path(data_raw_dir, 'max_diff_tnfa_bias.tsv')
  # file.remove(o_fn)
  if (rc(o_fn) || redo) {
    features <- rownames(targets::tar_read('kallisto_5029'))
    max_diff_tnfa_bias <- furrr::future_map_dfr(features,
      extract_max_diff_tnfa_bias)
    max_diff_tnfa_bias <- max_diff_tnfa_bias %>%
      reshape2::dcast(gene ~ duration,
        value.var = 'TNFa_bias') %>%
        dplyr::rename_with(
          ~paste0('max_diff_TNFa_bias', .x),
          matches('\\d+'))
    readr::write_tsv(max_diff_tnfa_bias, o_fn)
  } else {
    max_diff_tnfa_bias <- readr::read_tsv(o_fn, show_col_types = F)
    rownames(max_diff_tnfa_bias) <- max_diff_tnfa_bias$gene
    # max_diff_tnfa_bias['EPSTI1', ]
  }

  o_fn <- file.path(data_raw_dir, 'tp_max_diff_from_us.tsv')
  # file.remove(o_fn)
  if (rc(o_fn) || redo) {
    features <- rownames(targets::tar_read('kallisto_5029'))
    tp_max_diff_from_us <- furrr::future_map_dfr(features,
      extract_max_diff_stats)
    readr::write_tsv(tp_max_diff_from_us, o_fn)
  } else {
    tp_max_diff_from_us <- readr::read_tsv(o_fn, show_col_types = F)
  }

  # all_consistent_genes <- gs_data %>%
  #   dplyr::filter(timecor_median >= .75) %>%
  #   dplyr::filter(Amax >= 3) %>%
  #   pull(gene)

  features <- rownames(targets::tar_read('kallisto_5029'))
  if (T) {
    o_fn <- file.path(data_raw_dir,
      glue('time_informativeness{unconvincing_lfc}.tsv'))
    if (rc(o_fn) || redo) {
      # features <- all_consistent_genes
      time_informativeness <-
        furrr::future_map_dfr(features,
          determine_time_informativeness)
      readr::write_tsv(time_informativeness, o_fn)
    } else {
      time_informativeness <- readr::read_tsv(o_fn, show_col_types = F)
    }
  } else {
    time_informativeness <- tibble('gene' = features)
  }

  if (TRUE) {
    o_fn <- file.path(data_raw_dir,
      glue('response_gene_monotonicity.tsv'))
    # file.remove(o_fn)
    if (rc(o_fn) || redo) {
      # features <- all_consistent_genes
      reporter_monotonicity <-
        furrr::future_map_dfr(features, extract_monotonicity)
          # lookup_data = targets::tar_read('kallisto_5029'))
      readr::write_tsv(reporter_monotonicity, o_fn)
    } else {
      reporter_monotonicity <-
        readr::read_tsv(o_fn, show_col_types = F)
    }
  } else {
    reporter_monotonicity <- tibble('gene' = features)
  }


  ## Form a table of potential selection criteria and the TNFa_bias
  ## characteristic, which is a score for how strong the gene responds
  ## to TNFa in comparison to IFNy. Values close to 1 (0) indicate a
  ## gene that mostly responds to TNFa (IFNy).
  in_vitro_bulk_kallisto_cpm <- targets::tar_read('kallisto_5029')
  selection_crit_table <-
    tibble(
      'gene' = names(max_t),
      'geneset' = 'none',
      included = F,
      'bulk_detectable' = names(max_t) %nin% names(drop),
      'max_t' = max_t,
      'max_beta' = max_betas,
      'sigma2' = sigma_table$sigma,
      'Amean' = efit$Amean[names(max_t)],
      'Amax' = apply(log2(in_vitro_bulk_kallisto_cpm + 1), 1, max) %>%
        { .[names(max_t)] },
      'Amin' = apply(log2(in_vitro_bulk_kallisto_cpm + 1), 1, min) %>%
        { .[names(max_t)] },
      'Arange' = apply(in_vitro_bulk_kallisto_cpm, 1,
        function(x) diff(range(x))) %>%
        { .[names(max_t)] },
      'vitro_GDR' = unlist(in_vitro_GDR[match(names(max_t),
          in_vitro_GDR$gene),
        'vitro_GDR']),
      'vivo_GDR' = unlist(in_vivo_GDR[match(names(max_t),
          in_vivo_GDR$gene),
    'vivo_GDR']),
      'synergy_score' =
        with(as.data.frame(beta_M),
              (`100 ng/ml IFNy 10 ng/ml TNFa` -
               `100 ng/ml IFNy` -
               `10 ng/ml TNFa`) /
              (`100 ng/ml IFNy` + `10 ng/ml TNFa`)),
      # 'TNFa_bias' =
      #   ifelse(synergy_score > .2, NA_real_,
      #     with(as.data.frame(beta_M),
      #         `10 ng/ml TNFa` /
      #         (`100 ng/ml IFNy` + `10 ng/ml TNFa`))),
      'synergy_score_n' = ifelse(synergy_score > 1, synergy_score,
        1/synergy_score),
      'TNFa_bias' =
          with(as.data.frame(beta_M),
              abs(`10 ng/ml TNFa`) /
              (abs(`100 ng/ml IFNy`) + abs(`10 ng/ml TNFa`)))
      # 'IFNy_bias' = 1 - TNFa_bias
    ) %>%
    left_join(
      dplyr::rename(concentration_consistency_patterns,
        message_cp = message),
      by = 'gene') %>%
    left_join(
      dplyr::rename(combo_max_GE_test, message_cm = message),
      by = 'gene') %>%
    # left_join(feature_patterns_wide, by = 'gene') %>%
    # dplyr::mutate(TNFa_bias_max_diff =
    #   TNFa_max_diff_from_us /
    #   (IFNy_max_diff_from_us + TNFa_max_diff_from_us)) %>%
    # dplyr::mutate(synergy_score_max_diff =
    #   (comb_max_diff_from_us - IFNy_max_diff_from_us -
    #     TNFa_max_diff_from_us) /
    #   (IFNy_max_diff_from_us + TNFa_max_diff_from_us)) %>%
    left_join(stimulus_correlation, by = 'gene') %>%
    left_join(tp_max_diff_from_us, by = 'gene') %>%
    left_join(max_diff_tnfa_bias, by = 'gene') %>%
    left_join(time_informativeness, by = 'gene') %>%
    left_join(reporter_monotonicity, by = 'gene') %>%
    { . }

  ## V1 of manually selected genes
  genesets <- c('whitelist', 'expected_monoreporter_curated',
    'synergy', 'dubious', 'blacklist')
  for (gs in genesets) {
    idxs <- selection_crit_table$gene %in% read_geneset(gs)
    selection_crit_table$geneset[idxs] <- gs
  }
  selection_crit_table$geneset <-
    selection_crit_table$geneset %>%
    factor(levels = c('none', genesets))

  selection_crit_table$geneset <- selection_crit_table$geneset
  selection_crit_table$max_diff_from_us <-
    with(selection_crit_table,
      pmax(max_diff_ifn_from_us, max_diff_tnf_from_us))
  selection_crit_table$min_diff_from_us <-
    with(selection_crit_table,
      pmin(max_diff_ifn_from_us, max_diff_tnf_from_us))

  if (F) {
    apply_hard_thresholding()
  }

  return(selection_crit_table)
}


apply_hard_thresholding <- function() {
  ## Filtering thresholds used below
  min_t <- 7
  min_t <- 5
  min_t <- 3.2
  min_beta <- .6
  min_Amax <- 2
  min_max_diff_from_us <- 2
  min_max_diff_from_us <- 1.4
  # min_GDR <- .01
  min_GDR <- .0001

  dplyr::filter(selection_crit_table,
    geneset %in% c('whitelist', 'expected_monoreporter_curated')) %>%
    dplyr::summarize(across(c(max_beta, max_t, Amax, max_diff_from_us), min))

  ## The set of rules used to determine which genes get to stay
  selection_crit_table$effect_size_gs_sufficient <-
    with(selection_crit_table,
      max_t >= min_t & max_betas >= min_beta |
      (max_diff_from_us >= min_max_diff_from_us & Amax >= 2) |
      diff_max_mono_reporter != 'none'
    )

  selection_crit_table$sigma_gs_sufficient <- with(selection_crit_table,
    sigma2 >= 2
  )

  selection_crit_table$sc_GDR_gs_sufficient <-
    with(selection_crit_table,
    vitro_GDR >= min_GDR & vivo_GDR >= min_GDR)

  # selection_crit_table$Amean_gs_sufficient <- with(selection_crit_table,
  #   Amean >= .05)

  selection_crit_table$Amax_gs_sufficient <- with(selection_crit_table,
    Amax >= min_Amax)

  selection_crit_table$congruent_timepoints_gs_sufficient <-
    with(selection_crit_table,
      ((message_cp == 'OK' & n_congruent_timepoints >= 3) |
      (min_max_congruence >= 3 | max_max_congruence >= 3)) |
      diff_max_mono_reporter != 'none'
    )

  ## The set of rules used to determine which genes get to stay
  selection_crit_table$included <<- with(selection_crit_table,
    effect_size_gs_sufficient &
      sigma_gs_sufficient &
      # congruent_timepoints_gs_sufficient &
      sc_GDR_gs_sufficient &
      Amax_gs_sufficient
  )
}


find_cells <- function(
    so = exp6369_sc,
    grp = 'Unstimulated in vitro',
    duration = 2
  ) {
  grp_bools <- so@meta.data$stim_group == grp
  time_bools <- so@meta.data$duration == duration
  colnames(so@assays$SCT[, grp_bools & time_bools])
}


make_scatter <- function(dtf = merged_b, x, y) {
  ggplot(dtf, aes_string(x = x, y = y)) +
    geom_hline(
      yintercept=0, linetype="dashed",
      color = "red", size=.2) +
    geom_vline(
      xintercept=0, linetype="dashed",
      color = "red", size=.2) +
    geom_point(alpha = .5, color = 'grey20') +
    geom_abline(
      intercept=0, slope = 1, linetype="dashed",
      color = "red", size=.3)
}


tally_df <- function(dtf, ...) {
  # group_vars <-
  dplyr::group_by(dtf, !!group_vars) %>%
    dplyr::summarize(N = n()) %>%
    dplyr::mutate(perc = 100 * N / sum(N))
}


load_sc_data <- function(run_magic = F, redo = F) {
  if (!exists('so_6489') || redo) {
    so_6489 <- readRDS(compile_so_fns('6489')[['filtered']])
    # colnames(so_6489@meta.data)
    # str(so_6489@meta.data$condition_name)
    so_6489@meta.data$condition_name <-
      so_6489@meta.data$condition_name %>% factor()
    levs <- naturalsort::naturalsort(levels(so_6489@meta.data$condition_name))
    # levs <- c(levs[(21:24)], levs[-(17:24)], levs[17:20])
    levs <- c(levs[17:20], levs[-(17:20)])
    so_6489@meta.data$condition_name %<>% factor(levels = levs)
    so_6489 <<- so_6489
    if (run_magic) {
      so_6489_m <<- magic(so_6489, genes=VariableFeatures(so_6489))
      so_6489_m@active.assay <<- 'MAGIC_SCT'
    }
  }

  if (!exists('so_5310') || redo) {
    so_5310 <- readRDS(compile_so_fns('5310')[['filtered']])
    so_5310@meta.data$condition_name %<>% factor()
    levs <- c(levels(so_5310@meta.data$condition_name)[7],
      levels(so_5310@meta.data$condition_name)[-7])
    so_5310@meta.data$condition_name %<>% factor(levels = levs)
    so_5310 <<- so_5310
    if (run_magic) {
      so_5310_m <<- magic(so_5310, genes=VariableFeatures(so_5310))
      so_5310_m@active.assay <<- 'MAGIC_SCT'
    }
  }

  if (!exists('so_6369') || redo) {
    so_6369 <- readRDS(compile_so_fns('6369')[['filtered']])
    so_6369@meta.data$condition_name %<>% factor()
    so_6369 <<- so_6369
    if (run_magic) {
      so_6369_m <<- magic(so_6369, genes=VariableFeatures(so_6369))
      so_6369_m@active.assay <<- 'MAGIC_SCT'
    }
  }
}


filter_gene_reproducibility_score <- function(...)
  UseMethod('filter_gene_reproducibility_score')


filter_gene_reproducibility_score.integer <- function(level = NULL) {
  if (is.null(level)) return(NULL)

  allowed_genes <-
    targets::tar_read(gene_reproducibility_score) %>%
    dplyr::filter(score_group <= .env[['level']]) %>%
    dplyr::pull(gene)

  return(allowed_genes)
}
# filter_gene_reproducibility_score(1L)
# filter_gene_reproducibility_score(2L)


filter_gene_reproducibility_score.data.frame <- function(
  dtf, cn = 'gene', level = NULL) {

  if (is.null(level)) return(dtf)

  if (!cn %in% colnames(dtf)) {
    rlang::warn(glue::glue('Could not find column: {cn}'))
    return(dtf)
  }

  allowed_genes <-
    targets::tar_read(gene_reproducibility_score) %>%
    dplyr::filter(score_group <= .env[['level']]) %>%
    dplyr::select(gene)

  dplyr::inner_join(dtf, allowed_genes)
}


load_genes <- function(genelist, so = NULL) {
  if (stringr::str_detect(genelist, '^vst') && !is.null(so)) {
    N_genes <- stringr::str_replace(genelist, 'vst', '') %>%
      as.numeric()
    ## TODO Might wanna have a look at this.
    if (length(so@assays$SCT@SCTModel.list) > 1) {
      so <- SCTransform(so, vars.to.regress = NULL)
    }
    so <- Seurat::FindVariableFeatures(
      so, nfeatures = N_genes, selection.method = 'vst')
    return(VariableFeatures(so))
  } else if (!genelist %in% c('', 'all')) {
    return(read_geneset(glue::glue('{genelist}_genes')))
  } else if (!is.null(so)) {
    return(detectable_genes(so))
  } else {
    return(NULL)
  }
}


set_var_genes <- function(
  so,
  assay = DefaultAssay(so),
  GDR_thresh = NULL,
  filter_gene_reproducibility = NULL,
  genes = NULL,
  genelist = NULL) {

  detectable_genes <- detected_genes(so, assay = assay)

  if (!is.null(genes) && length(genes) > 0) {
    detectable_genes <- intersect(detectable_genes, genes)
  }

  if (!is.null(genelist)) {
    detectable_genes <- intersect(detectable_genes,
      load_genes(genelist = genelist, so = so)
    )
  }

  if (!is.null(filter_gene_reproducibility)) {
    rep_genes <-
      filter_gene_reproducibility_score(filter_gene_reproducibility)
    detectable_genes <- intersect(detectable_genes, rep_genes)
  }

  ## Apply GDR filtering; require gene to be detected in at least x
  ## % of cells in one or more conditions
  if (!is.null(GDR_thresh) && GDR_thresh > 0) {
    GDR_passing_genes <- get_GDR_table(
      extract_sc(so), thresh = GDR_thresh,
      group_var = 'condition_name'
    )
    detectable_genes <- intersect(detectable_genes, GDR_passing_genes)
  }

  VariableFeatures(so) <- detectable_genes
  return(so)
}


identify_non_reducible_genes <- function(
  # exp_combs = c('5029', '6369', '6434', '5310'),
  # sc_mode = 'pseudobulk_high_fidelity',
  # genelist = 'informativeV15'
  genelist = NULL) {

  source('~/MirjamHoekstra/R/init.R')
  so_P <- 
    list(
      '5029' = tar_read(bulk_5029_so),
      '6434' = tar_read(bulk_6434_so),
      '5310' = tar_read(sc_pseudobulk_5310),
      '6369' = tar_read(sc_pseudobulk_6369)
    ) %>%
    normalize_sample_list(
      GDR_thresh = NULL,
      split_5310 = F,
      merge_experiments = T,
      genelist = NULL
    )

  M_orig <- as.matrix(GetAssayData(so_P, 'counts', assay = 'SCT'))
  stopifnot(!any(is.na(M_orig)))
  stopifnot(all(M_orig >= 0))
  library(preprocessCore)
  M_corr <- preprocessCore::normalize.quantiles(M_orig)
  dimnames(M_corr) <- dimnames(M_orig)
  M_corr <- subset_feats(M_corr, read_geneset(genelist))

  id_vars = c('tnf_conc', 'ifn_conc', 'duration')
  meta <- fill_in_tnf(so_P@meta.data) %>%
    dplyr::mutate(sample_idx = 1:n())
  dup_settings <-
    meta[which(duplicated(meta[, id_vars])), id_vars] %>%
    set_rownames(NULL) %>%
    na.omit() %>%
    dplyr::distinct() %>%
    # head(n = 5) %>%
    { . }

  stats_comb <-
    1:nrow(dup_settings) %>%
    furrr::future_map_dfr(function(i) {
    # purrr::map_dfr(function(i) {
      l_meta <- dplyr::right_join(meta, dup_settings[i, ],
        by = colnames(dup_settings))
      id_cols <- l_meta[, c('exp', 'condition_name')]
      stopifnot(all(!duplicated(id_cols)))
      # source('~/MirjamHoekstra/R/init.R')
      stats <- compute_gene_stats(M_corr[, l_meta$sample_idx])
      stats$i <- i
      stats$N_pseudo_reps <- nrow(l_meta)
      stats <- bind_cols(stats, dup_settings[i, ])
      # stats$q9q1om_b <- cut(stats$q9q1om, c(0, .5, 1, 1.5, 2, 4))
      return(stats)
    })

  dup_settings$i <- 1:nrow(dup_settings)
  dup_settings$condition_name <-
    purrr::map_chr(1:nrow(dup_settings), function(i) {
      with(dup_settings[i, ],
        stringr::str_replace(glue::glue('\\
          {make_flag(duration)}\\
          {make_flag(tnf_conc)}\\
          {make_flag(ifn_conc)}'
        ), '^-', '')
      )
    })
  stats_comb <- dplyr::left_join(stats_comb, dup_settings)

  # save_target(stats_comb, 'non_reproducible_genes')
  # tst <- tar_read(non_reproducible_genes)
  return(stats_comb)
}
# res <- identify_non_reducible_genes()
# save_target(res, 'non_reproducible_genes')


#' Cluster genes and define a set of gene modules, genes acting
#' similarly across the input (bulk) samples
#'
#'
cluster_genes <- function(genes, k, 
  so = NULL, 
  dist_f = 'pearson', 
  clust_method = 'complete', 
  trans = function(x) log2(x + 1),
  # trans = identify,
  min_var = .5) {

  if (is.null(so)) {
    so <- tar_read(bulk_5029_so)
    so <- so[, so@meta.data$sn_dilution == '0']
  }

  M <- 
    so2M(so, datatype = 'counts') %>%
    subset_feats(genes) %>%
    trans() %>%
    { .[which(apply(., 1, var) > min_var), ] } %>%
    { t(scale(t(.))) } %>%
    na.omit() %>%
    { . }

  gco <- gen_clust_object(t(M), 
    dist_f = dist_f,
    clust_method = clust_method
  )
  c_a <- gen_tree_split(gco, cluster_k = k) 
  gs <- 
    purrr::map(unique(c_a), function(cl) {
      names(which(c_a == cl))
    }) %>%
    setNames(paste0('GM', 1:k))

  return(list(clust_object = gco, c_a = c_a, genesets = gs))
}
