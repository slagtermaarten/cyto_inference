comp_levs_lu = list(
  '5310' = list(
    tibble(
      'stim_group' = c('Unstimulated in vitro', '10 ng/ml TNFa')
    ),
    tibble(
      'stim_group' = c('Unstimulated in vitro', '100 ng/ml IFNy 10 ng/ml TNFa')
    ),
    tibble(
      'stim_group' = c('Unstimulated in vitro', '100 ng/ml IFNy')
    ),
    tibble(
      'stim_group' = c('100 ng/ml IFNy', '100 ng/ml IFNy 10 ng/ml TNFa')
    ),
    tibble(
      'stim_group' = c('10 ng/ml TNFa', '100 ng/ml IFNy 10 ng/ml TNFa')
    ),
    tibble(
      'stim_group' = c('Exposed to T-cells in vivo', 'Unexposed in vivo')
    )
  ),
  '6369' = list(
    tibble('ifn_conc' = c(100, 0), 'duration' = c(24, 24)),
    tibble('ifn_conc' = c(100, 0), 'duration' = c(6, 6)),
    tibble('ifn_conc' = c(100, 100), 'duration' = c(2, 24))
  ),
  '6489' = list(
    tibble('sn_dilution' = c(1, 1/625), 'duration' = c(24, 24)),
    tibble('sn_dilution' = c(1, 0), 'duration' = c(24, 24)),
    tibble('sn_dilution' = c(1, 1), 'duration' = c(2, 24)),
    tibble('sn_dilution' = c(1, 1), 'duration' = c(2, 24)),
    tibble('sn_dilution' = c(1, 1), 'duration' = c(6, 24)),
    tibble('sn_dilution' = c(1, 1), 'duration' = c(12, 24))
  ),
  '6493' = list(
    tibble(
      'stim_group' = c('Exposed to T-cells in vivo', 
        'Unexposed in vivo'),
      'duration' = c(44, 44)
    ),
    tibble(
      'stim_group' = c('Exposed to T-cells in vivo', 
        'Unexposed in vivo'),
      'duration' = c(16, 16)
    ),
    tibble(
      'stim_group' = c('Exposed to T-cells in vivo', 
        'Exposed to T-cells in vivo'),
      'duration' = c(44, 16)
    )
  ),
  '6600' = list(
    tibble(
      'stim_group' = c('Unstimulated in vitro', '10 ng/ml TNFa')
    ),
    tibble(
      'stim_group' = c('Unstimulated in vitro', '100 ng/ml IFNy 10 ng/ml TNFa')
    ),
    tibble(
      'stim_group' = c('10 ng/ml TNFa', '100 ng/ml IFNy 10 ng/ml TNFa')
    )
  ),
  '6601' = list(
    tibble(
      'stim_group' = c('Unexposed in vivo', '10 ng/ml TNFa'),
      'duration' = c(6, 6)
    ),
    tibble(
      'stim_group' = c('Unexposed in vivo', '100 ng/ml IFNy 10 ng/ml TNFa'),
      'duration' = c(6, 6)
    ),
    tibble(
      'stim_group' = c('10 ng/ml TNFa', '100 ng/ml IFNy 10 ng/ml TNFa'),
      'duration' = c(6, 6)
    ),
    tibble(
      'stim_group' = c('Unexposed in vivo', '100 ng/ml IFNy'),
      'duration' = c(24, 24)
    ),
    tibble(
      'stim_group' = c('Unexposed in vivo', '100 ng/ml IFNy'),
      'duration' = c(48, 48)
    ),
    tibble(
      'stim_group' = c('100 ng/ml IFNy', '100 ng/ml IFNy'),
      'duration' = c(24, 48)
    )
  )
)


extract_varying_comp <- function(comp_levs) {
  cn_v <-
    names(which(map_lgl(comp_levs, ~length(unique(.x)) > 1)))
  if (is.null(cn_v)) {
    rlang::abort('NULL cn_v')
  } else if (length(cn_v) != 1) {
    rlang::warn('Multiple varying columns in comp_levs')
    return(NULL)
  } else {
    return(cn_v)
  }
}


determine_comp_name <- function(
  comp_levs = tibble(
    'sn_dilution' = c(1, 1/625), 
    'duration' = c(2, 24)
  ), 
  ...) {
  cn_v <- extract_varying_comp(comp_levs)
  if (is.null(cn_v)) return(NULL)
  cn_nv <- setdiff(colnames(comp_levs), cn_v)

  id_part <- if (cn_v != 'condition_name') {
    glue::glue('{print_names_nc[cn_v]}: ') 
  } else {
    ''
  }
  comp_part <- out <- glue::glue('{id_part}\\
    {paste(comp_levs[[cn_v]], collapse=\' vs. \')}')

  if (length(cn_nv) > 0) {
    out <- glue::glue('{comp_part} at {print_names_nc[cn_nv]} {comp_levs[[cn_nv]][1]}')
  } else {
    out <- comp_part
  }

  return(out)
}
# determine_comp_name()


compare_gs_between_levs <- function(gs, comp_levs, dtf, 
  return_type = 'auc', ...) {
  cn <- colnames(comp_levs)

  if (any(!c(cn, gs) %in% colnames(dtf))) {
    return(NULL)
  }

  cn_v <- extract_varying_comp(comp_levs)

  ## Ensure comp_levs and gs are of same type
  char_cn <- names(map_chr(comp_levs, class) == 'character')
  fac_cn <- names(which(map_chr(dtf[, char_cn, drop = F], class) == 'factor'))
  if (length(fac_cn) > 0) {
    for (cn_l in fac_cn) {
      comp_levs[[cn_l]] <- factor(comp_levs[[cn_l]], 
        levels = levels(dtf[[cn_l]]))
    }
  }
  
  t_dat <-
    dtf %>%
    # dplyr::mutate(across(any_of(char_cn), factor_to_numeric)) %>%
    dplyr::right_join(comp_levs, by = cn) %>%
    dplyr::mutate(across(any_of(cn), factor)) %>%
    { . }

  t_dat[['y']] <- t_dat[[cn_v]]
  t_dat[['x']] <- t_dat[[gs]]
  t_dat <- t_dat[, c('x', 'y'), drop = F]

  if (F) {
     
    p <- ggplot(t_dat, aes(y, x)) + 
      geom_boxplot() +
      geom_point()
    print_plot_eval(print(p),
      width = 17.4, height = 10,
      filename = file.path(Sys.getenv('img_dir'),
        glue::glue('test.pdf')))

  }

  if (return_type == 'auc') {
    # out <- pROC::roc(t_dat, y, x)$auc
    out <- tryCatch(
      yardstick::roc_auc(t_dat, y, x)$.estimate, 
      error = function(e) { NA_real_ })  %>%
      { 
        if (!is.na(.) && . < .5) 1 - . else .
      }
    return(out)
  } else if (return_type == 'plot') {
    co <- yardstick::roc_curve(t_dat, y, x)
    lab <- round(yardstick::roc_auc(t_dat, y, x)$.estimate, 2) %>%
      { if (. < .5) 1 - . else . } %>%
      # { glue::glue('PR AUC: {.}') }
      { glue::glue('AUROC: {.}') }
    p <- autoplot(co) + annotate_npc(lab, .9, .9, hjust = 1)
    return(p)
  }
}


filter_so_score <- function(so_score, cn_mode) {
  group_var <- group_var_cn_mode[[cn_mode]]

  if ('Ag' %in% colnames(so_score@meta.data)) {
    so_score <- add_Ag_condition_name(so_score, new_cn = group_var)
  }

  if (cn_mode == 'all') {
  } else if (cn_mode == 'simple') {
  } else if (cn_mode == 'vivo') {
    so_score <- so_score[, so_score@meta.data[['sample_origin']] == 'in_vivo']
  } else if (cn_mode == 'vitro_simple') {
    so_score <- so_score[, so_score@meta.data[['sample_origin']] == 'in_vitro']
  } else if (cn_mode == 'vitro_extended') {
    so_score <- so_score[, so_score@meta.data[['sample_origin']] == 'in_vitro']
  # } else if (cn_mode == 'sc_digestion') {
  #   so_score <- so_score[, so_score@meta.data[['sc_digestion']] == 'Ag-']
  } else if (cn_mode == 'Ag-') {
    so_score <- so_score[, so_score@meta.data[['Ag']] == 'Ag-']
  }
  return(so_score)
}


group_var_cn_mode <- list(
  'all' = 'condition_name',
  'simple' = 'cn_simple',
  'vivo' = 'condition_name',
  'vitro_simple' = 'cn_simple',
  'vitro_extended' = 'condition_name',
  'sc_digestion' = 'condition_name',
  'Ag' = 'condition_name'
)


cns_cn_mode <- list(
  'all' = c('stim_group', 'duration', 'frozen', 'sc_digestion', 'Ag'),
  'simple' = c('stim_group', 'duration', 'Ag'),
  'vivo' = c('stim_group', 'duration'),
  'vitro_simple' = c('stim_group', 'duration'),
  'vitro_extended' = c('stim_group', 'duration', 'frozen',
    'sc_digestion'),
  'sc_digestion' = c('stim_group', 'duration', 'sc_digestion'),
  'Ag' = c('stim_group', 'Ag', 'duration')
)


compare_conditions_by_gs <- function(
  dtf, gs, comp_levs = comp_levs_lu[[experiment]]) {

  out <-
    tidyr::expand_grid(
      gs = gs,
      comp_levs = comp_levs
    ) %>%
    dplyr::mutate(comp_name =
      map(comp_levs, determine_comp_name)) %>%
    dplyr::filter(!sapply(comp_name, is.null)) %>%
    dplyr::mutate(AUROC =
      pmap_dbl(., compare_gs_between_levs, dtf = dtf)) %>%
    { 
      postfix <- extract_common_postfix(.$gs)
      .$gs <- stringr::str_replace(.$gs, postfix, '')
      .
    }
  return(out)
}

