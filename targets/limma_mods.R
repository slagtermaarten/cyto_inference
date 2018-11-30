limma_mod_table <- tribble(
    ~experiment, ~lf,
    '5029', rlang::sym('limma_model_5029'),
    # 6434, rlang::sym(limma_model_6434),
    # '5310', rlang::sym('limma_model_5310'),
    '6369', rlang::sym('limma_model_6369'),
    '6489', rlang::sym('limma_model_6489'),
    '6493', rlang::sym('limma_model_6493'),
    '6600', rlang::sym('limma_model_6600'),
    '6601', rlang::sym('limma_model_6601')
  ) %>%
  dplyr::mutate(i = 1:n()) %>%
  dplyr::mutate(terms_obj = 
      rlang::syms(glue::glue('limma_mod_cyto_{experiment}'))) %>%
  dplyr::mutate(pb_obj =
    case_when(
      experiment %in% sc_e ~
      rlang::syms(glue::glue('sc_pseudobulk_{experiment}')),
      experiment %in% bulk_e ~
      rlang::syms(glue::glue('bulk_{experiment}_so'))
    )
  ) %>%
  { . }

limma_model_targets <- tarchetypes::tar_map(
  names = experiment,
  values = limma_mod_table,
  unlist = T,

  tar_target(limma_mod_cyto, {
    limma_mod <- lf(so = pb_obj)
    extract_terms(limma_mod, max_p = 1, min_estimate = .00)
  }, iteration = 'list')
)
