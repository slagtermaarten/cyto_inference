## Group meta vars into classes, such that meta vars that always go
## together can be more easily invoked
param_types <- list(
  'reference' = c('reference'),
  'expression' = c('query', 'genelist', 'include_variable_feats'),
  'scaling' = c('GDR', 'Z_scale'),
  'neighbourhood' = c('k', 'd'),
  'min_neighbourhood_size' = c('min_neighbourhood_size'),
  'DA_type' = c('DA_type'),
  'transact' =
    c('kernel', 'N_source_components', 'N_target_components', 'N_PV'),
  # , 'report' = c('min_early_var', 'D_gamma')
  'kpca' = c('kernel', 'N_NLPC'),
  # , 'report' = c('min_early_var', 'D_gamma')
  'report' = c()
)


order_types <- function(types)
  ## Sort the types in the most logical order, let object param_types
  ## define that order, also (silently) remove requested types that
  ## aren't recognized.  The first argument of intersect() determines
  ## order of result in R
  intersect(names(param_types), types)


#' To a parameter table ('source_dtf'), add object references that are characterized by
#' the meta var types in 'types'. Column names will be added to the
#' parameter table for each variable in \code{source_name}, appended
#' by an ID string of all meta variables
#'
#'
add_obj_ref <- function(
  source_dtf, source_name,
  types = c('expression', 'neighbourhood'), id = NULL) {

  types <- order_types(types)
  lid <- suppressWarnings(
    rlang::flatten_chr(param_types[types])) %>%
    intersect(colnames(source_dtf))
  postfix <- apply(source_dtf[, lid], 1, paste0, collapse = '_')
  postfix <- as.character(glue::glue('_{postfix}'))
  postfix <- stringr::str_replace_all(postfix, ' ', '')
  postfix <- stringr::str_replace_all(postfix, '_0\\.00_', '_0_')
  postfix <- stringr::str_replace_all(postfix, '_0\\.0_', '_0_')
  if (F && any(stringr::str_detect(source_name, 'titration'))) {
    stringr::str_subset(source_name, 'titration')
    browser()
  }

  if (any(stringr::str_detect(postfix, ' '))) {
    rlang::warn('Spaces detected in postfix')
    if (interactive()) {
      browser()
    }
  }

  for (vn in source_name) {
    source_dtf[[glue::glue('{vn}_obj')]] <-
      rlang::syms(glue::glue('{vn}{postfix}'))
  }

  if (!is.null(id)) {
    source_dtf[[glue::glue('{id}_id')]] <- postfix
  }

  return(source_dtf)
}


types_to_names <- function(
  types = c('expression', 'neighbourhood')) {
  rlang::flatten_chr(unname(param_types[order_types(types)]))
}


subset_param_grid <- function(
  source_dtf,
  id_col_types = c('expression', 'neighbourhood',
    'min_neighbourhood_size')) {
  undesired <-
    types_to_names(setdiff(names(param_types), id_col_types)) %>%
    ## If overlap between id_col_types, ensure that desired columns
    ## are included in both desired and undesired types, remain
    ## included
    setdiff(types_to_names(id_col_types))
  ## Negatively select undesired columns such that unrecognized
  ## columns will remain unaltered
  source_dtf %>%
    dplyr::select(-any_of(undesired)) %>%
    dplyr::distinct(across(any_of(types_to_names(id_col_types))),
      .keep_all = T) %>%
    { . }
}
# subset_param_grid(PB_grid)
# subset_param_grid(PB_grid, c('reference', 'expression'))

add_static_obj_names <- function(dtf) {
  dtf %>%
    dplyr::mutate(ref_so =
      ifelse(stringr::str_detect(reference, '^\\d{4}$'),
        rlang::syms(glue::glue('bulk_{reference}_so')),
        rlang::syms('bulk_comb_so')
      )
    ) %>%
    dplyr::mutate(q_sc_so =
      rlang::syms(glue::glue('filtered_cleaned_so_{query}'))) %>%
    dplyr::mutate(q_pb_so =
      rlang::syms(glue::glue('sc_pseudobulk_{query}'))) %>%
    # dplyr::mutate(PB_report_id = glue::glue('
    #     {reference}_{query}_{genelist}\\
    #     _{N_source_components}_{N_target_components}\\
    #     _{D_gamma}_{100*round(min_early_var, 2)}')) %>%
    # dplyr::mutate(PB_report_name = glue::glue('
    #     PB_TRANSACT_report_{PB_report_id}')) %>%
    # dplyr::mutate(param_idx = 1:n()) %>%
    { . }
}


rm_var_parts <- function(
  v, 
  grid_dtfs = list(targets_env$NH_t_grid, targets_env$NH_e_grid), 
  retain_id = F) {

  ## Compile a list of all variable parts
  strings <-
    purrr::map(grid_dtfs, function(grid_dtf) {
      id_vars <- stringr::str_subset(colnames(grid_dtf), '_id$')
      purrr::map(id_vars, function(id_var) {
        unique(grid_dtf[[id_var]])
      })
    }) %>%
    unlist(recursive = T) %>%
    stringr::str_replace('_$', '') %>%
    stringr::str_replace('^_', '') %>%
    unique() %>%
    { .[base::order(nchar(.), ., decreasing = T)] } %>%
    { . }

  if (retain_id) {
    out <- 
      purrr::map(v, function(.x) {
        hits <- which(stringr::str_detect(.x, strings))
        if (length(hits) == 0) {
          return(.x)
        } else {
          return(strings[hits[1]])
        }
      }) %>%
      unlist()
  } else {
    out <- v
    for (cs in strings) {
      out <- stringr::str_replace(out, cs, '')
    }
  }

  out <- stringr::str_replace(out, '_$', '')

  return(out)
}


iso_var_parts <- function(
  v, 
  grid_dtfs = list(
    NH_t_grid, NH_e_grid,
    targets_env$NH_t_grid, targets_env$NH_e_grid)) {

  stopifnot(!all(map_lgl(grid_dtfs, is.null)))
  rm_var_parts(v, grid_dtfs = grid_dtfs, retain_id = T)
}


find_var_part <- function(
  v = 'comb_6601_informativeV15_0_TRUE_30_3_10_15_5_3', 
  grid_dtfs = list(NH_grid), retain_id = F) {
   
  strings <-
    seq_along(grid_dtfs) %>%
    purrr::map_dfr(function(i) {
      grid_dtf <- grid_dtfs[[i]]
      id_vars <- stringr::str_subset(colnames(grid_dtf), '_id$')
      purrr::map_dfr(maartenutils::auto_name(id_vars), function(id_var) {
        unique(grid_dtf[[id_var]]) %>%
        stringr::str_replace('_$', '') %>%
        stringr::str_replace('^_', '') %>%
        { tibble(var = id_var, string = .) } %>%
        { . }
      }) %>%
      dplyr::mutate(i = i) %>%
      # tibble::deframe() %>%
      { . }
    }) %>%
    # unlist(recursive = T) %>%
    # unique() %>%
    # { .[base::order(nchar(.), ., decreasing = T)] } %>%
    { . }

  return(strings[which(v == strings$string), ])
}


subset_var_part <- function(var_part, grid_dtfs = list(NH_grid)) {
  bm <- find_var_part(var_part, grid_dtfs = grid_dtfs)
  grid_dtf <- grid_dtfs[[bm$i]]
  grid_dtf[which(grid_dtf[[bm$var]] == glue::glue('_{var_part}')), ] %>%
    # dplyr::select(matches('id')) %>%
    { . }
}


subset_grid_dtf <- function(
  name = paste0("NH_gamma_titration_embedding_comb_6369", 
    "_informativeV15_FALSE_0_FALSE_10_100_0_10_100_10"),
  grid_dtf = NH_t_grid, 
  remove_obj_refs = T,
  maintain_query = F) {

  ## Find column name to look in
  cn <- recover_garbled_string(name, colnames(grid_dtf)) %>%
    unname() %>%
    pluck(1)
  out <- grid_dtf[which(as.character(grid_dtf[[cn]]) == name), ]
  out[[cn]] <- as.character(out[[cn]])
  if (remove_obj_refs) {
    remove_idxs <- 
      which(
        map(out, class) == 'list' |
        stringr::str_detect(colnames(out), '_id') |
        stringr::str_detect(colnames(out), '_name')
      )
    if (maintain_query) {
      remove_idxs <- remove_idxs %>%
        { .[names(.) != cn] }
    }
    out <- out[-remove_idxs]
  }
  return(out)
}


gen_NH_grid <- function(
  extra_settings = NULL,
  check_target_name = FALSE) {

  base_dtf <-
    tidyr::expand_grid(
      tidyr::expand_grid(
        # reference = c('5029', '6434', 'comb'),
        reference = c('comb'),
        N_source_components = c(15L, 10L, 15L),
        query = c('6369', '6489', '6493', '6600', '6601')
      ) %>%
      bind_rows(tibble(
        reference = c('5029'),
        N_source_components = c(15L),
        query = c('6489')
      )),
      tibble(
        N_target_components = c(20L, 15L, 10L, 5L),
        N_PV = N_target_components,
        d = N_target_components
      ),
      tibble(
        k = c(10L, 30L, 30L, 30L),
        min_neighbourhood_size = c(0L, 0L, 10L, 30L)
      ),
      include_variable_feats = FALSE,
      genelist = c(
        'informativeV15',
        # 'informativeV15_and_sc_var',
        # 'informativeV15_monotonic',
        # 'informativeV15_monoreporter',
        NULL
      ),
      # tibble(
      #   GDR = c(0, .1, .1),
      #   Z_scale = c(FALSE, TRUE, FALSE),
      # ),
      tibble(
        GDR = c(0, 0, .25),
        Z_scale = c(FALSE, TRUE, TRUE)
      ),
      # k = c(30L, 100L),
      # min_neighbourhood_size = c(0L)
      # D_gamma = 10,
      # min_early_var = c(0, .75, .9)
    ) %>%
    { . }

  if (!is.null(extra_settings)) {
    for (cn in
      intersect(colnames(base_dtf), colnames(extra_settings))) {
      base_dtf[[cn]] <- NULL
    }
    base_dtf <- dplyr::distinct(base_dtf)
    base_dtf <- base_dtf %>% tidyr::expand_grid(extra_settings)
  }

  NH_grid <-
    base_dtf %>%
    dplyr::distinct(.keep_all = T) %>%
  # _{round(min_early_var, 2)}\\
  # _{D_gamma}
  # dplyr::mutate(NH_report_id = glue::glue('
  #     {reference}_{query}_{genelist}\\
  #     _{round(GDR, 2)}\\
  #     _{Z_scale}\\
  #     _{k}_{d}\\
  #     _{min_neighbourhood_size}\\
  #     _{N_source_components}_{N_target_components}_{N_PV}')) %>%
    add_static_obj_names() %>%
    add_obj_ref(
      source_name = c('NH_i_feats', 'NH_agg_neighbourhoods', 
        'NH_test_DA'),
      types = c('expression', 'neighbourhood'),
      id = 'milo'
    ) %>%
    add_obj_ref(
      source_name = c('NH_so'),
      types = c('expression', 'neighbourhood',
        'min_neighbourhood_size'),
      id = 'query'
    ) %>%
    add_obj_ref(
      source_name = c(
        'NH_harmonized_Ms',
        'NH_pca_screes',
        'NH_expression_clustering',
        'Nhood_stats'),
      types = c('reference', 'expression', 'scaling', 'neighbourhood',
        'min_neighbourhood_size'),
      id = 'int_preproc'
    ) %>%
    add_obj_ref(
      source_name = c(
        # 'NH_gamma_titration_cosine',
        # 'NH_optimal_log_gamma',
        'NH_gamma_titration_embedding',
        # 'NH_score_error',
        'NH_Nhood_error',
        'NH_Nhood_gaussian_error',
        'NH_TRANSACT_reference_reconstruction',
        'ann_NH_TRANSACT_reference_reconstruction',
        'NH_dtf'
      ),
      types = c('reference', 'expression', 'scaling', 'neighbourhood',
        'min_neighbourhood_size', 'transact', 'report'),
      id = 'transact'
    ) %>%
    dplyr::mutate(NH_report_name =
      glue::glue('NH_TRANSACT_report{transact_id}')) %>%
    # dplyr::filter(!(query %in% c('6600', '6601') & N_target_components > 5)) %>%
    { . }

  ## Since target instantation is burdensome, minimize it if possible
  if (check_target_name &&
      exists('target_name') &&
      stringr::str_detect(target_name, 'NH_TRANSACT')) {
    if (!target_name %in% NH_grid$NH_report_name) {
      print(recover_garbled_string(target_name,
          NH_grid$NH_report_name))
      rlang::abort('target_name not found')
    }
    NH_grid <- NH_grid %>%
      dplyr::filter(NH_report_name == target_name)
  }

  detected_string <-
    NH_grid$NH_optimal_log_gamma_obj %>%
    purrr::map_chr(~as.character(.x)) %>%
    stringr::str_detect(' ') %>%
    any()
  if (detected_string) {
    rlang::abort('Spaces detected')
  }

  return(NH_grid)
}
