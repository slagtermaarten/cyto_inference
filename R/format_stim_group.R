order_stim_group <- function(...) UseMethod('order_stim_group')

order_stim_group.character <-
order_stim_group.factor <- function(v) {
  v <- stringr::str_replace(v, '^\\s', '')
  v <- stringr::str_replace(v, '\\s$', '')

  # source('~/MirjamHoekstra/R/init.R')
  bsa <- 
    harmonize_bind_rows(
      tar_read(sample_annotation_exp5029), 
      tar_read(sample_annotation_exp6434), 
      tar_read(sc_6489_sample_annotation)
    ) %>%
    # data.table::rbindlist(
    #   list(
    #     tar_read(sample_annotation_exp5029), 
    #     tar_read(sample_annotation_exp6434), 
    #     tar_read(sc_6489_sample_annotation)), 
    #   fill = TRUE) %>%
    dplyr::mutate(duration = as.numeric(duration)) %>%
    dplyr::mutate(duration = factor_to_numeric(duration)) %>%
    dplyr::mutate(sn_dilution = as.numeric(sn_dilution)) %>%
    dplyr::mutate(tnf_conc = as.numeric(tnf_conc)) %>%
    dplyr::mutate(ifn_conc = as.numeric(ifn_conc)) %>%
    dplyr::mutate(sn_dilution = ifelse(sn_dilution > 1, 1/sn_dilution,
        sn_dilution)) %>%
    dplyr::mutate(sn_dilution = ifelse(is.na(sn_dilution), 0,
        sn_dilution)) %>%
    dplyr::mutate(ifn_conc = ifelse(is.na(ifn_conc), 0,
        ifn_conc)) %>%
    dplyr::mutate(tnf_conc = ifelse(is.na(tnf_conc), 0,
        tnf_conc)) %>%
    dplyr::arrange(ifn_conc, tnf_conc, sn_dilution, duration) %>%
    # order_condition_name() %>%
    # dplyr::arrange(condition_name)
    { . }
  # bsa[, c('stim_group', 'duration', 'ifn_conc', 'tnf_conc', 'sn_dilution')]
  # bsa$sn_dilution

  mouse_levs <- c('Mix PBS', 'Ag-GAS + T', 'Mix + T')
  b_ord <- c(as.character(unique(bsa$stim_group)), mouse_levs)

  front_matter <- c('Unstimulated in vitro', 'Unstimulated')
  if (!'Unstimulated in vitro' %in% v &&
      'Unexposed in vivo' %in% v) {
    front_matter <- c(front_matter, 'Unexposed in vivo')
  }
  front_matter <- intersect(front_matter, v)
  recog <- intersect(b_ord, v)
  recog <- c(front_matter, setdiff(recog, front_matter)) 

  back_matter <- 
    unique(stringr::str_subset(as.character(v), 'vivo')) %>%
    setdiff(front_matter)
  if (length(back_matter) > 0) {
    T_cell <- 
      unique(stringr::str_subset(as.character(v), 'T-cell')) %>%
      naturalsort::naturalsort()
    non_T_cell <- 
      naturalsort::naturalsort(setdiff(back_matter, T_cell))
    back_matter <- c(non_T_cell, T_cell)
    recog <- setdiff(recog, back_matter)
  }

  unrecog <- setdiff(v, recog) %>% setdiff(back_matter)

  levs <- c(recog, unrecog, back_matter)


  return(factor(v, levels = levs))
}


order_stim_group.data.frame <-
order_stim_group.tbl_df <- function(dtf, cn='stim_group') {
  if (!cn %in% colnames(dtf)) {
    return(dtf)
  } else {
    dtf[[cn]] <- order_stim_group(dtf[[cn]])
    return(dtf)
  }
}


order_stim_group.Seurat <- function(so) {
  so@meta.data <- order_stim_group(so@meta.data)
  return(so)
}


recover_stim_group <- function(...) UseMethod('recover_stim_group')


recover_stim_group.character <- function(v) {
  stringr::str_replace(v, ' -.*$', '')
}


recover_stim_group.factor <- function(v) {
  recover_stim_group(as.character(v))
}


recover_stim_group.tbl_df <- 
recover_stim_group.data.frame <- function(dtf) {
  if ('stim_group' %in% colnames(dtf) && 
    !all(is.na(dtf$stim_group))) {
    return(dtf)
  }

  if (!'condition_name' %in% colnames(dtf)) {
    rlang::warn('Neither stim_group nor condition_name found')
    return(dtf)
  }

  dtf$stim_group <- recover_stim_group(dtf$condition_name)

  return(dtf)
}


order_condition_name <- function(...)
  UseMethod('order_condition_name')


order_condition_name.data.frame <-
order_condition_name.tbl_df <- function(dtf) {
  if (!any(c('stim_group', 'condition_name') %in% colnames(dtf))) {
    return(dtf)
  }
  dtf <- recover_stim_group(dtf)
  # dtf$stim_group
  dtf <- order_stim_group(dtf)
  condition_name_levels <- dtf %>%
    dplyr::arrange(stim_group, duration) %>%
    dplyr::select(stim_group, duration) %>%
    unique() %>%
    add_condition_name() %>%
    pull(condition_name)
  dtf$condition_name <- 
    order_factor(
      dtf$condition_name, 
      exp_levels = condition_name_levels
    )
  # dtf$condition_name
  return(dtf)
}


order_condition_name.Seurat <- function(so) {
  so@meta.data <- order_condition_name(so@meta.data)
  return(so)
}


order_duration <- function(...)
  UseMethod('order_duration')


order_duration.factor <- 
order_duration.integer <- 
order_duration.character <- 
order_duration.numeric <- function(v) {
  v <- as.character(v)
  v[is.na(v)] <- 'Unknown'
  levs <- naturalsort::naturalsort(unique(v))
  factor(v, levels = levs)
}


order_duration.data.frame <-
order_duration.tbl_df <- function(dtf) {
  dur_vars <- stringr::str_subset(colnames(dtf), 'duration')
  if (length(dur_vars) == 0L) return(dtf)
  dplyr::mutate(dtf, across(any_of(dur_vars), order_duration))
}


order_duration.Seurat <- function(so) {
  so@meta.data <- order_duration(so@meta.data)
  return(so)
}


order_concentration <- function(...) UseMethod('order_concentration')


order_concentration.factor <- 
order_concentration.integer <- 
order_concentration.character <- 
order_concentration.numeric <- function(v) {
  v <- as.numeric(as.character(v))
  v[is.na(v)] <- 0
  v
}


order_concentration.data.frame <-
order_concentration.tbl_df <- function(dtf) {
  dplyr::mutate(dtf, across(matches('_conc'), order_concentration))
}



add_condition_name <- function(...) UseMethod('add_condition_name')


add_condition_name.data.frame <-
add_condition_name.tbl_df <- function(dtf) {
  duration_part <- with(dtf,
    ifelse(is.na(duration) | duration == 'Unknown', '', 
      glue::glue(' - {as.character(duration)}h')))
  dtf$condition_name <-
    as.character(glue::glue('{dtf$stim_group}{duration_part}'))
  return(dtf)
}


add_condition_name.Seurat <- function(so) {
  so@meta.data <- add_condition_name(so@meta.data)
  return(so)
}


#' Add the condition_name variable to a SC experiment metatable table
#'
#'
add_sc_condition_name <- function(sa) {
  sa <- sa %>% dplyr::select(-matches('condition_name'))
  if (!'frozen' %in% colnames(sa)) sa$frozen <- NA

  sa$condition_name <-
    sa %>%
    dplyr::mutate(duration = ifelse(is.na(duration), '',
        glue::glue(' - {as.character(duration)}h'))) %>%
    dplyr::mutate(sc_digestion = ifelse(is.na(sc_digestion), '',
        ifelse(sc_digestion, ' - SC digest', ''))) %>%
    dplyr::mutate(frozen = ifelse(is.na(frozen), '',
        ifelse(frozen, ' - frozen', ''))) %>%
    dplyr::mutate(condition_name =
      glue::glue('{stim_group}{duration}{sc_digestion}{frozen}')) %>%
    dplyr::pull(condition_name)

  if (all(grepl('SC digest', sa$condition_name)))
    sa$condition_name <-
      gsub(' - SC digest$', '', sa$condition_name)

  if (F && any(duplicated(sa$condition_name))) {
    for (lev in names(table(sa$condition_name))) {
      temp_v <- sa$condition_name[sa$condition_name == lev] 
      sa$condition_name[sa$condition_name == lev] <-
        paste0(temp_v, ' - mouse ', seq_along(temp_v))
    }
  }

  return(sa)
}


format_stim_group <- function(dtf) {
  if (!'stim_group' %in% colnames(dtf)) return(dtf)
  old_class <- class(dtf)
  old_dimnames <- dimnames(dtf)
  setDT(dtf)

  dtf <- dtf %>%
    dplyr::mutate(tnf_rank = factor_to_rank(tnf_conc)) %>%
    dplyr::mutate(ifn_rank = factor_to_rank(ifn_conc)) %>%
    dplyr::mutate(sn_rank = factor_to_rank(sn_dilution)) %>%
    { . }

  ## Format stim group for all samples
  dtf$stim_group <- dtf$stim_group %>%
    { gsub('ifng', 'IFNy', .) } %>%
    { gsub('sn', 'SN', .) } %>%
    { gsub('ng_ml', 'ng/ml', .) } %>%
    { gsub('1_2', '1/2', .) } %>%
    { gsub('1_200', '1/200', .) } %>%
    { gsub('1_20000', '1/20000', .) } %>%
    maartenutils::var_to_label(reps = c('in vivo' = 'in vivo')) %>%
    { gsub('in Vitro', 'in vitro', .) } %>%
    { gsub('in Vivo', 'in vivo', .) }

  ## Correct labels for some single cells
  dtf[grepl('^\\d$', stim_group),
     stim_group := ifelse(is.na(tnf_conc) & !is.na(ifn_conc),
                           sprintf('%.0f ng/ml IFNy',
                                   factor_to_numeric(ifn_conc)),
                           stim_group)]
  dtf[grepl('^\\d$', stim_group),
      stim_group := ifelse(!is.na(tnf_conc) & !is.na(ifn_conc),
                           sprintf('%.0f ng/ml IFNy %.0f ng/ml TNFa',
                                   factor_to_numeric(tnf_conc),
                                   factor_to_numeric(ifn_conc)),
                           stim_group)]
  dtf[grepl('^\\d$', stim_group),
      stim_group := ifelse(!is.na(tnf_conc) & is.na(ifn_conc),
                           sprintf('%.0f ng/ml TNFa',
                                   factor_to_numeric(tnf_conc)),
                           stim_group)]

  dtf[stim_group == 'unstim', stim_group := 'Unstimulated in vitro']
  dtf[order(-sn_rank, -ifn_rank, -tnf_rank), stim_group]
  levs <- dtf[order(sn_rank, ifn_rank, tnf_rank),
              c(unique(stim_group), 'Unknown')]
  dtf[, stim_group := factor(stim_group, levels = levs)]
  class(dtf) <- old_class
  rownames(dtf) <- old_dimnames[[1]]
  return(dtf)
}


#' Extract all levels from an object that is not necessarily a factor
#'
#'
e_levels <- function(v, allow_unknown = T) {
  if (is.null(v) || all(is.na(v)))
    return(NULL)
  if (is.factor(v)) {
    levs <- levels(v)
    tab <- table(v)
    levs <- intersect(levs, names(tab[tab > 0]))
  } else {
    levs <- naturalsort::naturalsort(unique(v))
  }
  levs <- setdiff(levs, NA)
  if (!allow_unknown) {
    levs <- setdiff(levs, c('Unknown', 'unknown'))
  }
  return(levs)
}


check_condition_name_consistency <- function(dtf) {
  dtf %>%
    dplyr::select(exp, condition_name, stim_group, duration) %>%
    dplyr::mutate(recov_sg = recover_stim_group(condition_name)) %>%
    dplyr::filter(as.character(stim_group) != recov_sg) %>%
    dplyr::group_by(exp) %>%
    dplyr::summarize(n()) %>%
    { nrow(.) == 0 }
}


separate_duration <- function(...)
  UseMethod('separate_duration')


separate_duration.data.frame <- 
  separate_duration.tbl_df <- function(dtf) {

  if (maartenutils::null_dat(dtf)) return(dtf)

  dtf <- numerify_regressors(dtf)

  for (conc_v in c('tnf_conc', 'ifn_conc', 'sn_dilution')) {
    cyto <- stringr::str_replace(conc_v, '_(conc|dilution)', '')
    if (!conc_v %in% colnames(dtf)) {
      next()
    }
    duration_v <- glue::glue('{cyto}_duration')
    if (duration_v %in% colnames(dtf)) {
      next()
    }
    dtf[[duration_v]] <- dtf$duration
    dtf[[duration_v]][which(is.na(dtf[[conc_v]]))] <- NA_real_
    dtf[[duration_v]][maartenutils::eps(dtf[[conc_v]], 0)] <- 0.0
    # idxs <- which(
    #   is.na(dtf[[conc_v]]) | maartenutils::eps(dtf[[conc_v]], 0)
    # )
    # dtf[[duration_v]][idxs] <- NA_real_
  }
  return(dtf)
}


separate_duration.Seurat <- function(so) {
  so@meta.data <- separate_duration(so@meta.data)
  return(so)
}
