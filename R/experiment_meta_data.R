#' Add two variables to a metadata table
#'
#' Use sample_name, tnf_conc, ifn_conc, sn_dilution, duration
#' Create 'stim_group' and (through color_by_stim_group) 'col'
#'
merge_sample_annotation <- function(
  dtf, expo = 2, merge_meta = TRUE,
  exp = NULL, 
  sa = targets::tar_read(sample_annotation_exp5029)) {

  # stopifnot('sample_name' %in% colnames(dtf))
  original_sample_order <- dtf$sample_name
  setDT(dtf)
  dtf <- dtf[, !duplicated(colnames(dtf)), with = F]

  if (merge_meta) {
    if (is.null(sa)) {
      sa <- .__NAMESPACE__.$lazydata$sample_annotation
    }
    dtf <- maartenutils::controlled_merge(
      dtf, sa, by_cols = 'sample_name')
  }

  if (!all(dtf$sample_name == original_sample_order)) {
    setkey(dtf, sample_name)
    dtf <- dtf[original_sample_order]
  }

  if (F) {
    dtf <- dtf %>%
      dplyr::mutate(tnf_rank = factor_to_rank(tnf_conc)) %>%
      dplyr::mutate(ifn_rank = factor_to_rank(ifn_conc)) %>%
      dplyr::mutate(sn_rank = factor_to_rank(sn_dilution)) %>%
      # dplyr::mutate(
      # stim_group = gsub('\\d{1,2}h_(.*)', '\\1', sample_name)) %>%
      as.data.table
  } else {
    if (is.null(dtf$stim_group)) browser()
  }

  # dtf[order(tnf_conc), .(tnf_conc, tnf_rank)]
  # dtf[, .(tnf_conc, tnf_rank)]
  # if (F && !is.null(exp) && exp == '5310') {
  #   ## Overwrite stim group for single cells
  #   dtf[hash_tag == '1', stim_group := 'Exposed to T cells in vivo']
  #   dtf[hash_tag == '2', stim_group := 'Unstimulated in vivo']
  #   dtf[hash_tag == '3', stim_group := 'Unstimulated in vitro']
  #   dtf[hash_tag == '4', stim_group := '100_ng_ml_ifng']
  #   dtf[hash_tag == '5', stim_group := '10_ng_ml_tnfa']
  #   dtf[hash_tag %in% c('6', '7'), stim_group := '100_ng_ml_ifng_10_ng_ml_tnfa']
  #   dtf[hash_tag %in% c('5'), tnf_rank := 4]
  #   dtf[hash_tag %in% c('6', '7'), tnf_rank := 4]
  #   dtf[hash_tag %in% c('6', '7'), ifn_rank := 4]
  # }

  dtf <- format_stim_group(dtf)
  # dtf[, levels(stim_group)]

  dtf[is.na(duration), duration := 'Unknown']
  # if (!is.null(exp) && exp == '5310') {
  #   dtf[hash_tag %in% c('1', '2'), duration := 'Unknown']
  #   dtf[hash_tag %in% c('3', '4', '5', '6', '7'), duration := '24']
  # }
  dtf$duration[is.na(dtf$duration)] <- 'Unknown'
  dtf$stim_group[is.na(dtf$stim_group)] <- 'Unknown'
  dtf$sample_type[is.na(dtf$sample_type)] <- 'sc'
  dtf$sample_name[is.na(dtf$sample_name)] <- 'Unknown'

  dtf[pmax(tnf_rank, ifn_rank, sn_rank, na.rm = T) <= 1 &
      sample_origin == 'in_vitro',
    stim_group := 'Unstimulated in vitro']

  dtf <- color_by_stim_group(dtf, expo = expo)
  return(dtf)
}


cleanup_bulk_sample_names <- function(sn) {
  sn %>% 
    unique() %>%
    tolower() %>%
    { gsub('(\\d)ng', '\\1_ng', .) } %>%
    { gsub('0_1', '0.1', .) } %>%
    { gsub('0_01', '0.01', .) } %>%
    { gsub('tnf(?!a)', 'tnfa', ., perl = T) } %>%
    { gsub('ifn(?!g)', 'ifng', ., perl = T) } %>%
    { . }
}


#' Extract meta data from the way Mirjam and the GCF have named the
#' bulk samples and return a \code{data.frame} with the various
#' metadata fields in separate columns
#'
#'
extract_bulk_meta_data_from_sample_names <- function(sn) {
  ## Clean up sample names first
  sn <- cleanup_bulk_sample_names(sn)
  duration <- factor(as.numeric(gsub('(\\d{1,2})h_.*$', '\\1', sn)))
  tnf_conc <- gsub('\\d{1,2}h.*_(\\d{1,2}\\.*\\d*)_ng_ml_tnfa$', 
    '\\1', sn)
  tnf_conc[tnf_conc == sn] <- 0
  tnf_conc <- factor(as.numeric(tnf_conc))
  ifn_conc <- gsub('\\d{1,2}h.*_(\\d{1,2}\\.*0*1*)_ng_ml_ifn.*', 
    '\\1', sn)
  ifn_conc[ifn_conc == sn] <- 0
  ifn_conc <- factor(as.numeric(ifn_conc))
  sn_dilution <- gsub('\\d{1,2}h.*_1_(\\d{1,5})_sn$', '\\1', sn)
  sn_dilution[sn_dilution == sn] <- 0
  sn_dilution <- factor(sn_dilution, levels=rev(c(2, 200, 20000, 0)))
  stim_groups <- gsub('_', ' ', sn) %>%
    { gsub('ng ml', 'ng/ml', .) } %>%
    { gsub('(\\d) (\\d+) sn', '\\1/\\2 SN', .) } %>%
    { gsub('tnfa', 'TNFa', .) } %>%
    { gsub('tnf', 'TNFa', .) } %>%
    { gsub('ifn[y|g]', 'IFNy', .) } %>%
    { gsub('unstim', 'Unstimulated in vitro', .) } %>%
    { gsub('^\\d{1,2}h ', '', .) } %>%
    { gsub(' 0 ng/ml TNFa', '', .) } %>%
    { . }
  if (!all(stim_groups %in% stim_group_levels)) {
    browser()
    stim_groups[!stim_groups %in% stim_group_levels]
  }
  sample_annotation <- data.frame(
    'sample_name' = sn, 
    'duration' = include_na_lev(duration, el = 48),
    stim_group = stim_groups,
    'sn_dilution' = include_na_lev(sn_dilution), 
    'tnf_conc' = include_na_lev(tnf_conc), 
    'ifn_conc' = include_na_lev(ifn_conc),
    sample_origin = include_na_lev(v = 'in_vitro', 
      el = c('in_vivo', 'intratumoral_injection', 'ex_vivo_in_vitro')),
    sample_type = include_na_lev('bulk', el = 'sc')
  )
  rm(stim_groups)
  sample_annotation <- format_stim_group(sample_annotation)
  sample_annotation <- add_condition_name(sample_annotation)
  return(sample_annotation)
}


#' Update the meta data of an object by HTO (combination)
#'
#'
update_meta <- function(dtf, experiment = dtf$exp[1] %||% NULL) {
  stopifnot(!is.null(experiment))

  stopifnot(is.data.frame(dtf))
  dtf <- dplyr::select(dtf, -matches('condition_name'))
  rn <- rownames(dtf)
  sa <- get(glue::glue('sc_{experiment}_sample_annotation'))

  if ('condition_index' %in% colnames(dtf) &&
      'condition_i' %nin% colnames(dtf)) {
    dtf[['condition_i']] <- dtf[['condition_index']]
  }

  if ('hash_tag' %in% colnames(dtf)) {
    sa <- dplyr::select(sa, hash_tag, condition_name)
    dtf <- left_join(dtf, sa, by = 'hash_tag')
  } else if ('condition_i' %in% colnames(dtf)) {
    dtf$condition_name <-
      dplyr::pull(sa, condition_name)[dtf$condition_i]
    # mean(is.na(dtf$condition_name))
  } else if ('HT1' %in% colnames(dtf)) {
    sa <- dplyr::select(sa, HT1, HT2, condition_name)
    dtf <- left_join(dtf, sa, by = c('HT1', 'HT2'))
    # dtf[, c('HT1', 'HT2')]
    dtf %>%
      dplyr::filter(is.na(condition_name)) %>%
      dplyr::select(HT1, HT2)
  } else {
    stop('Unexpected structure of meta data')
  }

  dtf$condition_name <- factor(dtf[['condition_name']],
    levels = unlist(sa[['condition_name']]))
  # if (any(duplicated(levels(sa$condition_name)))) browser()
  if (mean(is.na(dtf$condition_name)) > .75) browser()
  rownames(dtf) <- rn

  return(dtf)
}
