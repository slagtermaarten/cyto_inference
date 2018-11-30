#' Reward 'differing' neighbors, but do not penalize identical
#' neighbors
#'
#'
std_mismatch <- function(cn, levs) {
  levs <- maartenutils::auto_name(levs)
  list(
    'cn' = 'exp',
    'subs_table' = 
      map(levs, function(lev) {
        map(levs, ~ifelse(.x != lev, 1, 0))
      })
  )
}


gen_gs_func <- function(subs_table) {
  stopifnot(!is.null(subs_table))

  ## Complete substitution table if necessary
  obs_names <- names(subs_table)
  if ('default' %in% obs_names) {
    ## Expected names, infer from the fields in the 'default' entry
    exp_names <- setdiff(names(subs_table[['default']]), obs_names)
    additions <- map(auto_name(exp_names), ~subs_table[['default']])
    subs_table <- subs_table[setdiff(obs_names, 'default')] %>%
      c(additions)
  }

  gs_func <- function(e, o, w) {
    stopifnot(length(o) == length(w))
    # unique(o[!o %in% names(subs_table)])
    scores <- vapply(seq_along(o), function(i) {
      subs_table[[e]][[o[i]]] * w[i]
    }, numeric(1))
    return(sum(scores))
  }
  # error_func('5310', c('6434', '5310'), c(.4, 4))

  return(gs_func)
}


gs_subs_tables <- list(
  'sample_type_mixing' =
    list(
      'cn' = 'sample_type',
      'subs_table' = list(
        'sc' = list('sc' = 0, 'bulk' = 1),
        'bulk' = list('bulk' = 0, 'sc' = 1)
      )
    ),
  'sc_contact' = list(
     'cn' = 'sample_type',
     'subs_table' = list('default' = list('sc' = 1, 'bulk' = 0))
  ),
  'vitro_sc_contact' = list(
     cn = 'group_id',
     subs_table = list(
       'default' =
         list('sc-in_vitro' = 1, 'sc-in_vivo' = 0, 'bulk-in_vitro' = 0)
      )
    ),
  'vivo_sc_contact' = list(
    cn = 'group_id',
    subs_table = list('default' =
       list('sc-in_vitro' = 0, 'sc-in_vivo' = 1, 'bulk-in_vitro' = 0)
      )
    ),
  'experiment' =
    std_mismatch(
      cn = 'exp',
      levs = c('5029', '5310', '6369', '6434', '6489', '6493')
    ),
  'plain_duration' =
    std_mismatch(cn = 'duration', levs = c(2, 6, 12, 24)),
  'weighted_duration' =
    list(
      cn = 'duration', 
      subs_table = list(
        '2' = list('2' = 0, '6' = 1, 
          '12' = 2, '24' = 3, 'Unknown' = 0),
        '6' = list('2' = 1, '6' = 0, 
          '12' = 1, '24' = 2, 'Unknown' = 0),
        '12' = list('2' = 2, '6' = 1, 
          '12' = 0, '24' = 1, 'Unknown' = 0),
        '24' = list('2' = 3, '6' = 2, 
          '12' = 1, '24' = 0, 'Unknown' = 0)
      )
    )
)

compute_gs <- function(so, gs_mod = 'sc_contact') {
  ## Setup all requirements first
  SNN_name <- setdiff(str_extract(names(so@graphs), '.*snn'), NA)
  M <- so@graphs[[SNN_name]]
  if (is.null(M)) stop('No integrated SNN available')
  ## Extract colname to look for values in meta data
  cn <- gs_subs_tables[[gs_mod]][['cn']]
  subs_table <- gs_subs_tables[[gs_mod]][['subs_table']]
  f <- gen_gs_func(subs_table = subs_table)
  ## Extract the meta data column
  y <- so@meta.data[[cn]]

  print(table(y))

  ## Compute graph scores
  map_dbl(1:nrow(M), function(.x) {
    w <- M[.x, ]
    f(e = y[.x], o = y, w = w) / sum(w, na.rm = T)
  })
}
