suppressPackageStartupMessages(library(dplyr))

pacman::p_load('targets')
pacman::p_load('tarchetypes')

library(targets)

targets::tar_config_set('script' = '~/MirjamHoekstra/_targets.R')
# targets::tar_config_set('script' = '_targets.R')

# targets:::tar_option_get('storage')
# targets:::tar_option_get('error')
targets:::tar_option_set(
  packages = pkgs,
  storage = 'worker', 
  retrieval = 'worker',
  # cue = tar_cue(mode = 'never'),
  error = ifelse(interactive(), 'stop', 'continue'), 
  # error = 'abridge', 
  # error = 'stop', 
  # workspace_on_error = ifelse(interactive(), F, T)
  workspace_on_error = T
  # workspace_on_error = F
)

do_par = F
do_par = T
if (!interactive() && exists('do_par') && do_par) {
  pacman::p_load('future')
  pacman::p_load('future.callr')
  plan(callr)
  # pacman::p_load('future')
  # plan(multicore)
}


tar_traceback_raw <- function(name) {
  cl <- match.call()
  cl[[1L]] <- quote(tar_traceback)
  cl[[2L]] <- name
  eval.parent(cl)
}


tar_error <- function(name) {
  targets::tar_meta() %>%
    dplyr::filter(name %in% .env[['name']]) %>%
    dplyr::pull(error)
}


get_last_errored_target <- function() {
  last_errored_target <- 
    tryCatch({
      tar_progress() %>%
      dplyr::filter(progress == 'errored') %>%
      pull(name) %>%
      { .[length(.)] }
    }, error = function(e) { c() }) 

  if (length(last_errored_target) == 0) {
    rlang::warn('Did not find an errored target')
    return(invisible(NULL))
  } else {
    return(last_errored_target)
  }
}


last_target_error <- function() {
  last_errored_target <- get_last_errored_target()
  print(tar_error(last_errored_target))
  print(tar_traceback_raw(last_errored_target))
}
# last_target_error()


find_targets <- function(pattern) {
  stringr::str_subset(tar_meta()$name, pattern = pattern)
}


show_targets_meta <- function(pattern = 'harmonized_Ms') {
  tar_dtf <- tar_meta()
  idxs <- which(tar_dtf$name %in% find_targets(pattern))
  if (length(idxs) == 0) {
    rlang::warn('No targets were found')
  } else {
    return(tar_dtf[idxs, ])
  }
}


show_targets_progress <- function(pattern = 'harmonized_Ms') {
  tar_dtf <- 
    tar_progress(fields = c('type', 'progress', 'mtime', 'time')) %>% 
    { . }
  idxs <- which(tar_dtf$name %in% find_targets(pattern))
  if (length(idxs) == 0) {
    rlang::warn('No targets were found')
  } else {
    return(tar_dtf[idxs, ])
  }
}


smart_progress <- function(by_var_part = TRUE, pattern = NULL) {
  byv <- ifelse(by_var_part, 'variable parts', 'constants parts')
  cat('Printing by', byv, '\n')

  wide_progress <- 
    tar_progress(fields = c('type', 'progress', 'mtime', 'time')) %>% 
    { . }

  if (!is.null(pattern)) {
    # pattern = '6493'
    wide_progress <- 
      wide_progress %>%
      dplyr::filter(stringr::str_detect(name, pattern))
  }

  wide_progress <-
    wide_progress %>%
    dplyr::mutate(stem = 
      rm_var_parts(name, retain_id = by_var_part)) %>%
    dplyr::mutate(progress = factor(progress, 
        levels = c('skipped', 'started', 'built', 'errored'))) %>%
    # dplyr::mutate(tNA = map_lgl(name, ~test_NA(tar_read_raw(.x)))) %>%
    dplyr::group_by(stem, progress) %>%
    dplyr::summarize(total = n()) %>%
    tidyr::pivot_wider(
      names_from = progress, 
      values_from = total,
      values_fill = 0
    ) %>%
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    { . }

  if ('skipped' %in% colnames(wide_progress)) {
    if ('started' %in% colnames(wide_progress)) {
      wide_progress <- wide_progress %>%
        dplyr::filter(!(skipped == total & started == 0L))
    } else {
      wide_progress <- wide_progress %>%
        dplyr::filter(!(skipped == total))
    }
  }

  if ('errored' %in% colnames(wide_progress)) {
    wide_progress$frac_error <- 
      with(wide_progress, errored / total)
    wide_progress <- arrange(wide_progress, total, frac_error)
  }

  # print(wide_progress, n=400L, width=200L)
  # print(as.data.frame(wide_progress))
  return(wide_progress)
}


load_targets <- function() {
  library(dplyr)
  fo <- maartenutils::gen_file_overview(
    '~/MirjamHoekstra/targets',
    include_full = T, pat = '*\\.R$') %>%
    dplyr::arrange(short_fn)
  # targets_env <- new_environment()
  targets_env <- rlang::env()
  for (fn in fo$full_fn) { 
    if (interactive()) message(fn)
    tryCatch({
      suppressWarnings(source(fn, local = targets_env))
    }, error = function(e) { 
      print(paste0('Problem loading ', fn)); print(e) 
    }) 
  }
  stopifnot(length(ls(targets_env)) > 5L)
  # targets_env$sc_metas_targets
  # for (fn in fo$full_fn) { source(fn) }
  return(targets_env)
}
# targets_env <- load_targets()
# ls(targets_env)


if (interactive()) {
  fo <- maartenutils::gen_file_overview('~/MirjamHoekstra/targets',
    include_full = T, pat = '*\\.R$') %>%
    dplyr::arrange(short_fn)
  do_load_targets <- 
    !exists('targets_env') || 
    (exists('last_load_time') && any(fo$mtime > last_load_time))
} else {
  do_load_targets <- TRUE
}


if (do_load_targets) {
  # tryCatch({
  # }, error = function(e) { print('Could not load targets'); print(e) }) 
  targets_env <- load_targets()
  last_load_time <- Sys.time()
  ls(targets_env)
}


if (F) {
  fo <- 
    maartenutils::gen_file_overview('~/MirjamHoekstra/targets',
      include_full = T, pat = '*\\.R$') %>%
    dplyr::arrange(short_fn)
  for (fn in fo$full_fn) { 
    tryCatch({
      source(fn) 
    }, error = function(e) { 
      print(paste0('Problem loading ', fn)); print(e) 
    }) 
  }
}

# stopifnot(!maartenutils::null_dat(tar_read(bulk_comb_so)))
