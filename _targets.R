setwd('~/MirjamHoekstra')
suppressPackageStartupMessages(library(dplyr))

source('~/MirjamHoekstra/R/targets_init.R')
source('~/MirjamHoekstra/R/init.R')

if (F) {
  ls(pattern = '_targets$', targets_env) %>%
    purrr::map(~get(.x, env = targets_env)) %>%
    unlist(recursive = F) %>%
    { .[which(map_lgl(., ~!'tar_stem' %in% class(.x)))] } %>%
    { . }
}

## All objects ending with _targets are assumed to contained
## tar_target objects or lists thereof
# ls(pattern = '.*', targets_env)
# targets_env$NH_M_targets_K
ls(pattern = '_targets$', targets_env) %>%
  setdiff(c(
    'combined_experiment_targets',
    'exp_groupings',
    'ML_targets',
    NULL
  )) %>%
  purrr::map(~get(.x, env = targets_env)) %>%
  purrr::discard(is.null) %>%
  tarchetypes::tar_hook_outer(
    hook = { source('~/MirjamHoekstra/R/init.R'); .x },
    names = everything()
  )
