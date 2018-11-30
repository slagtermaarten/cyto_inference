setwd('~/MirjamHoekstra')
suppressPackageStartupMessages(library(dplyr))

source('~/MirjamHoekstra/R/targets_init.R')
source('~/MirjamHoekstra/R/init.R')

# targets::tar_config_set('script' = '~/MirjamHoekstra/_targets.R')
# pkgs <- c('magrittr', 'dplyr', 'Seurat', 'glue', 'purrr', 'ggplot2',
#   'patchwork', 'maartenutils')

# # targets:::tar_option_get('storage') targets:::tar_option_get('error')
# targets:::tar_option_set(
#   packages = pkgs,
#   storage = 'worker',
#   retrieval = 'worker',
#   # cue = tar_cue(mode = 'never'),
#   error = 'continue',
#   # error = 'abridge',
#   # error = 'stop',
#   workspace_on_error = F
# )
# tar_manifest()

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
    # 'Milo_targets',
    # 'NH_targets',
    # 'NH_M_targets',
    # 'TRANSACT_NH_targets',
    'combined_experiment_targets',
    'PRECISE_targets',
    'GS_ratio_targets',
    'GS_ratio_plot_targets',
    'exp_groupings',
    'scVI_gpu_prep',
    'scVI_targets',
    'harmony_targets',
    'ML_targets',
    NULL
  )) %>%
  purrr::map(~get(.x, env = targets_env)) %>%
  purrr::discard(is.null) %>%
  tarchetypes::tar_hook_outer(
    hook = { source('~/MirjamHoekstra/R/init.R'); .x },
    names = everything()
  )
# list(
#   bulk_targets,
#   gene_selection_targets,
#   sc_metas,
#   plotting_targets,
#   seurat_object_targets,
#   seurat_object_characterization_targets,
#   # combined_experiment_targets,
#   # PRECISE_targets,
#   # GS_ratio_targets,
#   # GS_ratio_plot_targets,
#   # exp_groupings,
#   # scVI_gpu_prep,
#   # scVI_targets,
#   # harmony_targets,
#   # ML_targets,
#   neighbourhood_targets,
#   TRANSACT_PB_targets,
#   TRANSACT_NH_targets,
#   NULL
# ) 
