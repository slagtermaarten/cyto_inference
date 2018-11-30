source('~/MirjamHoekstra/.Rprofile')


# python_loc <- '/DATA/users/m.slagter/miniconda3/envs/r4/bin/python3'
# reticulate::use_condaenv('r4', required=TRUE)
# reticulate::py_config()

# pacman::p_load(pkgs, character.only = T)
# devtools::install(file.path('~/libs', 'maartenutils'))
suppressWarnings(suppressMessages(suppressPackageStartupMessages({
  pacman::p_load(pkgs, character.only = TRUE)

  if (F) {
    devtools::document('~/libs/DISTINCT', quiet = TRUE)
  }
  devtools::load_all('~/libs/DISTINCT')
  # Sys.getenv()
  source(file.path(Sys.getenv('r_dir'), 'helpers.R'))
  source(file.path(Sys.getenv('r_dir'), 'define_target.R'))
  source(file.path(Sys.getenv('r_dir'), 'targets_init.R'))
  source(file.path('targets', '1-param_grid.R'))
  source(file.path(libs_dir, 'result_cacher.R'))
  source(file.path(Sys.getenv('r_dir'), 'S3.R'))
  # source(file.path(Sys.getenv('r_dir'), 'symbol_capping.R'))
  # library(envnames)
  source(file.path(Sys.getenv('r_dir'), 'format_stim_group.R'))
  source(file.path(Sys.getenv('r_dir'), 'assess_clusters.R'))
  source(file.path(Sys.getenv('r_dir'), 'genesets.R'))
  # if (!'cytokineinference' %in% loadedNamespaces()) {
  #   if (T) {
  #     devtools::load_all(package_dir)
  #   } else {
  #     ## Devtools can be so horribly slow (and why?!)! Try to avoid
  #     ## devtools::load_all if possible Please don't judge me Hadley
  #     R_files <- list.files(package_r_dir, pattern = '\\.R', full.names = T)
  #     for (f in R_files) {
  #       source(f)
  #     }
  #   }
  # }
  # sessionInfo()
  # source(file.path(r_dir, 'read_data.R'))
  # source(file.path(r_dir, 'merge_annotation.R'))
  # source(file.path(r_dir, 'pca.R'))
  # source(file.path(r_dir, 'assess_embedding.R'))
  # source(file.path(r_dir, 'GPR_estimator.R'))
  # source(file.path(r_dir, 'classify_transcriptome.R'))
  # source(file.path(r_dir, 'graph_classification.R'))
  source(file.path(Sys.getenv('r_dir'), 'RunTMM.R'))
  source(file.path(Sys.getenv('r_dir'), 'experiment_meta_data.R'))
  source(file.path(Sys.getenv('r_dir'), 'explore_clusters.R'))
  # source(file.path(analysis_dir, 'transact', 'transact_funcs.R'))
  source(file.path(Sys.getenv('r_dir'), 'plotting.R'))
  # source(file.path(Sys.getenv('r_dir'), 'gsea.R'))
  source(file.path(Sys.getenv('r_dir'), 'data_integration.R'))
  source(file.path(Sys.getenv('r_dir'), 'preprocessing_sc.R'))
  source(file.path(Sys.getenv('r_dir'), 'time_informativeness.R'))
  source(file.path(Sys.getenv('r_dir'), 'limma.R'))
  source(file.path(Sys.getenv('r_dir'), 'mouse.R'))
  source(file.path(Sys.getenv('r_dir'), 'select_genes.R'))
  source(file.path(Sys.getenv('r_dir'), 'geneset_scores.R'))
  source(file.path(Sys.getenv('r_dir'), 'geneset_plots.R'))
  source(file.path(Sys.getenv('r_dir'), 'ML.R'))
  source(file.path(Sys.getenv('r_dir'), 'ML_plot.R'))
  source(file.path(Sys.getenv('r_dir'), 'R_PRECISE.R'))
  # source(file.path(Sys.getenv('r_dir'), 'TRANSACT.R'))
  source(file.path(Sys.getenv('r_dir'), 'kPCA.R'))
  source(file.path(Sys.getenv('r_dir'), 'milo.R'))
  # source(file.path(Sys.getenv('r_dir'), 'query2ref.R'))
  source(file.path(Sys.getenv('r_dir'), 'AUC_comp_gs.R'))
  # source(file.path(Sys.getenv('r_dir'), 'compare_TRANSACT_settings.R'))

  # source(file.path(Sys.getenv('r_dir'), 'scvi_hyperopt_result_analysis.R'))
  # source(file.path(Sys.getenv('r_dir'), 'tidymodels_helpers.R'))
  # source(file.path(Sys.getenv('r_dir'), 'tidymodels_regression.R'))
  # source(file.path(Sys.getenv('r_dir'), 'scVI.R'))
  # source(file.path(Sys.getenv('r_dir'), 'harmony.R'))
  # source(file.path(Sys.getenv('r_dir'), 'latent_distance.R'))
  # source(file.path(Sys.getenv('r_dir'), 'binarize_scores.R'))
  # source(file.path(Sys.getenv('r_dir'), 'splicing.R'))
})))

Sys.setenv(HDF5_USE_FILE_LOCKING='FALSE')
# Sys.getenv('HDF5_USE_FILE_LOCKING')


source('~/libs/maartenutils/R/string_formatting.R')
library(tibble)
