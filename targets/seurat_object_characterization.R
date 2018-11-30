gen_so_c_table <- function(experiments = sc_e) {
  tibble::tibble(
    experiment = experiments,
    f_so = rlang::syms(glue::glue('filtered_so_{experiment}')),
    oc_markers = rlang::syms(glue::glue('characterize_outlying\\
        _clusters_step_{experiment}')),
    fc_so = rlang::syms(glue::glue('filtered_cleaned_so\\
        _{experiment}')),
    fci_so = rlang::syms(glue::glue('filtered_cleaned_informative\\
        _so_{experiment}')),
    fcmi_so = rlang::syms(glue::glue('filtered_cleaned_mono_\\
        informative_so_{experiment}')),
    fc_magic_so = rlang::syms(
      glue::glue('filtered_cleaned_magic_so_{experiment}'))
  )
}


seurat_object_characterization_targets <- list(
  tarchetypes::tar_map(
    values = gen_so_c_table(c(sc_e, '6743')),
    names = experiment,

    tar_target(simple_UMAPs, {
      p_width <- ifelse(experiment == '5310', 8.7, 17.4)
      s_var <- ifelse(experiment == '5310', NULL, 'duration')
      p1 <- DimPlot(
        fc_so, group.by = 'stim_group', label.size = 2,
        split.by = s_var,
        label = FALSE,
        raster = TRUE) +
        # scale_colour_stim_group(so@meta.data) +
        theme_cyto_inf() +
        guides(colour = guide_legend(ncol = 2))
      print_plot_eval(plot(p1),
        width = p_width, height = 10,
        filename = file.path(Sys.getenv('img_dir'),
          glue::glue('exp{experiment}_UMAP.pdf')))

      p1 <- DimPlot(
        fci_so, group.by = 'stim_group', label.size = 2,
        split.by = s_var,
        label = FALSE,
        raster = TRUE) +
        # scale_colour_stim_group(so@meta.data) +
        theme_cyto_inf() +
        guides(colour = guide_legend(ncol = 2))
      print_plot_eval(plot(p1),
        width = p_width, height = 10,
        filename = file.path(Sys.getenv('img_dir'),
          glue::glue('exp{experiment}_RF_genes_UMAP.pdf')))

      of <- file.path(Sys.getenv('img_dir'),
          glue::glue('exp{experiment}_mono_UMAP.pdf'))
      p1 <- DimPlot(
        fcmi_so, group.by = 'stim_group', label.size = 2,
        split.by = s_var,
        label = FALSE,
        raster = TRUE) +
        # scale_colour_stim_group(so@meta.data) +
        theme_cyto_inf() +
        guides(colour = guide_legend(ncol = 2))
      print_plot_eval(plot(p1),
        width = p_width, height = 10,
        filename = of)
      return(of)
    }, format = 'file'),

    tar_target(read_count_summary, {
      gen_so_read_count_summary(fc_so)
    }),

    tar_target(condition_count_barplot_step, {
      condition_count_barplot(
        sos = list(
          'vanilla' = vanilla_HTO_ann_so,
          "HTO_QC" = HTO_QC_so,
          'lib_QC' = filtered_so,
          'outlier_QC' = filtered_cleaned_so
        ),
        experiment = experiment
      )
    }, format = 'file'),

    tar_target(scree_plot, { plot_scree(fc_so) }, format = 'file'),

    # tar_target(scree_plot_magic, plot_scree(fc_magic_so), format = 'file'))

    tar_target(characterize_outlying_clusters_step, {
      characterize_outlying_clusters(f_so)
    }),

    tar_render(
      name = oc_characterization_report,
      path = file.path(rmd_dir, 'outlying_cluster_report.Rmd'),
      output_format = 'all',
      output_dir = Sys.getenv('reports_dir'),
      output_file = paste0('outlying_cluster_report-experiment=', 
        experiment),
      clean = TRUE,
      params = list(
        'experiment' = experiment,
        'oc_markers' = characterize_outlying_clusters_step,
        'f_so' = f_so
      )
    ),

    tar_target(HT_identifiability_filtered, {
      plot_hashtag_identifiability(
        so = f_so,
        o_fn = gen_HT_identifiability_fn(
          experiment = Project(f_so),
          id = 'filtered'
        )
      )
    }, format = 'file'),

    tar_target(HT_identifiability_filtered_cleaned, {
      if (ncol(fc_so) < ncol(f_so)) {
        plot_hashtag_identifiability(
          so = fc_so,
          o_fn = gen_HT_identifiability_fn(
            experiment = Project(fc_so),
            id = 'filtered_cleaned'
          )
        )
      }
    }, format = 'file')
  ),

  tarchetypes::tar_map(
    values = gen_so_c_table(sc_e),
    names = experiment,

    tar_target(HT_identifiability_filtered_cleaned_informative, {
      so <- fci_so
      plot_hashtag_identifiability(
        so, o_fn = gen_HT_identifiability_fn(
          experiment = Project(fci_so),
          id = 'filtered_cleaned_informative'
        )
      )
    }, format = 'file'),

    tar_target(HT_identifiability_filtered_cleaned_mono_informative, {
      plot_hashtag_identifiability(fcmi_so,
        o_fn = gen_HT_identifiability_fn(
          experiment = Project(fcmi_so),
          id = 'filtered_cleaned_mono_informative')
      )
    }, format = 'file')

    # , tar_target(milo_step, {
    #   source(r_dir, 'milo.R')
    #   run_milo(so = filtered_cleaned_so)
    # }, pattern = map(filtered_cleaned_so), format = 'file'),

    # , tar_target(milo_step_i, {
    #   source(file.path(r_dir, 'milo.R'))
    #   run_milo(so = filtered_cleaned_informative_so, 
    #     fn_app = 'informative')
    # # }, pattern = map(filtered_cleaned_informative_so), 
    # }, pattern = slice(map(filtered_cleaned_informative_so), 2), 
    # format = 'file'),

    # , tar_target(milo_step_im, {
    #   source(r_dir, 'milo.R')
    #   run_milo(so = filtered_cleaned_mono_informative_so, 
    #     fn_app = 'mono_informative')
    # }, pattern = map(filtered_cleaned_mono_informative_so), 
    # format = 'file'),

  )

  # , tar_target(experiment_descriptive_table, {
  #   rlang::syms(c(
  #       glue::glue('bulk_{bulk_e}_so'), 
  #       glue::glue('filtered_cleaned_so_{sc_e}'))) %>%
  #   purrr::map_dfr(function(so) {
  #     so <- read_preproc_experiments(
  #       experiments = experiment, 
  #       sc_mode = 'SCT'
  #     )[[1]]
  #   })
  # })
)
