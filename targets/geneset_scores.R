gs_param_table <-
  tibble(experiment = sc_e) %>%
  dplyr::mutate(
    i = 1:n()
  ) %>%
  tidyr::expand_grid(
    cn_mode = c('all', 'simple', 'vivo', 'vitro_simple',
      'vitro_extended', 'sc_digestion', 'Ag'),
    p_values = c('reference')
  ) %>%
  dplyr::mutate(
    # all_univariate_gs_comps =
    #   rlang::syms(glue::glue('all_univariate_gs_comps_{report_cp}')),
    report_cp = glue::glue('{experiment}_{cn_mode}_{p_values}'),
    unfiltered_so = rlang::syms(glue::glue('HTO_QC_so_{experiment}')),
    so = rlang::syms(
      glue::glue('filtered_cleaned_informative_so_{experiment}')),
    so_score_obj = rlang::syms(glue::glue('so_score_{experiment}')),
    all_fns_obj = rlang::syms(glue::glue('all_fns_{report_cp}')),
    GS_p_dat_obj = rlang::syms(glue::glue('GS_p_dat_{report_cp}')),
    selected_univariate_gs_comps_obj =
      rlang::syms(
        glue::glue('selected_univariate_gs_comps_{report_cp}')),
    all_univariate_gs_comps_obj =
      rlang::syms(glue::glue('all_univariate_gs_comps_{report_cp}')),
    GS_contrast_heatmap_obj =
      rlang::syms(glue::glue('GS_contrast_heatmap_{report_cp}'))
  ) %>%
  dplyr::filter(!(experiment != '5310' &
    cn_mode %in% c('vivo', 'vitro_simple', 'vitro_extended',
      'sc_digestion'))) %>%
  dplyr::filter(!(experiment != '6369' & cn_mode %in% c('frozen'))) %>%
  dplyr::filter(!(experiment != '6493' & cn_mode %in% c('Ag'))) %>%
  # dplyr::left_join(
  #   dplyr::select(NH_t_grid, experiment = query, matches('agg|DA')),
  #   by = 'experiment'
  # ) %>%
  { . }


GS_targets <- tarchetypes::tar_map(
  names = experiment,
  unlist = T,
  values = dplyr::distinct(gs_param_table, experiment, .keep_all = T),

  tar_target(
    name = so_score,
    command = {
      so_s <- ComputeGeneSetScores(
        so = so,
        weighting_funcs = all_weighting_funcs['unweighted'],
        genesets = all_genesets,
        tV = apply(tM_human, 1, max),
        simplify_names = F
      )
      so_s <-
        min_max_scale_assay(so_s, assay = 'GS') %>%
        order_condition_name() %>%
        order_duration() %>%
        simplify_condition_name(
          new_cn = 'cn_simple',
          rep_regex = ' - SC digest| - frozen'
        )
      return(so_s)
    },
    iteration = 'list'
  ),

  NULL
)


GS_report_targets <- tarchetypes::tar_map(
  names = all_of(c('experiment', 'cn_mode', 'p_values')),
  unlist = T,
  values = gs_param_table,

  tar_target(
    name = GS_p_dat,
    command = {
      # attach(gs_param_table[1, ])
      cns <- cns_cn_mode[[cn_mode]]
      # so_score <- tar_read(so_score_obj[[1]])
      p_dat <- create_GS_p_dat(so = so_score_obj, meta_vars = cns)
      ann_p_dat <- create_ann_p_dat(p_dat, cns = cns)

      ## Get the violins to print in the same order as the annotation
      for (cn in colnames(ann_p_dat)) {
        p_dat[[cn]] <- factor(
          p_dat[[cn]],
          levels = unique(ann_p_dat[[cn]])
        )
        ann_p_dat[[cn]] <- factor(
          ann_p_dat[[cn]],
          levels = unique(ann_p_dat[[cn]])
        )
      }
      ann_p_dat$grp_i <- factor(1:nrow(ann_p_dat),
        levels = nrow(ann_p_dat):1)
      # p_dat$grp_i <- factor(p_dat$grp_i, levels = rev(levels(p_dat$grp_i)))
      p_dat <- inner_join(p_dat, ann_p_dat)
      return(p_dat)
    }
  ),

  tar_target(
    name = selected_univariate_gs_comps,
    command = {
      all_fns <-
        rownames(so_score_obj[['GS']]) %>%
        stringr::str_replace('-', '_') %>%
        intersect(colnames(GS_p_dat))

      out <- compare_conditions_by_gs(
        dtf = GS_p_dat,
        gs = all_fns, 
        comp_levs = comp_levs_lu[[experiment]]
      )

      return(out)
    }
  ),

  tar_target(
    name = all_univariate_gs_comps,
    command = {

      all_fns <-
        rownames(so_score_obj[['GS']]) %>%
        stringr::str_replace('-', '_') %>%
        intersect(colnames(GS_p_dat))

      all_comps <- 
        tidyr::expand_grid(
          cn1 = levels(GS_p_dat$condition_name),
          cn2 = levels(GS_p_dat$condition_name)
        ) %>%
        dplyr::filter(cn1 != cn2) %>%
        purrr::pmap(function(cn1, cn2) {
          tibble('condition_name' = c(cn1, cn2))
        })

      if (interactive() && !test_rendering())
        source(file.path(r_dir, 'AUC_comp_gs.RA'))
      out <- compare_conditions_by_gs(
        dtf = GS_p_dat,
        gs = all_fns, 
        comp_levs = all_comps
      )

      return(out)
    }
  ),

  tar_target(
    name = GS_violins,
    command = {

      if (interactive())
        source('~/MirjamHoekstra/R/init.R')

      all_fns <- rownames(so_score_obj[['GS']]) %>%
        stringr::str_replace('-', '_')  %>%
        intersect(colnames(GS_p_dat))

      fns <- all_fns %>%
        c(.[1], .) %>%
        purrr::imap(function(fn, i) {
          mk_ann_GS_vln(
            GS_p_dat = GS_p_dat,
            experiment = experiment,
            fn = fn,
            p_values = p_values,
            cn_mode = cn_mode,
            plot_lgd = i == 1
          )
        })

      return(fns[[length(fns)]])
    }
  ),

  # tar_target(
  #   name = violin_legends,
  #   command = {
  #     p_dat <- create_GS_p_dat(so_score)
  #     ann_p_dat <- create_ann_p_dat(
  #       p_dat = p_dat,
  #       cns = c('stim_group', 'duration')
  #     )
  #     fns <- make_violin_legends(ann_p_dat, 'normal')

  #     ann_p_dat <- create_ann_p_dat(
  #       p_dat = p_dat,
  #       cns = c('stim_group', 'duration', 'frozen', 'sc_digestion')
  #     )
  #     fns <- make_violin_legends(ann_p_dat, 'full')

  #     if (experiment == '5310') {
  #       cns <- c('stim_group', 'sample_origin')
  #       ann_p_dat <- create_ann_p_dat(cns) %>%
  #         dplyr::nest_by(sample_origin)
  #       map2(ann_p_dat$data, ann_p_dat$sample_origin,
  #         ~make_legends(.x, .y))
  #     }

  #     return(fns[1])
  #   },

  #   format = 'file'
  # ),

  tarchetypes::tar_render(
    name = fig4_report,
    path = file.path(rmd_dir, 'fig4.Rmd'),

    output_format = 'all',
    output_dir = Sys.getenv('reports_dir'),
    output_file = paste0('fig4-', experiment, '-', cn_mode, '-',
      p_values),
    clean = TRUE,

    params = list(
      'experiment' = experiment,
      'so' = so,
      'GS_p_dat' = GS_p_dat_obj,
      'cn_mode' = cn_mode,
      'all_univariate_gs_comps' = all_univariate_gs_comps_obj,
      'report_cp' = report_cp,
      'p_values' = p_values
    )
  ),

  # tarchetypes::tar_render(
  #   name = GS_report,
  #   path = file.path(rmd_dir, 'GS_report.Rmd'),

  #   output_format = 'all',
  #   output_dir = Sys.getenv('reports_dir'),
  #   output_file = paste0('GS_report-', experiment, '-', cn_mode, '-',
  #     p_values),
  #   clean = TRUE,

  #   params = list(
  #     experiment = experiment,
  #     cn_mode = cn_mode,
  #     p_values = p_values,
  #     so_score = so_score_obj,
  #     all_fns = all_fns_obj,
  #     GS_p_dat = GS_p_dat_obj,
  #     univariate_gs_comps = univariate_gs_comps_obj,
  #     GS_contrast_heatmap = GS_contrast_heatmap_obj
  #   )
  # ),

  NULL
)
