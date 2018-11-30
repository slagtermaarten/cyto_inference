#' ROC separation test of two numeric vectors
sep_test <- function(x, y) {
  out <- yardstick::roc_auc_vec(
    truth = factor(c(rep('a', length(x)), rep('b', length(y)))),
    estimate = c(x, y)
  )
  if (out < .5) {
    out <- 1 - out
  }
  return(list(p.value = out))
}

# sep_test(x = runif(10, .5, 1), y = runif(10, .8, 1))
# sep_test(x = runif(10, .5, 1), y = runif(10, .5, 1))


make_all_tgs_scatter_plots <- function(...)
  UseMethod('make_all_tgs_scatter_plots')


make_all_tgs_scatter_plots.data.frame <-
  function(
    dtf, 
    colour_var = 'duration',
    facet_var = NULL,
    comps = tribble(
      ~x_var, ~y_var,
      'IFNy', 'TNFa',
      # 'IFNy', 'synergy',
      # 'TNFa', 'synergy',
      # 'IFNy', 'IFNy.late',
      # 'IFNy', 'IFNy.plateau',
      # 'limma_MR1', 'limma_MR2',
      # 'limma_MR3', 'limma_MR4',
      # 'limma_MR4', 'limma_MR5',
      # 'TNFa.6h.max', 'TNFa.early',
      'IFNy.plateau', 'IFNy.late',
      # 'TNFa', 'TNFa.early',
      'TNFa.early', 'TNFa.late'
    ),
    add_centroids = TRUE,
    add_marginal = TRUE,
    marginal_group_color = TRUE,
    marginal_group_fill = FALSE,
    balanced_axes = FALSE,
    point_alpha = 1
  ) {
  pacman::p_load('ggrastr')

  if (!'p_title' %in% colnames(comps)) {
    comps$p_title <- rep('', nrow(comps))
  }

  plots <- comps %>%
    purrr::pmap(function(x_var, y_var, p_title) {
      medians <- dtf %>%
        dplyr::group_by_at(colour_var) %>%
        dplyr::summarize(across(all_of(c(x_var, y_var)), median))
      # p <- purrr::exec(plot_scatter_cor, dtf = dtf,
      #   x_var = x_var, y_var = y_var) +
      p <- ggplot(dtf, 
        aes_string(x = x_var, y = y_var, colour = colour_var)) +
        geom_point(size=.5, alpha=point_alpha) + 
        # ggrastr::rasterise(geom_point(size=.5, alpha=point_alpha), 
        #   dpi = 300) +
        # geom_smooth(method = 'loess', se = T, color = 'black')
        geom_smooth(method = 'loess', se = T)
      if (add_centroids) {
        p <- p + 
          geom_point(
            data = medians, size = 5,
            show.legend = F, colour = 'white',
            ) +
          geom_point(
            data = medians, size = 4,
            # shape = 1,
            show.legend = F, colour = 'black',
            ) +
          geom_point(
            data = medians, size = 3.2,
            show.legend = F, colour = 'white',
            ) +
          geom_point(
            data = medians, size = 3,
            show.legend = F
            ) +
          # geom_text(
          #   mapping = aes_string(x = x_var, y = y_var,
          #     label = colour_var),
          #   colour = 'white',
          #   data = medians, size = 2.4, show.legend = F
          # ) +
          geom_text(
            mapping = aes_string(x = x_var, y = y_var,
              label = colour_var),
            colour = 'black',
            data = medians, size = 2.1, show.legend = F)
      }
      p <- p + ggtitle(p_title)
      p <- p + guides(
        colour = guide_legend(override.aes = list(size=3, alpha = 1))
      )
      if (!is.null(facet_var) && !add_marginal) {
        p <- p + facet_wrap(facet_var)
      }
      if (balanced_axes) {
        ranges <- range(c(dtf[[x_var]], dtf[[y_var]]))
        p <- p + 
          scale_x_continuous(limits = ranges) +
          scale_y_continuous(limits = ranges)
      }
      if (add_marginal) {
        p <- add_marginal_density(p, 
          group_fill = marginal_group_fill,
          group_color = marginal_group_color
        )
      }
      return(p)
    })

  return(plots)
}


add_marginal_density <- function(p, group_color = F, group_fill = F) {
  test_fn <- file.path(Sys.getenv('img_dir'),
    glue::glue('test.pdf'))
  print_plot_eval({
    p <<- ggExtra::ggMarginal(p, type = 'density',
      groupColour = group_color, groupFill = group_fill)
  }, width = 10, height = 10,
  filename = test_fn)
  file.remove(test_fn)
  return(p)
}


make_all_tgs_scatter_plots.Seurat <- function(so,
  stim_group = NULL,
  duration = NULL) {
  if (is.null(stim_group) + is.null(duration) != 1)
    stop('Stim group OR duration, not both')
  active_var <- ifelse(is.null(stim_group),
    'duration', 'stim_group')
  colour_var <- setdiff(c('duration', 'stim_group'), active_var)
  plots <- tribble(~x_var, ~y_var,
    'IFNy', 'TNFa',
    'IFNy', 'synergy',
    'TNFa', 'synergy',
    'IFNy', 'IFNy_late',
    'IFNy', 'IFNy_plateau',
    'IFNy_plateau', 'IFNy_late',
    'TNFa', 'TNFa_early',
    'TNFa_early', 'TNFa_late') %>%
    purrr::pmap(function(x_var, y_var) {
      # install.packages('ggExtra')
      p_dat <- subset_stim_extract_gs(so,
        stim_group = stim_group, duration = duration)
      purrr::exec(plot_scatter_cor, dtf = p_dat,
        x_var = x_var, y_var = y_var, colour_var = colour_var) %>%
        ggExtra::ggMarginal(type = 'density',
          groupColour = T, groupFill = T)
    })
  return(plots)
}


plot_ratio_HM <- function(obj) {
  if (class(obj) == 'SummarizedExperiment') {
    library(SummarizedExperiment)
    M <- t(assays(obj)$ratio)
    ra <- as_tibble(colData(obj))
    experiment <- metadata(obj)$experiment
  } else if (class(obj) == 'data.frame') {
    M <-
      obj %>% dplyr::select(where(is.numeric)) %>%
      as.matrix()
    ra <-
      obj %>% dplyr::select(where(~!is.numeric(.x)))
    experiment <- obj$exp[1]
  }

  row_ann_f <-
    if (any(c('stim_group', 'condition_name') %in% colnames(ra)))
      gen_sample_annotation_HM else function(...) NULL

  H <- gen_HM(
    M,
    ra = ra,
    value_name = 'Gene set ratio',
    col_ann_f = function(...) NULL,
    column_title = glue::glue('Gene set ratios (log2FC)'),
    row_title = glue::glue('Samples (exp {obj$exp[1]})'),
    show_row_dend = T,
    show_column_names = T,
    row_ann_f = row_ann_f
  )

  # {attr(obj, \'fn_mod\') %||% \'\'}\\
  o_fn <- file.path(exp_plot_dir(experiment),
    glue::glue('raw_ratios-{experiment}.pdf'))
  print_plot_eval(
    draw(H, merge_legend = T, newpage = F),
    filename = o_fn)
}


plot_vln <- function(so, fn, assay = NULL,
  # sig_star_ref = "Unstimulated in vitro - 24h",
  sig_star_ref = NULL,
  add_median = T,
  group_var = 'condition_name', ...) {

  p <- VlnPlot(so, features = fn, assay = assay,
    group.by = group_var,
    combine = F, ...)[[1]] +
    theme_cyto_inf(legend.position = 'none') +
    ylab('Gene set score') +
    xlab('') +
    scale_fill_condition_name(so) +
    ggtitle(fn)

  if (add_median) {
    p <- p + stat_summary(fun = median, color = 'grey90')
  }

  if (!is.null(sig_star_ref)) {
    p <- p +
      ggpubr::stat_compare_means(
        label = "p.signif", method = "t.test",
        ref.group = sig_star_ref,
        label.x.npc = .5,
        label.y.npc = .1,
        colour = 'red'
      )
      # ggpubr::stat_compare_means(
      #   label = "p.signif", method = "t.test",
      #   ref.group = "100 ng/ml IFNy 10 ng/ml TNFa - 24h",
      #   label.x.npc = .5,
      #   label.y.npc = .9,
      #   colour = 'red'
      # )
  }

  return(p)
}


plot_vln <- function(
  p_dat,
  fn = 'TNFa',
  # sig_star_ref = "Unstimulated in vitro - 24h",
  sig_star_ref = NULL,
  add_median = T,
  p_y_loc = quantile(p_dat[[as.character(fn)]], .98) * 1.01,
  p_values = 'reference',
  fill_var = 'stim_group',
  # violin_alpha = .8,
  violin_alpha = 1,
  group_var = 'condition_name', ...) {

  pacman::p_load('ggsignif')
  pacman::p_load('ggbeeswarm')

  p_dat[['y']] <- p_dat[[fn]]
  p_dat[['x']] <- p_dat[[group_var]]
  if (is.factor(p_dat[['x']])) {
    p_dat[['x']] <- droplevels(p_dat[['x']])
  }
  p_dat[['fill_var']] <- p_dat[[fill_var]]

  p <- ggplot(p_dat, mapping = aes(x = x, y = y, fill = fill_var)) +
    theme_cyto_inf(legend.position = 'none') +
    ylab('Gene set score') +
    xlab('') +
    geom_violin(alpha = violin_alpha) +
    geom_jitter(alpha = 1, size = .1, width = .5) +
    # ggbeeswarm::geom_beeswarm(alpha = .4, size = .1) +
    coord_flip() +
    ggtitle(fn)

  if (fill_var == 'stim_group') {
    p <- p + scale_fill_stim_group(p_dat)
  }

  if (add_median) {
    p <- p + stat_summary(fun = median, color = 'grey90')
  }

  cns <- levels(p_dat[['x']])

  if (p_values == 'all') {
    comps <-
      tidyr::expand_grid(
        cn1 = seq_along(cns),
        cn2 = seq_along(cns)
      ) %>%
      dplyr::filter(cn2 > cn1) %>%
      { . }
  } else if (p_values == 'reference') {
    comps <-
      tibble(
        # cn1 = 1,
        # cn2 = setdiff(seq_along(cns), 1)
        cn1 = length(cns),
        cn2 = setdiff(seq_along(cns), length(cns))
      ) %>%
      { . }
  } else {
    comps <- NULL
  }

  if (!is.null(comps) && nrow(comps) > 0) {
    filter_sig_comps <- function(comps) {
      comps %>%
        dplyr::mutate(test = purrr::pmap(., function(cn1, cn2) {
          test_GS_comp(
            p_dat = p_dat,
            fn = as.character(fn),
            group_var = as.character(group_var),
            cn_i = cns[cn1], cn_j = cns[cn2]
          )
        })) %>%
        tidyr::unnest(test) %>%
        dplyr::filter(p.value <= .05) %>%
        dplyr::select(cn1, cn2) %>%
        { . }
    }

    if (stringr::str_detect(p_values, 'sig_only')) {
      comps <- filter_sig_comps(comps)
    }

    comps_p <-
      comps %>%
      dplyr::mutate(across(everything(), ~as.character(cns)[.x])) %>%
      purrr::pmap(~c(.x, .y)) %>%
      { . }

    if (F) {
      p <- p +
        geom_signif(
          comparisons = comps_p,
          colour = 'black',
          test = 'ks.test',
          # y_position = max(p_dat[[as.character(fn)]]) * 1.01,
          # y_position = quantile(p_dat[[as.character(fn)]], .6) * 1.01,
          # y_position = quantile(p_dat[[as.character(fn)]], .9) * 1.01,
          y_position = quantile(p_dat[[as.character(fn)]], .98) * 1.01,
          na.rm = TRUE,
          map_signif_level = TRUE,
          margin_top = 0.00,
          step_increase = 0.02,
          tip_length = 0.00,
          textsize = 3
        )
    } else if (F) {
      pacman::p_load('ggpubr')
      # browser()
      p <- p + 
        ggpubr::stat_compare_means(
          ref.group = "1",
          # comparisons = comps_p,
          # label = "p.signif",
          # test = 'ks.test'
          # method = sep_test
          method = 'sep_test'
        )
    } else {
      roc_p_vals <- map_dbl(comps_p, function(x) {
        y <- p_dat[[as.character(fn)]]
        sep_test(
          y[which(p_dat[['x']] == x[1])],
          y[which(p_dat[['x']] == x[2])]
        )$p.value %>%
        round(2)
      })

      if (stringr::str_detect(p_values, 'reference')) {
        p <- p +
          geom_text(
            data = tibble(
              x = map_chr(comps_p, ~.x[2]), 
              p = roc_p_vals),
            mapping = aes(
              x = x, y = p_y_loc, label = p
            ), colour = 'black', size = 3
          )
      } else {
        p <- p +
          geom_signif(
            comparisons = comps_p,
            colour = 'black',
            test = sep_test,
            y_position = p_y_loc,
            na.rm = TRUE,
            # map_signif_level = function(x) {
            #   # if (x < .7)
            #   # round(x * 10)
            #   round(x, 2)
            # },
            # margin_top = 0.00,
            step_increase = 0.05,
            tip_length = 0.005,
            textsize = 3
          )
      }
    }
  }

 return(p)
}


plot_GS_scatter <- function(
  so,
  assay = stringr::str_subset(names(so@assays), 'GS')[2],
  x_var = 'IFNy',
  y_var = 'IFNy-24',
  colour_var = 'IFNy') {

  var_rename <- function(x) str_replace(x, '-', '')
  stopifnot(c(x_var, y_var, colour_var) %in% rownames(so[[assay]]))
  DefaultAssay(so) <- assay
  p_dat <-
    FetchData(so,
      vars = c(x_var, y_var, colour_var, 'condition_name')) %>%
    rename_with(~str_replace(.x, 'gstime_', '')) %>%
    rename_with(~str_replace(.x, 'gscomps_', '')) %>%
    rename_with(~var_rename(.x))

  p <- p_dat %>%
    ggplot(
      aes_string(
        x = glue::glue('{var_rename(x_var)}'),
        y = glue::glue('{var_rename(y_var)}'),
        colour = var_rename(colour_var))) +
    geom_point() +
    theme_cyto_inf(legend.position = 'bottom') +
    scale_colour_viridis_c() +
    facet_wrap(~condition_name, nrow = 3) +
    ggtitle('')

  maartenutils::plot_panel_layout(
    plots = list(p),
    plot_direct = test_rendering(),
    labels = NULL, nrow = 1,
    ncol = 1,
    filename = file.path(img_dir,
      glue::glue('time_geneset_scores_scatters\\
        {make_flag(assay)}\\
        {make_flag(x_var)}{make_flag(y_var)}\\
        {make_flag(colour_var)}.pdf'))
  )
}


GS_umap_panel <- function(query_obj, project_on_ref = F) {
  ref_obj <- tar_read(ratio_bulk)
  if (project_on_ref) {
    ref_umap <- se2tibble(ref_obj) %>%
      add_umap(column_selector = matches('_vs_')) %>%
      order_duration()
    query_umap <- se2tibble(query_obj) %>%
      add_umap(column_selector = matches('_vs_'),
        umap_obj = attr(ref_umap, 'umap_obj')) %>%
      order_duration()
    umap_emb <- dplyr::bind_rows(ref_umap, query_umap)
  } else {
    umap_emb <-
      bind_rows(
        se2tibble(ref_obj)
        , se2tibble(query_obj) %>%
          dplyr::mutate(
            duration = ifelse(!is.na(duration),
              as.character(duration), 'Unknown'))
        # , se2tibble(pb_query_obj) %>%
        #   dplyr::mutate(duration = 'Unknown') %>%
        #   dplyr::mutate(stim_group = 'Unknown') %>%
        #   { . }
      ) %>%
      add_umap(column_selector = matches('_vs_'))
  }

  umap_emb <- order_duration(umap_emb)

  plot_list <-
    c('stim_group',
      stringr::str_subset(colnames(umap_emb), '_vs')) %>%
    maartenutils::auto_name() %>%
    map(~plot_l_umap(dtf = umap_emb, colour_var = .x)) %>%
    map(~.x + guides(shape = 'none') + ggtitle(''))
    # ggpmisc::geom_segment(
    #   aes(x = 5, y = 30, xend = 3.5, yend = 25),
    #   arrow = arrow(length = unit(0.5, "cm")))

  ## UMAP arrows
  # g <- lineGrob(y=0, height = unit(0.1, "npc"), vjust=0,
  #               gp=gpar(fill="black"))
  # p <- p + annotation_custom(g,
  #   xmin = min(umap_emb$UMAP1), xmax=Inf,
  #   ymin = min(umap_emb$UMAP2), ymax=Inf
  # )
  # i + geom_segment(aes(x = 5, y = 30, xend = 3.5, yend = 25),
  #                   arrow = arrow(length = unit(0.5, "cm")))
  o_fn <- file.path(exp_plot_dir(metadata(query_obj)$experiment),
    glue::glue('mixing_based_on_ratios\\
      -{metadata(query_obj)$experiment}\\
      {make_flag(project_on_ref)}.pdf'))
  # print_plot_eval(
  #   print(patchwork::wrap_plots(plot_list,
  #       guides = 'auto', ncol = 2)),
  #   filename = o_fn,
  #   width = 18, height = 25)
  maartenutils::plot_panel_layout(
    plot_list, filename = o_fn,
    width = 18, height = 25, ncol = 2, nrow = 2)
}


GS_umap_panel_2 <- function(query_obj, facet_ncol = NULL) {
  library(SummarizedExperiment)
  experiment <- metadata(query_obj)$experiment

  umap_emb <-
    bind_rows(
      se2tibble(query_obj) %>%
        dplyr::mutate(
          duration = ifelse(!is.na(duration),
            as.character(duration), 'Unknown'))
    ) %>%
    apply_data_filtering(quos(Ag == 'Ag-')) %>%
    add_umap(column_selector = matches('_vs_')) %>%
    dplyr::mutate(duration =
      factor(duration, levels =
        naturalsort::naturalsort(unique(duration)))) %>%
    order_condition_name()

  if (experiment %in% c('6489', '6369', '6493') &&
      is.null(facet_ncol)) {
    facet_ncol <- 4
  }

  plot_list <-
    c('stim_group', stringr::str_subset(colnames(umap_emb), '_vs_')) %>%
    maartenutils::auto_name() %>%
    purrr::map(~plot_l_umap(dtf = umap_emb, colour_var = .x)) %>%
    purrr::map(~.x + guides(shape = 'none') + ggtitle('')) %>%
    purrr::map_at(-1,
      ~.x + facet_wrap(~condition_name, ncol = facet_ncol))
    # ggpmisc::geom_segment(
    #   aes(x = 5, y = 30, xend = 3.5, yend = 25),
    #   arrow = arrow(length = unit(0.5, "cm")))

  ## UMAP arrows
  # g <- lineGrob(y=0, height = unit(0.1, "npc"), vjust=0,
  #               gp=gpar(fill="black"))
  # p <- p + annotation_custom(g,
  #   xmin = min(umap_emb$UMAP1), xmax=Inf,
  #   ymin = min(umap_emb$UMAP2), ymax=Inf
  # )
  # i + geom_segment(aes(x = 5, y = 30, xend = 3.5, yend = 25),
  #                   arrow = arrow(length = unit(0.5, "cm")))

  plot_panel_layout(plot_list,
    filename = file.path(img_dir,
      glue::glue('mixing_based_on_ratios2-{experiment}.pdf')),
    width = 18, height = 25, nrow = 2, ncol = 1)
}


make_all_gs_scatters <- function(...)
  UseMethod('make_all_gs_scatters')


make_all_gs_scatters.character <- function(sc_experiment) {
  so <- tar_read(
    filtered_cleaned_so,
    branch = experiment2index(sc_experiment))[[1]]
  make_all_gs_scatters(so)
}


make_all_gs_scatters.Seurat <- function(so) {
  if ('Ag' %in% colnames(so@meta.data)) {
    so <- so[, so@meta.data$Ag == 'Ag-']
  }
  experiment <- unique(so@meta.data$exp)
  stopifnot(!is.null(experiment) || !is.na(experiment))

  so <- order_duration(so)
  so <- order_condition_name(so)

  clean_string <- function(x) x %>%
    stringr::str_replace_all(' ', '_') %>%
    stringr::str_replace_all('/', '_') %>%
    stringr::str_replace_all('-', '_') %>%
    tolower()

  fns <-
    levels(so@meta.data$stim_group) %>%
    purrr::map_chr(function(stim_group) {
      plots <- make_all_tgs_scatter_plots(
        so = so, stim_group = stim_group, duration = NULL)
      stim_group <- clean_string(stim_group)
      o_fn <- file.path(exp_plot_dir(experiment),
        glue::glue('mirjam_time_geneset_scores\\
          {make_flag(experiment)}\\
          {make_flag(stim_group)}.pdf'))
      plot_panel_layout(
        plots = plots,
        plot_direct = test_rendering(),
        labels = NULL, nrow = 3, ncol = 2,
        filename = o_fn
      )
      return(o_fn)
    })

  fns <-
    levels(so@meta.data$duration) %>%
    purrr::map_chr(function(duration) {
      plots <- make_all_tgs_scatter_plots(
        so = so, stim_group = NULL, duration = duration)
      duration <- clean_string(duration)
      o_fn <- file.path(exp_plot_dir(experiment),
          glue::glue('mirjam_time_geneset_scores\\
            {make_flag(experiment)}\\
            {make_flag(duration)}.pdf'))
      plot_panel_layout(
        plots = plots,
        plot_direct = test_rendering(),
        labels = NULL, nrow = 3, ncol = 2,
        filename = o_fn
      )
      return(o_fn)
    })

  return(fns[1])
}


create_GS_p_dat <- function(so, assay = 'GS', meta_vars = c()) {
  meta_vars <- intersect(
    union(meta_vars, c('condition_name', 'cn_simple',
      'sc_digestion', 'frozen', 'Ag', 'Ag-', 'Ag+',
      'sample_origin', 'ifn_conc', 'sn_dilution',
      'percent.mt', 'tnf_conc', 'duration')),
    colnames(so@meta.data)
  )

  p_dat <-
    cbind(t(so[[assay]][,]), FetchData(so, meta_vars)) %>%
    as.data.frame() %>%
    order_stim_group() %>%
    order_condition_name() %>%
    order_duration()
  if (all(is.na(p_dat$frozen)))
    p_dat$frozen <- NULL

  colnames(p_dat) <- stringr::str_replace_all(
    colnames(p_dat), '-', '_')
  return(p_dat)
}


create_ann_p_dat <- function(p_dat, cns) {
  ann_p_dat <-
    p_dat %>%
    dplyr::distinct(across(any_of(cns)),
      .keep_all = F) %>%
    dplyr::arrange(across(any_of(cns))) %>%
    set_rownames(NULL) %>%
    dplyr::mutate(y = 1:n()) %>%
    { . }

  for (cn in colnames(ann_p_dat)) {
    if (length(unique(ann_p_dat[[cn]])) == 1) {
      ann_p_dat[[cn]] <- NULL
    }
  }

  ann_p_dat$y <- NULL

  return(ann_p_dat)
}


extract_ann_p_dat <- function(p_dat, cns) {
  p_dat %>%
    # dplyr::arrange(across(any_of(c('duration', 'stim_group')))) %>%
    dplyr::distinct(across(any_of(c(cns, 'i')))) %>%
    dplyr::arrange(across(everything()))
}


make_violin_legends <- function(
  ann_p_dat, experiment, ann_type = 'normal', plot_direct = FALSE) {
  ann_p_dat$y <- 1:nrow(ann_p_dat)

  rev_y <- scale_y_continuous(
    name = '', breaks = c(), expand = c(0, 0), limits = rev)
  rev_y <- scale_y_reverse(
    name = '', breaks = c(), expand = c(0, 0))

  make_generic <- function(fill_var = 'sc_digestion') {
    ggplot(ann_p_dat, aes_string(x = '1', y = 'y', fill = fill_var)) +
      geom_tile(color = "white", lwd = 0, linetype = 1) +
      coord_fixed() +
      scale_fill_discrete(name = 'SC digestion') +
      scale_x_continuous(name = '', breaks = c(), expand = c(0, 0)) +
      rev_y +
      theme(plot.margin = unit(c(5, 0, 0, 0), 'mm'),
        legend.position = 'bottom') +
      guides(fill = guide_legend(ncol = 1))
  }

  if ('duration' %in% colnames(ann_p_dat) &&
      (T || experiment %in% c('6369', '6489', '6493', '6601'))) {
    a1 <- ggplot(ann_p_dat, aes_string(x = '1', y = 'y', fill = 'duration')) +
      geom_tile(color = "white", lwd = 0, linetype = 1) +
      coord_fixed() +
      scale_fill_duration(ann_p_dat) +
      scale_x_continuous(name = '', breaks = c(), expand = c(0, 0)) +
      rev_y +
      theme(plot.margin = unit(c(5, 0, 0, 0), 'mm'),
        legend.position = 'bottom') +
      guides(fill = guide_legend(ncol = 1))
  } else {
    a1 <- NULL
  }

  if ('stim_group' %in% colnames(ann_p_dat) &&
      (T || experiment %in% c('5310', '6369', '6489', '6493',
          '6600', '6601'))) {
    a2 <- ggplot(ann_p_dat, aes_string(x = '1', y = 'y',
        fill = 'stim_group')) +
      geom_tile(color = "white", lwd = 0, linetype = 1) +
      coord_fixed() +
      scale_fill_stim_group(ann_p_dat) +
      scale_x_continuous(name = '', breaks = c(), expand = c(0, 0)) +
      rev_y +
      theme(plot.margin = unit(c(5, 0, 0, 0), 'mm'),
        legend.position = 'bottom') +
      guides(fill = guide_legend(ncol = 1))
  } else {
    a2 <- NULL
  }

  if ('sc_digestion' %in% colnames(ann_p_dat)) {
    a3 <- make_generic(fill_var = 'sc_digestion')
  } else {
    a3 <- NULL
  }

  if ('frozen' %in% colnames(ann_p_dat)) {
    a4 <- make_generic(fill_var = 'frozen')
  } else {
    a4 <- NULL
  }

  if ('Ag' %in% colnames(ann_p_dat)) {
    a5 <- make_generic(fill_var = 'Ag')
  } else {
    a5 <- NULL
  }

  plots <- list(a1, a2, a3, a4, a5) %>% purrr::discard(is.null)

  if (plot_direct) {
    print_plot_eval(print(wrap_plots(plots, nrow = 1)),
      width = 17.4, height = .5 * nrow(ann_p_dat) + 2,
      filename = file.path(Sys.getenv('img_dir'),
        glue::glue('test_ann_{ann_type}.pdf')))
  } else {
    return(plots)
  }
}


test_GS_comp <- function(p_dat, fn, group_var, cn_i, cn_j) {
  idxs_i <- which(p_dat[[group_var]] == cn_i)
  idxs_j <- which(p_dat[[group_var]] == cn_j)
  if (length(idxs_i) == 0 || length(idxs_j) == 0)
    return(NULL)
  # tidy(wilcox.test(p_dat[[fn]][idxs_i], p_dat[[fn]][idxs_j]))
  pacman::p_load('broom')
  tidy(ks.test(p_dat[[fn]][idxs_i], p_dat[[fn]][idxs_j]))
}


mk_ann_GS_vln <- function(p_dat, cn_mode, experiment,
  fn = 'IFNy', 
  p_values = 'reference_sig_only',
  out_dir = file.path(Sys.getenv('img_dir'),
    glue::glue('violins_{experiment}')),
  plot_lgd = FALSE, return_mode = 'file', 
  fn_add = '',
  p_y_loc = quantile(p_dat[[fn]], 1, na.rm = T),
  ...) {

  group_var <- group_var_cn_mode[[cn_mode]]
  if (any(!c(group_var, fn) %in% colnames(p_dat))) {
    return(NULL)
  }
  cns <- cns_cn_mode[[cn_mode]]
  if (cn_mode == 'sc_digestion') {
    p_dat <- p_dat %>%
      dplyr::filter(stim_group == '100 ng/ml IFNy 10 ng/ml TNFa')
  } else if (cn_mode == 'frozen') {
    ## Complete later
    p_dat <- p_dat %>%
      dplyr::filter(stim_group == '100 ng/ml IFNy 10 ng/ml TNFa')
  }
  ann_p_dat <- extract_ann_p_dat(p_dat, cns = cns) %>%
    # order_stim_group() %>%
    dplyr::arrange(across(any_of('stim_group')))

  if (!'grp_i' %in% colnames(p_dat)) {
    # p_dat$grp_i <- NULL
    ann_p_dat$grp_i <- factor(1:nrow(ann_p_dat), 
      levels = nrow(ann_p_dat):1)
    p_dat <- p_dat %>% left_join(ann_p_dat, 
      by = setdiff(colnames(ann_p_dat), 'grp_i'))
  }

  legends <- make_violin_legends(
    ann_p_dat = ann_p_dat,
    experiment = experiment
  )
  if (F && interactive()) {
    pe <- patchwork::wrap_plots(
      legends, nrow = 1, guides = 'collect'
    )
    print_plot_eval(print(pe),
      # width = fw * 2.54, height = fh * 2.54,
      force_file = TRUE,
      filename = file.path(out_dir,
        glue::glue('testA.pdf')))
  }

  p <-
    plot_vln(
      p_dat = p_dat,
      fn = fn,
      group_var = 'grp_i',
      fill_var = 'stim_group',
      assay = 'GS',
      p_values = p_values,
      p_y_loc = p_y_loc,
      ...
    ) +
    scale_x_discrete(breaks = c()) +
    theme()
  if (F && interactive()) {
    print_plot_eval(print(p),
      width = 12, height = 10,
      filename = file.path(Sys.getenv('img_dir'),
        glue::glue('test_{fn}.pdf')))
  }

  if (plot_lgd) {
    plots <- legends
    pe <- patchwork::wrap_plots(plots, nrow = 1)
  } else {
    plots <- c(legends, list(p))
    plots <- purrr::map(plots, ~.x + theme(legend.position = 'none'))
    # sapply(plots, class)
    pe <- patchwork::wrap_plots(
      plots,
      nrow = 1,
      # guides = 'collect',
      widths = rep(.075, length(legends)) %>%
        { c(., 1-sum(.)) }
    )
    attr(pe, 'N_legends') <- length(legends)
  }

  if (return_mode == 'p') {
    return(pe)
  } else if (return_mode == 'file') {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    of <- print_plot_eval(
      print(pe),
      height = 10, width = 20,
      force_file = TRUE,
      filename = file.path(out_dir,
        glue::glue('vln_plot_{cn_mode}{fn_add}\\
          _{fn}_{make_flag(p_values)}{make_flag(plot_lgd)}.pdf'))
    )
    return(of)
  }
}


univariate_gs_auc_plots <- function(
  univariate_gs_comps, 
  expected_conditions = unique(c(univariate_gs_comps$cn1,
    c(univariate_gs_comps$cn2))),
  o_fn = file.path(out_dir,
    glue::glue('auroc_all_conditions{make_flag(report_cp)}.pdf'))) {

  fns <- unique(univariate_gs_comps$gs)

  auc_plots <- map(fns, function(feat) {
    p_dat <-
      univariate_gs_comps %>%
      dplyr::filter(gs == feat) %>%
      dplyr::filter(cn1 %in% expected_conditions) %>%
      dplyr::filter(cn2 %in% expected_conditions) %>%
      # dplyr::bind_rows(
      #   dplyr::rename(univariate_gs_comps, cn1 = cn2, cn2 = cn1)
      # ) %>%
      # dplyr::mutate(
      #   cn1 = droplevels(factor(cn1,
      #   levels = levels(GS_p_dat$condition_name)))
      # ) %>%
      # dplyr::mutate(
      #   cn2 = droplevels(factor(cn2,
      #   levels = levels(GS_p_dat$condition_name)))
      # ) %>%
      dplyr::mutate(
        cn1 = droplevels(factor(cn1,
        levels = expected_conditions))
      ) %>%
      dplyr::mutate(
        cn2 = droplevels(factor(cn2,
        levels = rev(expected_conditions)))
      ) %>%
      dplyr::arrange(cn1, cn2) %>%
      dplyr::filter((length(expected_conditions) - as.integer(cn2) + 1) >= as.integer(cn1)) %>%
      { . }
    p <-
      ggplot(data = p_dat, aes(x = cn1, y = cn2,
          fill = AUROC, 
          # label = signif(AUROC, 3))) +
          label = round(AUROC, 2))) +
      geom_tile() +
      geom_text(data = dplyr::filter(p_dat, AUROC >= .75), 
        hjust = .5, size = 2, color = 'grey10') +
      geom_text(data = dplyr::filter(p_dat, AUROC < .75), 
        hjust = .5, size = 2, color = 'grey90') +
      scale_x_discrete(name = '', expand = c(0, 0)) +
      scale_y_discrete(name = '', expand = c(0, 0)) +
      # scale_fill_viridis_c(
      #   direction = -1L,
      #   limits = c(0.5, 1),
      #   breaks = seq(0, 1, by = .1)) +
      scale_fill_continuous(
        limits = c(0.5, 1),
        breaks = seq(0, 1, by = .1)) +
      rotate_x_labels(45) +
      ggtitle(feat)
    if (F) {
      pwidth <- 10
      pheight <- 10
      of <- file.path(
        out_dir,
        paste0('GS_contrast_heatmap_', report_cp, '-', feat, '.pdf')
      )
      print_plot_eval(print(p),
        width = pwidth,
        height = pheight,
        filename = of
      )
    }

    return(p)
  })

  plot_panel_layout(
    plots = auc_plots,
    ncol = 1, nrow = 2,
    clear_redundant_legends = F,
    width = length(expected_conditions) * 1.1 + 3, 
    height = 2*(length(expected_conditions) * 1.1 + 3), 
    filename = o_fn, labels = NULL
  )

  return(invisible(o_fn))
}


compute_all_pairwise_medians <- function(
  GS_p_dat,
  comp_vars = names(which(map_lgl(GS_p_dat, is.numeric))),
  expected_conditions = unique(GS_p_dat$condition_name)) {
  
  out <- 
    tidyr::expand_grid(
      gs = comp_vars,
      cn1 = factor(expected_conditions),
      cn2 = factor(expected_conditions)
    ) %>%
    dplyr::filter(as.integer(cn1) > as.integer(cn2)) %>%
    { . }

  tr <- function(x) log2(median(x, na.rm = T) + 1)

  stats <- 
    purrr::pmap_dfr(out, function(gs, cn1, cn2) {
      v <- GS_p_dat[[gs]]
      wc_test <- broom::tidy(wilcox.test(
        v[GS_p_dat$condition_name == cn1],
        v[GS_p_dat$condition_name == cn2]))
      list('log2FC_med' = 
        tr(median(v[GS_p_dat$condition_name == cn1], na.rm = T)) - 
        tr(median(v[GS_p_dat$condition_name == cn2], na.rm = T))
      ) %>%
      c(wc_test)
    }) %>%
    dplyr::select(-method, -alternative) %>%
    dplyr::mutate(percentual_relative_to_smallest = 100*(2^abs(log2FC_med)-1)) %>%
    dplyr::mutate(med_fc = percentual_relative_to_smallest / 100 + 1)

  return(bind_cols(out, stats))
}

mean_expression_heatmaps <- function(GS_p_dat, o_fn) {
  all_fns <- names(which(map_lgl(GS_p_dat, is.numeric)))
  plots <-
    all_fns %>%
    map(function(feat) {
      # feat = 'CD274'
      # feat = 'WARS'
      # feat <- all_fns[1]

      t_dat <-
        # cbind(FetchData(so, feat), so@meta.data[,]) %>%
        GS_p_dat %>%
        dplyr::group_by_at(group_vars) %>%
        dplyr::summarize(
          feat = median(.data[[feat]]),
          perc_exp = mean(.data[[feat]] > 0)
          ) %>%
        { . }

      p <- ggplot(t_dat, aes(x = stim_group, y = duration,
          fill = feat,
          # size = perc_exp,
          # width = sqrt(perc_exp),
          # height = sqrt(perc_exp)
          label = round(feat, 2),
          )) +
      rotate_x_labels(45) +
      geom_tile() +
      geom_text(size = 3) +
      # scale_fill_viridis_c(direction = -1L, name = feat) +
      scale_x_discrete(name = 'Stimulus', expand = c(0, 0)) +
      scale_y_discrete(name = 'Duration', expand = c(0, 0)) +
      ggtitle(feat) +
      guides(fill = 'none')

    if (length(remaining_group_vars) > 0) {
      p <- p + facet_wrap(remaining_group_vars)
    }
    # print_plot_eval(print(p),
    #   width = 8.7, height = 10,
    #   filename = file.path(out_dir,
    #     glue::glue('expression_HM{make_flag(feat)}.pdf')))
    return(p)
  })

  plot_panel_layout(
    plots = plots,
    # ncol = ifelse(experiment %in% c('6369', '6493'), 3, 2),
    ncol = 2,
    nrow = 3,
    clear_redundant_legends = F,
    width = 17.4, height = 25, filename = o_fn, labels = NULL
  )
}

