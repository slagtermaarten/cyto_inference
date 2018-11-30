library(ggplot2)

source(file.path(r_dir, 'result_cacher.R'))

# pacman::p_load('showtext')
# showtext_auto()

if (F) {
  # pacman::p_load('extrafont')
  # font_import()
  # font_add_google("Lato", "lato")
  # fonts()
  # font_paths()
  # font_files()
  # loadfonts(device = "pdf")
  # fonts()
}

okabe <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
  '#D55E00', '#CC79A7')
options(ggplot2.discrete.fill = okabe)
options(ggplot2.discrete.colour = okabe)

heatmap_fs <- 6
heatmap_font <- NULL
heatmap_font <- 'Arial'
heatmap_font <- 'ArialMT'
heatmap_font <- 'Helvetica'
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt(
  heatmap_row_names_gp = gpar(fontsize = heatmap_fs,
    fontfamily = heatmap_font),
  heatmap_row_title_gp = gpar(fontsize = heatmap_fs,
    fontfamily = heatmap_font),
  heatmap_column_names_gp = gpar(fontsize = heatmap_fs,
    fontfamily = heatmap_font),
  heatmap_column_title_gp = gpar(fontsize = heatmap_fs,
    fontfamily = heatmap_font),
  legend_title_gp = gpar(fontsize = heatmap_fs,
    fontface = 'italic', fontfamily = heatmap_font),
  legend_labels_gp = gpar(fontsize = heatmap_fs,
    fontface = 'italic', fontfamily = heatmap_font),
  legend_border = NA,
  annotation_border = F,
  heatmap_border = F
)
ht_opt$HEATMAP_LEGEND_PADDING = unit(1, 'cm')

library(maartenutils)
theme_cyto_inf <- function(base_size = 8, ...) {
  maartenutils::theme_ms(
    base_size = base_size,
    base_family = heatmap_font
  ) +
  theme(
    strip.background = ggplot2::element_rect(
      fill = 'white', colour = 'grey80', size = 0.5),
    panel.border = ggplot2::element_rect(
      colour = 'grey20', fill = NA, size = .5,
      linetype = 'solid'),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    legend.key.width = unit(6, 'mm'),
  ) +
  theme(...)
}

ggplot2::theme_set(theme_cyto_inf())

cyto_inf_cols <- maartenutils::gen_color_vector(10)


#' Encode marginal data selection in a string, for use in filenames or plot
#' titles
#'
#'
marginal_sel_to_character <- function(marginal_sel) {
  items <- imap(unlist(as.list(marginal_sel)), ~ sprintf('%s_%.1f', .y, .x))
  paste0('-', items, collapse = '-')
}


plot_expression_dynamics <- result_cacher(
  min_mod_time = '2022-05-21 00:00',
  f = function(
    fn = 'CD274',
    p_dat = NULL,
    leave_out_sn = FALSE,
    axis_name = NULL,
    meta_data = tar_read(sample_annotation_exp5029),
    merge_cn = NULL,
    y_scale = 'linear',
    title = fn,
    feature_mode = NULL,
    lookup_data = targets::tar_read('kallisto_5029')) {

    # if (ncol(lookup_data) <= 70) browser()

    if (is.null(p_dat)) {
      if (is.null(merge_cn))
        find_shared_genes
      if (is.matrix(lookup_data) &&
          !'sample_name' %in% colnames(lookup_data)) {
        if (fn %nin% rownames(lookup_data))
          return(NULL)
        p_dat <-
          lookup_data[which(rownames(lookup_data) == fn), ] %>%
          unlist %>%
          tibble::enframe(merge_cn, 'gexp') %>%
          dplyr::left_join(meta_data, by = merge_cn)

      } else {
        p_dat <- lookup_data
        p_dat$gexp <- p_dat[[fn]]
        p_dat$duration <- factor(p_dat$duration,
          levels = naturalsort::naturalsort(unique(p_dat$duration)))
        print(setdiff(p_dat$stim_group, stim_group_levels))
        p_dat$stim_group <- factor(p_dat$stim_group,
          levels = stim_group_levels)
        stopifnot(!any(is.na(p_dat$stim_group)))
        p_dat <- dplyr::select(p_dat, gexp, duration, stim_group)
      }

      setDT(p_dat)
      # col_pal <- p_dat[, setNames(col, stim_group)]
      # setdiff(as.character(p_dat$stim_group), names(col_pal))
      # setdiff(names(col_pal), as.character(p_dat$stim_group))
    }

    if (leave_out_sn) {
      p_dat <- p_dat %>% dplyr::filter(!grepl('_sn', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('lympho', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('Lympho', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('IL-2', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('il.*2', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('IFNA1', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('ifna1', stim_group))
      p_dat <- p_dat %>% dplyr::filter(!grepl('IFNa1', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('TGFb', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('TGFB', sample_name))
      p_dat <- p_dat %>% dplyr::filter(!grepl('tgfb', sample_name))
    }

    if (!grepl('^PC\\d+$', fn) && !grepl('^MG\\d{1,2}$', fn)) {
      if (y_scale == 'log2') {
        p_dat <- dplyr::mutate(p_dat, gexp = log2(gexp + 1))
        if (missing(axis_name)) {
          axis_name <- axis_labels['log_gene_expression_cpm']
        }
      } else {
        if (missing(axis_name)) {
          axis_name <- axis_labels['gene_expression_cpm']
        }
      }
    }

    med_pt_dat <- p_dat %>%
      dplyr::group_by(duration, stim_group) %>%
      dplyr::summarize(gexp = median(gexp))

    if (is.null(p_dat$col) || any(is.na(p_dat$col))) {
      # col_pal <- gen_cyto_inf_col_scale(p_dat$stim_group)
      col_pal <- defensive_stim_col_scale(p_dat$stim_group)
    } else {
      col_pal <-
        distinct(p_dat, stim_group, .keep_all = T) %>%
        { with(., setNames(col, stim_group)) }
    }

    exp6434_test <-
      any(stringr::str_detect(meta_data$sample_name, '^6434_'))
    if (exp6434_test) {
      med_pt_dat <-
        dplyr::left_join(
          med_pt_dat,
          dplyr::select(meta_data, stim_group, tnf_conc),
          by = 'stim_group'
          ) %>%
      dplyr::mutate(tnf_conc = numeric2factor(tnf_conc))
    }

    p <- ggplot(med_pt_dat,
      aes_string(
        x = 'duration', y = 'gexp',
        color = 'stim_group', group = 'stim_group')) +
      geom_path() +
      geom_point() +
      # rotate_x_labels(rotate_labels=45) +
      scale_colour_manual(
        name = '', values = col_pal) +
      scale_y_continuous(
        name = axis_name, expand = c(0, 0.1)) +
      scale_x_discrete(
        name = axis_labels['duration'], expand = c(0, 0.1)) +
      ggtitle(title) +
      theme(legend.position = 'right', legend.direction = 'vertical')

    if (exp6434_test) {
      p <- p + facet_wrap(~tnf_conc)
    }
    return(p)
  },
  filename = function() {
    source(file.path(r_dir, 'result_cacher.R'))
    if (is.null(p_dat)) {
      fn <- file.path(data_dir, 'expression_dynamics_plots',
         glue::glue('{fn}{make_flag(leave_out_sn)}{y_scale}.rds'))
    } else {
      fn <- NULL
    }
    return(fn)
  }
)


#' Convenience function for plot_expression_dynamics
#'
#'
exp_dynamics_panel <- function(
  features = NULL,
  plot_titles = NULL,
  version = '',
  leave_out_sn = T,
  lookup_data = as.matrix(targets::tar_read('kallisto_5029')),
  meta_data = targets::tar_read(sample_annotation_exp5029),
  merge_cn = 'condition_name',
  y_scale = 'linear',
  feature_mode = NULL,
  # max_plots_per_doc = 47,
  max_plots_per_doc = 1e3,
  caching = T,
  nrow_per_doc = 4,
  ncol_per_doc = 3,
  plot_direct = maartenutils::test_rendering(),
  out_dir = img_dir,
  redo = F) {

  source(file.path(r_dir, 'result_cacher.R'))
  tryCatch(detach('package:rlang'), error = function(e) { })

  fm_check <- !missing(feature_mode) && !is.null(feature_mode)
  ld_check <- !missing(lookup_data) || !is.null(lookup_data)
  arg_count <- as.integer(fm_check) + as.integer(ld_check)
  force(redo)

  if (arg_count == 0) {
    lookup_data <- load_EE_feature_mat()
  } else if (fm_check) {
    lookup_data <-
      prep_embedding_for_plotting(feature_mode, redo = redo)
    lookup_data <- lookup_data[lookup_data$sample_type == 'bulk', ]
    ## 2021-06-28 11:50 BROKEN FOR NOW
  }

  if (is.null(features) || missing(features)) {
    features <- stringr::str_extract(colnames(lookup_data),
      '(MG|UMAP|scVI)_*\\d') %>%
      setdiff(NA)
  }

  missing_genes <- setdiff(features, rownames(lookup_data))
  if (length(missing_genes) > 0) {
    warning('Missing: ', paste(missing_genes, collapse = ', '))
    features <- intersect(features, rownames(lookup_data))
  }
  pca_mode <- all(stringr::str_detect(features, '^(PC|PV|CF)'))

  if (is.null(plot_titles))
    plot_titles <- features

  N_docs <- ceiling(length(features) / max_plots_per_doc)

  purrr::map(seq(N_docs), function(i) {
    subs_start <- max_plots_per_doc * (i-1) + 1
    subs_end <- min(length(features), max_plots_per_doc * i)
    l_features <- features[subs_start:subs_end]
    l_plot_titles <- plot_titles[subs_start:subs_end]

    plots <-
      tibble(fn = l_features, pt = l_plot_titles, i = i) %>%
      # furrr::future_pmap(function(fn, pt, i) {
      pmap(function(fn, pt, i) {
        tryCatch(plot_expression_dynamics(
          fn = fn,
          leave_out_sn = leave_out_sn,
          y_scale = y_scale,
          caching = caching && !pca_mode,
          meta_data = meta_data,
          merge_cn = merge_cn,
          lookup_data = subset_feats(lookup_data, fn),
          redo = redo) + ggtitle(pt),
          error = function(e) { print(e); NULL })
      })

    valid_bool <- purrr::map_lgl(plots,
      ~tryCatch('gg' %in% class(.x) ||
        !is.null(ggpubr::as_ggplot(ggpubr::get_legend(.x))),
      error = function(e) { NULL })
    )
    plots <- plots[valid_bool]

    if (length(plots) == 0) {
      plot_expression_dynamics(
        fn = l_features[1], leave_out_sn = F,
        y_scale = y_scale,
        caching = caching,
        lookup_data = lookup_data, redo = redo)
      return(NULL)
    }

    if (F) {
      cached_legend <-
        ggpubr::as_ggplot(ggpubr::get_legend(plots[[1]]))
    } else {
      cached_legend <- NULL
    }
    plots <- c(plots[1],
      map(plots, ~.x + theme(legend.position = 'none')))
    if (plot_direct) {
      o_fn <- NULL
    } else {
      o_fn <- file.path(out_dir,
        glue('feature_dynamics_{version}{make_flag(y_scale)}{make_flag(leave_out_sn)}_{i}.pdf'))
    }
    maartenutils::plot_panel_layout(
      plots = c(plots, list(cached_legend)),
      plot_direct = plot_direct,
      labels = NULL, nrow = nrow_per_doc, ncol = ncol_per_doc,
      filename = o_fn
    )
  })
}


stim_regexes <- list(
  tnf = '^(?<!\\d{3} ng/ml IFNy )(0.1|1|10) ng/ml TNFa$',
  ifn = '^(0.01|1|100) ng/ml IFNy(?! 10 ng/ml TNFa)$',
  combo = '^(0.01|1|10|100) ng/ml IFNy (0.1|1|10) ng/ml TNFa',
  sn = '.* SN',
  unstim = '(u|U)nstim',
  vivo = '.* in vivo$',
  injection_vivo = '^In vivo.*$',
  T_cell_exposed = '.*(e|E)xposed.*'
)
# v <- naturalsort::naturalsort(unique(dtf[[colour_var]]))
# tibble(value = v, selected = grepl(hl_stim, v, perl = T))


plot_UMAP <- function(dtf, ...) {
  UseMethod('plot_UMAP')
}


#' Plot UMAP low D rep of samples
#'
#'
plot_UMAP.data.frame <- function(dtf,
                                 stim_hl_group = NULL,
                                 hash_label_hl_group = NULL,
                                 shape_var = NULL,
                                 size_var = 'sample_type',
                                 # label_var = 'sample_name',
                                 colour_var = 'stim_group',
                                 pt_alpha = .1, pt_size = .5) {
  if (any(class(dtf) == 'integer')) browser()
  if (is.numeric(dtf$duration) || is.integer(dtf$duration)) {
    dtf$duration[is.na(dtf$duration) | dtf$duration == 'NA'] <-
      'Unknown'
  } else if (is.factor(dtf$duration)) {
    idxs <-
      as.numeric(dtf$duration) == which(is.na(levels(dtf$duration)))
    dtf$duration[idxs] <- 'Unknown'
  }
  dtf <- dtf[!grepl('SN', dtf$stim_group), ]

  if (!is.null(shape_var) && !shape_var %in% colnames(dtf)) {
    message('Skipping the plot for as it does not exist: ',
      shape_var)
    return(NULL)
  }

  if (!is.null(size_var) && !size_var %in% colnames(dtf)) {
    message('Skipping the plot for as it does not exist: ',
      size_var)
    return(NULL)
  }

  if (!is.null(colour_var) && !colour_var %in% colnames(dtf)) {
    message('Skipping the plot for as it does not exist: ',
      colour_var)
    return(NULL)
  }

  levs <- c('Unknown',
    sort(as.numeric(setdiff(unique(dtf$duration), 'Unknown'))))
  dtf$duration <- factor(dtf$duration, levels = levs)
  dtf$sample_type[is.na(dtf$sample_type)] <- 'Unknown'

  if (F && !is.null(stim_hl_group) && stim_hl_group == 'vivo')
    browser()

  if (F && !is.null(colour_var) && colour_var == 'clusters')
    browser()

  if (!is.null(stim_hl_group)) {
    ## For highlighting certain groups in a within experiment type
    ## analysis
    v <- unique(dtf$stim_group)
    # if (stim_hl_group == 'combo') browser()
    allowed <- grep(stim_regexes[[stim_hl_group]], v,
      perl = T, value = T)
    dtf$stim_group[!dtf$stim_group %in% allowed] <- 'background'
    table(dtf$stim_group)
    levs <- c('background',
      setdiff(sort(unique(dtf$stim_group)), 'background'))
    dtf$stim_group <- factor(dtf$stim_group, levels = levs)
    if (stim_hl_group == 'vivo') {
      rn <- rownames(dtf)
      # dtf <- dtf %>% dplyr::filter(stim_group != 'background')
    }
    dtf <- dtf %>% dplyr::arrange(stim_group)
    # rle(as.character(dtf$stim_group))
    col_group <- table(dtf$stim_group)
    if (length(col_group) == 1) {
      return(NULL)
    }
  } else {
    detected_levs <- unique(dtf$stim_group)
    ue_levs <- setdiff(detected_levs, stim_group_levels)
    if (length(ue_levs) > 0) {
      exps <- p_dat %>%
        dplyr::filter(stim_group %in% ue_levs) %>%
        pull(exp)
      dtf %>% group_by(stim_group) %>%
        summarize(N = n(), .groups = 'drop')
    }
    dtf$stim_group <-
      factor(dtf$stim_group, level = stim_group_levels)
  }

  if (!is.null(hash_label_hl_group) &&
      'hash_label' %in% colnames(dtf)) {
    ## For highlighting certain groups in a across experiment type
    ## analysis
    v <- unique(dtf$hash_label)
    allowed <- grep(hash_label_hl_group, v, perl = T, value = T)
    if (length(allowed) == 0) {
      message('No groups match hash_label_hl_group')
      browser()
    }
    dtf$hash_label[!dtf$hash_label %in% allowed] <- 'background'
    levs <- c('background',
              setdiff(sort(unique(dtf$hash_label)), 'background'))
    dtf$hash_label <- factor(dtf$hash_label, levels = levs)
    dtf <- dtf %>% arrange(hash_label)
    col_group <- dtf %>%
      dplyr::group_by(hash_label) %>%
      dplyr::summarize(N = n(), .groups = 'drop')
    if (nrow(col_group) == 1) {
      return(NULL)
    }
  }

  if (F && !is.null(colour_var) && colour_var == 'clusters')
    browser()

  discrete_col_scale <-
    is.character(dtf[[colour_var]]) ||
    is.factor(dtf[[colour_var]])
  if (!is.null(colour_var) && discrete_col_scale) {
    stopifnot(colour_var %in% colnames(dtf))
    if (discrete_col_scale) {
      levs <- setdiff(unique(dtf[[colour_var]]), NA)
      bg_lev <-
        as.character(levs[which(levs %in% c('Unknown', 'background'))])
      levs <- setdiff(levs, bg_lev)
      dtf[[colour_var]] <-
        factor(dtf[[colour_var]], levels = c(bg_lev, levs))
      stopifnot(!all(is.na(dtf[[colour_var]])))
      dtf <- dtf %>% dplyr::filter(!is.na(.data[[colour_var]]))
      dtf <- dtf %>% dplyr::arrange(.data[[colour_var]])
      if ('background' %in% dtf[[colour_var]] &&
          rle(as.character(dtf[[colour_var]]))[[2]][1] !=
        'background') {
        browser()
      }
    } else {
      # stop('Not implemented')
    }
  }

  empty_cols <- which(sapply(dtf, function(x) all(is.na(x))))
  if (length(empty_cols) > 1) {
    dtf <- dtf[, -empty_cols]
  }

  if (!is.null(size_var) && length(unique(dtf[[size_var]])) > 1) {
    if (size_var == 'cluster_representative') {
      dtf <- dtf %>% filter(!is.na(.data[[size_var]]))
      # dtf <- dtf %>% filter(!is.na(.data[[size_var]]))
      dtf[[size_var]][is.na(dtf[[size_var]])] <- FALSE
    }
    if (size_var %nin% colnames(dtf)) browser()
    dtf[[size_var]] <- factor(dtf[[size_var]])
  }

  if (!is.null(shape_var) && length(unique(dtf[[shape_var]]))) {
    dtf[['mylabel']] <- dtf[[shape_var]] %>%
      { .[. == 'Unknown'] <- 'U'; . } %>%
      as.factor()
  } else {
    dtf[['mylabel']] <- NA
  }

  features <- stringr::str_extract(colnames(dtf),
    '(MG|UMAP|scVI)_*\\d') %>%
    setdiff(NA)

  if (F && !is.null(colour_var)) {
    stopifnot(rle(as.character(dtf[[colour_var]]))[[2]][1] ==
      'background')
  }

  p <- ggplot(data = dtf,
    mapping = aes_string(x = features[1], y = features[2],
      shape = shape_var,
      label = 'mylabel',
      colour = colour_var,
      size = size_var)) +
    geom_point(alpha = pt_alpha, show.legend = T, size = pt_size) +
    theme(
      legend.position = 'right',
      legend.direction = 'vertical',
      legend.box.just = 'left'
    )

  if (!is.null(colour_var) &&
      !colour_var %in% c('N_UMI', 'clusters')) {
    p <- p + geom_text(alpha = .7, show.legend = F)
  }

  if (!is.null(colour_var) && colour_var == 'clusters') {
    p_dat$clusters <- factor(p_dat$clusters)
    cl_stats <- plyr::ldply(levels(p_dat$clusters), function(cl) {
      idxs <- which(p_dat$clusters == cl)
      A <- apply(p_dat[idxs, features], 2, mean, na.rm = T)
      SD <- apply(p_dat[idxs, features], 2, sd, na.rm = T)
      out <- tibble(cluster = cl, variable = features, A = A, SD = SD)
      return(out)
    },
    .parallel = parallel,
    .progress = ifelse(parallel, 'none', 'text'))

    cl_centroids <- reshape2::dcast(cl_stats, cluster ~ variable,
      value.var = 'A')
    p <- p + geom_text(data = cl_centroids,
      mapping = aes_string(x = features[1], y = features[2],
        label = 'cluster'),
      inherit.aes = F, alpha = .7, color = 'black')
  }

  if (!is.null(colour_var) && discrete_col_scale) {
    if (colour_var == 'stim_group') {
      if (is.null(stim_hl_group)) {
        levs <- setdiff(dtf$stim_group, 'background')
        # levels(dtf$stim_group)
        col_scale <- gen_cyto_inf_col_scale(levs)
        col_scale <- c('background' = 'grey80', col_scale)
      } else {
        levs <- levels(dtf$stim_group)[-1]
        cbPalette <- c('#E69F00', '#56B4E9', '#009E73',
                       '#F0E442', '#0072B2', '#D55E00', '#CC79A7')
        N <- length(cbPalette)
        N_levs <- length(levs)
        f <-
          circlize::colorRamp2(seq(0, 1, length.out = N), cbPalette)
        lev_cols <- f(seq(0, 1, length.out = N_levs)) %>%
          as.color_vector
        col_scale <- c('background' = 'grey80',
          setNames(lev_cols, levs))
      }
      p <- p + scale_colour_manual(name = '', values = col_scale)
    } else {
      # levs <- levels(dtf[[colour_var]])
      if (colour_var == 'hash_label') {
        # detected_exps <-
        #   setdiff(gsub('([^-]*)-.*', '\\1', dtf$hash_label),
        #     'background')
        # print(exp_count)
        if (length(exp_count) > 1) {
          base_cols <-
            maartenutils::gen_color_vector(
              names(exp_count), 'FantasticFox1')
          interpolated_cols <-
            purrr::imap(exp_count, function(.x, .y) {
              if (.x > 1) {
                pal_colors <- c(
                  maartenutils::darken(base_cols[.y], .5),
                  maartenutils::darken(base_cols[.y], 1.5))
              } else {
                pal_colors <- base_cols[.y]
              }
              colorRampPalette(pal_colors)(.x)
            })
          col_scale <- setNames(unlist(interpolated_cols), levs)
        } else {
          col_scale <- maartenutils::gen_color_vector(
            levs, 'FantasticFox1')
        }
      } else {
        col_scale <- maartenutils::gen_color_vector(arg = levs) %>%
          setNames(levs)
      }
      if (length(bg_lev) > 0) {
        col_scale[bg_lev] <- 'grey80'
      }
      p <- p + scale_colour_manual(values = col_scale)
    }
  }

  if (!is.null(size_var)) {
    if (size_var == 'sample_type') {
      p <- p + scale_size_discrete(
        name = 'Sample type', breaks = c('sc', 'bulk'),
        range = c(4, 1), labels = c('Single-cell', 'Bulk')
      )
    } else if (size_var == 'evenness') {
      p <- p + scale_size_continuous(name = 'Posterior peakness')
    } else if (size_var == 'N_features') {
      dtf[[size_var]] <- as.integer(dtf[[size_var]])
      p <- p + scale_size_continuous(
        name = 'Number of features',
        breaks = sort(unique(dtf[[size_var]]))
      )
    }
  }

  p <- p + gg_legend_alpha_cancel

  return(p)
}


plot_UMAP.Seurat <- function(object, ...) {
  cbind(object@reductions$umap@cell.embeddings, object@meta.data) %>%
    as.data.frame %>%
    plot_UMAP(...)
}


gen_col_fun <- function(M,
  grad_colors = c('skyblue4', 'grey95', 'tomato4'),
  LB = 0, UB = max(M, na.rm = T),
  direction = 1L) {
  if (is.null(M) || all(is.na(M))) return(NULL)
  if (all(as.vector(M) >= 0, na.rm = T)) {
    if (F) {
      grad_breaks <- c(0, max(as.vector(M), na.rm = T))
      if (length(grad_colors) == 3) {
        grad_colors <- grad_colors[c(2, 3)]
      }
      col_fun <- circlize::colorRamp2(
        colors = grad_colors, breaks = grad_breaks)
    } else {
      N_inter = 40L
      grad_colors <- viridis::magma(N_inter, direction = direction)
      grad_colors <- viridis::viridis(N_inter, direction = direction)
      col_fun <- circlize::colorRamp2(
        colors = grad_colors,
        breaks = seq(LB, UB, length.out = N_inter)
      )
    }
  } else {
    max_val <- max(abs(as.vector(M)), na.rm = T)
    grad_breaks <- c(-max_val, 0, max_val)
    col_fun <- circlize::colorRamp2(
      colors = grad_colors, breaks = grad_breaks)
    new_grad_breaks <- sort(c(0, range(M, na.rm = T)))
    col_fun <- circlize::colorRamp2(
      colors = col_fun(new_grad_breaks), breaks = new_grad_breaks)
    # browser()
    # max_val
    # col_softening_f <- 1
    # col_softening_f <- abs(grad_breaks / max_val) %>%
    #   { .[2] = 1; . }
    # grad_colors <- map_chr(seq_along(grad_colors),
    #   ~maartenutils::darken(
    #     grad_colors[.x],
    #     factor = col_softening_f[.x])
    # )
  }
  return(col_fun)
}



#' Gen heatmap of samples and expression
#'
#' @param M Feature matrix, [features x samples]
#'
gen_HM <- function(
  M,
  ca = NULL,
  ra = NULL,
  top_HM = NULL,
  base_size = heatmap_fs,
  show_row_names = (nrow(M) <= 50) &&
    (is.null(N_hl_genes) || N_hl_genes == 0),
  show_column_names = FALSE,
  # N_genes = max(nrow(M), 50),
  # N_hl_genes = ifelse(nrow(M) > 50, nrow(M) / 2, 0),
  N_genes = NULL,
  N_hl_genes = NULL,
  column_ordering = 1:ncol(M),
  # trans = function(x) log10(x + 1),
  trans = NULL,
  cluster_columns = T,
  cluster_rows = T,
  row_names_side = 'right',
  show_column_dend = T,
  column_annotation_name_side = 'right',
  show_row_dend = T,
  # value_name = 'Gene expression',
  # value_name = NULL,
  col_ann_f = gen_sample_annotation_HM,
  row_ann_f = gen_sample_annotation_HM,
  show_heatmap_legend = T,
  use_raster = F,
  raster_by_magick = F,
  col_fun = NULL,
  ann_ordering = NULL,
  show_ca_legend = T,
  show_ra_legend = T,
  highlight_genes = NULL,
  # clustering_distance = 'pearson') {
  clustering_distance_rows = 'euclidean', 
  clustering_distance_columns = 'euclidean', 
  ...) {

  dots <- list(...)

  library(ComplexHeatmap)
  if (!is.null(N_genes)) {
    M <- gene_filtering(M, N_genes = N_genes)
  }
  if (!is.null(trans))
    M <- trans(M)

  # M <- M[, naturalsort::naturalsort(colnames(M))]

  if (F && null_dat(ca)) {
    setDT(combined_sample_annotation)
    setkey(combined_sample_annotation, sample_name)
    ca <- combined_sample_annotation[colnames(M)]
    ca <- merge_sample_annotation(ca, merge_meta = F)
  }

  legend_params <- list(
    title_gp = gpar(fontsize = base_size, fontface = 'italic'),
    direction = 'horizontal',
    labels_gp = gpar(fontsize = base_size),
    grid_height = unit(3, 'mm'),
    grid_width = unit(3, 'mm')
  )

  if (!is.null(col_ann_f) && !maartenutils::null_dat(ca)) {
    CA <- col_ann_f(
      ca, heatmap_type = 'column',
      ann_ordering = ann_ordering,
      show_legend = show_ca_legend
    )
  } else {
    CA <- NULL
  }

  if (!is.null(dots$left_annotation)) {
    RA1 <- dots$left_annotation
    # RA1 <- RA1[rownames(M), ]
    dots$left_annotation <- NULL
  } else if (!is.null(row_ann_f) && 
             !maartenutils::null_dat(ra)) {
    RA1 <- row_ann_f(
      sa = ra, heatmap_type = 'row',
      ann_ordering = ann_ordering,
      show_legend = show_ra_legend)
  } else {
    RA1 <- NULL
  }

  if (!is.null(N_hl_genes) && N_hl_genes != 0 &&
      is.null(highlight_genes)) {
    highlight_genes <-
      apply(M, 1, function(x) var(x, na.rm = T)) %>%
      { rank(.) } %>% { .[. > (length(.) - N_hl_genes)] } %>%
      names
  }

  if (!is.null(dots$right_annotation)) {
    RA2 <- dots$right_annotation
  } else if (!is.null(highlight_genes)) {
    gene_highlight_idx <- match(highlight_genes, rownames(M))

    af <- anno_mark
    RA2 <- rowAnnotation(
      foo = af(
        at = gene_highlight_idx,
        labels = rownames(M)[gene_highlight_idx],
        labels_gp = gpar(fontsize = base_size, fontface = 'italic')
      )
    )
  } else {
    RA2 <- NULL
  }
  dots$right_annotation <- NULL

  if (is.null(col_fun))
    col_fun <- gen_col_fun(M)

  p_args <- list(
    col = col_fun,
    # name = value_name,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    clustering_distance_rows = clustering_distance_rows,
    clustering_distance_columns = clustering_distance_columns,
    row_title_gp = annotation_par,
    show_heatmap_legend = show_heatmap_legend,
    gap = unit(0, 'mm'),
    # row_gap = unit(0, 'mm'),
    # column_gap = unit(0, 'mm'),
    # col = circlize::colorRamp2(colors = c('grey95', col),
    #                            breaks = c(0, max(M, na.rm = T))),
    # row_names_gp = annotation_par,
    column_title_gp = annotation_par,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    show_column_dend = show_column_dend,
    show_row_dend = show_row_dend,
    show_row_names = show_row_names,
    row_names_side = row_names_side,
    show_column_names = show_column_names,
    use_raster = use_raster,
    raster_by_magick = raster_by_magick,
    left_annotation = RA1,
    right_annotation = RA2,
    top_annotation = CA,
    heatmap_legend_param = legend_params
  ) 
  p_args <- p_args[setdiff(names(p_args), names(dots))]
  p_args$cluster_columns
  
  MH <- purrr::exec(Heatmap, M[, column_ordering, drop=F], !!!p_args, !!!dots)

  if (F) {
    if (!is.null(top_HM)) {
      H <- top_HM %v% MH
    } else if (!is.null(RA1) || !is.null(RA2)) {
      H <- RA1 + MH + RA2
    } else {
      H <- MH
    }

    return(H)
  }

  return(MH)
}


color_by_stim_group <- function(dtf, expo = 2) {
  dtf$ifn_rank <- factor_to_rank(dtf$ifn_conc)
  dtf$tnf_rank <- factor_to_rank(dtf$tnf_conc)
  dtf$sn_rank <- factor_to_rank(dtf$sn_dilution)
  # table(dtf$tnf_rank)
  # table(dtf$ifn_rank)
  # table(dtf$sn_rank)

  rn <- rownames(dtf)
  setDT(dtf)
  # dtf[, 'col' data.table::`:=` '']
  dtf[, 'col' := '']
  ## Make color palette
  dtf[, .(ifn_rank, ifn_conc)]
  if (!all(is.na(dtf$ifn_rank)) && any(dtf$ifn_rank)) {
    ifn_cols <- seq(2, max(dtf$ifn_rank, na.rm = T)) %>%
      { 1 + (. / max(.))^expo } %>%
      { maartenutils::lighten('firebrick4', factor=.) } %>%
      setNames(seq(2, max(dtf$ifn_rank, na.rm = T))) %>%
      tibble::enframe('ifn_rank', 'col') %>%
      dplyr::mutate(ifn_rank = as.integer(ifn_rank))
    dtf <- maartenutils::controlled_merge(
      dtf, ifn_cols, by = 'ifn_rank')
  }

  if (!all(is.na(dtf$tnf_rank)) && any(dtf$tnf_rank)) {
    tnf_cols <- seq(2, max(dtf$tnf_rank, na.rm = T)) %>%
      { 1 + (. / max(.))^expo } %>%
      { maartenutils::lighten('dodgerblue4', factor=.) } %>%
      setNames(seq(2, max(dtf$tnf_rank, na.rm = T))) %>%
      tibble::enframe('tnf_rank', 'col') %>%
      dplyr::mutate(tnf_rank = as.integer(tnf_rank))
    dtf <- maartenutils::controlled_merge(
      dtf, tnf_cols, by = 'tnf_rank')
  }

  if (!all(is.na(dtf$sn_rank)) && any(dtf$sn_rank > 1)) {
    sn_cols <- seq(2, max(dtf$sn_rank, na.rm = T)) %>%
      { .5 + (. / max(.))^(expo) } %>%
      { maartenutils::lighten('darkolivegreen4', factor=.) } %>%
      rev %>%
      setNames(seq(2, max(dtf$sn_rank, na.rm = T))) %>%
      tibble::enframe('sn_rank', 'col') %>%
      dplyr::mutate(sn_rank = as.integer(as.character(sn_rank)))
    dtf <- maartenutils::controlled_merge(
      dtf, sn_cols, by = 'sn_rank')
  }
  # dtf[!is.na(sn_rank) & sn_rank != 1,
  #     'col' := sn_cols[as.numeric(sn_rank) - 1]]

  ## Do we need multple combo colors?
  combo_stim_groups <- dtf %>%
    dplyr::filter(ifn_rank > 1 & tnf_rank > 1) %>%
    dplyr::distinct(ifn_rank, tnf_rank, .keep_all = T) %>%
    dplyr::pull(stim_group)
  if (T || length(combo_stim_groups) > 1) {
    combo_cols <-
      maartenutils::auto_name(combo_stim_groups) %>%
      purrr::map_chr(function(csg) {
        ranks <- dtf %>%
          dplyr::filter(stim_group == csg) %>%
          dplyr::select(tnf_rank, ifn_rank) %>%
          dplyr::distinct(.keep_all = T)
        b_cols <- c(
          inner_join(ifn_cols, ranks[, 'ifn_rank'],
            by = 'ifn_rank')$col,
          inner_join(tnf_cols, ranks[, 'tnf_rank'],
            by = 'tnf_rank')$col) %>% unique
        if (length(b_cols) != 2) browser()
        blended_col <-
          circlize::colorRamp2(colors = b_cols, breaks = c(0, 1))(.5)
        if (F) {
          color_vector <- c(b_cols, 'test' = blended_col)
          class(color_vector) <- 'color_vector'
          plot(color_vector)
        }
        return(blended_col)
      }) %>% tibble::enframe('stim_group', 'col')
    if (F) {
      color_vector <- with(combo_cols,
        rlang::set_names(col, stim_group)) %>%
        maartenutils::attr_pass('class', 'color_vector')
      plot.color_vector(color_vector)
    }
    dtf <- maartenutils::controlled_merge(
      dtf, combo_cols, by = 'stim_group')
    stopifnot(combo_cols$col %in% dtf$col)
  } else {
    dtf[tnf_rank == 4 & ifn_rank == 4,
        'col' := 'darkorchid4']
  }

  dtf[stim_group == 'Exposed to T cells in vivo', 'col' := 'aquamarine']
  dtf[stim_group == 'Exposed to T-cells in vivo', 'col' := 'aquamarine']
  dtf[stim_group == 'Unstimulated in vivo' , 'col' := 'aquamarine4']
  dtf[stim_group == 'Unexposed in vivo' , 'col' := 'aquamarine4']
  dtf[stim_group == 'Unstimulated', 'col' := 'grey50']
  dtf[stim_group == 'Unstimulated in vitro', 'col' := 'grey50']

  ## i. All stim groups should have a color
  ## ii. All colors should occur with equal frequency (if bulk seq)
  if (dtf[col == '', .N] > 0 ||
    (nrow(dtf) <= 64 && dtf[, .N, col][, data.table::uniqueN(N) != 1])) {
    browser()
    dtf[col == dtf[, .N, col][N > 1, col]]
    table(dtf$col)
    dtf[col == '']
    # dtf[stim_group == '10 ng/ml TNFa', tnf_rank]
    # dtf[, table(stim_group)]
  }
  dtf <- as.data.frame(dtf)
  rownames(dtf) <- rn
  return(dtf)
}


gen_cyto_inf_col_scale_from_obs <- function(dtf) {
  dtf %>%
    dplyr::select(stim_group, matches('conc|dilution')) %>%
    unique() %>%
    color_by_stim_group() %>%
    dplyr::select(stim_group, col) %>%
    tibble::deframe()
}


gen_cyto_inf_col_scale <- function(
  subset_vals = NULL,
  mirjam_mode = TRUE) {
  primaries <- maartenutils::gen_color_vector(5)
  primaries <- c('firebrick4', 'dodgerblue4', 'darkorchid2',
    'darkolivegreen4', 'aquamarine4', 'orange') %>%
    maartenutils::col_to_hex()
  primaries[4] <- darken(primaries[4], factor = 1.8)

  col_scale <- c(
      'Unstimulated' = '#D6D6D6',
      'Unstimulated in vitro' = '#D6D6D6',
      'background' = 'grey90',
      setNames(maartenutils::gen_col_gradient_vector(
          maartenutils::lighten(primaries[1]), N = 4),
        c('0.01 ng/ml IFNy', '1 ng/ml IFNy',
          '10 ng/ml IFNy', '100 ng/ml IFNy')),
      # setNames(
      #     maartenutils::gen_col_gradient_vector(
      #       '#209CFF',
      #       # 'dodgerblue4',
      #       maartenutils::darken('dodgerblue4', 1.9),
      #       N = 3
      #   ), paste0(c(.1, '1', '10'),  ' ng/ml TNFa')),
      setNames(c('#73D1FF', '#408FFC', '#1149AE'),
        paste0(c(.1, '1', '10'),  ' ng/ml TNFa')),
      setNames(maartenutils::gen_col_gradient_vector(
          primaries[4], offset_col = 'darkolivegreen1', N = 8),
        c('1/20000 SN', '1/625 SN', '1/200 SN',
          '1/125 SN', '1/25 SN', '1/5 SN', '1/2 SN', 'SN')),
      setNames(primaries[5], 'Unexposed in vivo'),
      setNames(maartenutils::lighten(primaries[5]),
        'Exposed to T-cells in vivo'),
      setNames(primaries[6], '5000 U/ml IFNA1'),
      setNames('orange', '100 ng/ml IFNy 10 ng/ml TGFb')
    ) %>%
    # setNames(gen_col_gradient_vector(primaries[4], N = 2),
    #          c('Unexposed in vivo', 'Exposed to T-cells in vivo'))) %>%
    # c(comb_col_scale) %>%
    { . }

  if (mirjam_mode) {
    reps <- c(
      setNames(
        c('#EBB0A6', '#E37D75', '#DC2F3E', '#BA0402'),
        c('0.01 ng/ml IFNy', '1 ng/ml IFNy',
          '10 ng/ml IFNy', '100 ng/ml IFNy')
      ),

      setNames(
        # c('#73D1FF', '#408FFC', '#1149AE'),
        # c('#73b4ff', '#4053FC', '#1149AE'),
        c('#73b4ff', '#0f6ff5', '#1149AE'),
        c('0.1 ng/ml TNFa', '1 ng/ml TNFa', '10 ng/ml TNFa')),

      setNames(c('grey50', '#167700'), c('Unexposed in vivo',
          'Exposed to T-cells in vivo'))
    )
    # l1 <- c('a' = 4); l2 <- c('a' = 3, 'b' = 4);
    # res <- c(l2, l1); res[unique(names(res))]
    col_scale <- c(reps, col_scale) %>%
      { .[unique(names(.))] }
  }

  detected_combo_levs <-
    grep('IFNy.*TNFa', stim_group_levels, value = T)
  if (length(detected_combo_levs) > 0) {
    comb_col_scale <-
      detected_combo_levs %>%
      maartenutils::auto_name() %>%
      purrr::map_chr(function(csg) {
        ifn_groups <- stringr::str_subset(names(col_scale), 'IFNy')
        ifn_match <- recover_garbled_string(csg, ifn_groups)[[1]]
        ifn_col <- col_scale[ifn_match]
        tnf_groups <- stringr::str_subset(names(col_scale), 'TNFa')
        tnf_match <- recover_garbled_string(csg, tnf_groups)[[1]]
        tnf_col <- col_scale[tnf_match]
        circlize::colorRamp2(colors = c(ifn_col, tnf_col), breaks = c(0, 1))(.5)
      })
    col_scale <- c(col_scale, comb_col_scale)
  }

  if (F) {
    col_scale['In vivo 10 ng/ml TNFa'] <-
      maartenutils::darken(col_scale['10 ng/ml TNFa'])
    col_scale['In vivo 100 ng/ml IFNy'] <-
      maartenutils::darken(col_scale['100 ng/ml IFNy'])
    col_scale['In vivo 100 ng/ml IFNy 10 ng/ml TNFa'] <-
     maartenutils::darken(col_scale['100 ng/ml IFNy 10 ng/ml TNFa'])
  }
  col_scale['Ag-GAS + T'] <- 'grey20'
  col_scale['Mix + T'] <- col_scale['Exposed to T-cells in vivo']
  col_scale['Mix PBS'] <- 'grey60'

  if (mirjam_mode) {
    col_scale['100 ng/ml IFNy 10 ng/ml TNFa'] <- '#802A8A'
  }
  # col_scale['Ag-GAS + T'] <- 'green'
  # col_scale['Mix + T'] <- 'blue'
  # col_scale['Mix PBS'] <- 'yellow'

  col_scale['5000 U/ml IFNA1'] <- 'cyan'
  col_scale['10 ng/ml TGFb'] <- 'lightgoldenrod'
  col_scale['2 ng/ml lymphotoxin A1 B2'] <- 'chartreuse'
  col_scale['6000 U/ml IL-2'] <- 'deeppink2'

  if (!is.null(subset_vals)) {
    subset_vals <- naturalsort::naturalsort(unique(subset_vals))
    subset_vals <- as.character(unique(subset_vals))
    missing_vals <- setdiff(subset_vals, names(col_scale)) %>%
      setdiff('')
    if (!all(is.na(missing_vals)) && length(missing_vals) > 0) {
      message('Missing from col_scale: ',
        paste(missing_vals, collapse = ', '))
    }
    col_scale <- col_scale %>% { .[intersect(names(.), subset_vals)] }
  }

  col_scale <- maartenutils::as.color_vector(col_scale)

  if (F) {
    plot.color_vector <- function(cv) {
      old_par <- par()
      par(mar = rep(0, 4), plt = c(0, 1, 0, 1), oma = c(0, 0, 0, 0))
      on.exit(par(old_par))
      plot(NA, xlim = c(0, 1), ylim = c(0, length(cv)),
          axes = F, xlab = "", ylab = "")
      for (i in 1:length(cv)) {
        polygon(y = c(i - 1, i, i, i - 1), x = c(0, 0, 1, 1),
          col = cv[i])
        if (!is.null(names(cv))) {
          text(x = 0.5, y = i - 0.5, labels = names(cv)[i])
        }
      }
      invisible()
    }
    print_plot_eval({ plot.color_vector(col_scale) },
      width = 9, height = 10,
      filename = file.path(Sys.getenv('img_dir'),
        glue::glue('test_col_scale.pdf')))
  }

  return(col_scale)
}



defensive_stim_col_scale <- function(v) {
  # col <- gen_cyto_inf_col_scale(sa$stim_group)
  elevs <- e_levels(v)

  alt_col <- gen_cyto_inf_col_scale(v)

  if (F) {
    static_col <- tar_read(stim_cols)[elevs]
    static_col <- static_col[!is.na(static_col)]
    # static_col <- static_col[static_col != '#BE6186FF']
    static_col <- static_col[!is.na(names(static_col))]
    # static_col <- static_col[stringr::str_subset(names(static_col), 'SN|(?<!(IFNy)TNF)', negate = T)]
    # static_col <- static_col[stringr::str_subset(names(static_col), 'SN', negate = T)]
    # stringr::str_subset(names(static_col), 'IFNy.*TNFa$')
    static_col <- static_col[stringr::str_subset(names(static_col), 'IFNy.*TNFa$')]
    # static_col <- static_col[stringr::str_subset(names(static_col), 'SN|TNF', negate = T)]

    ## Overrride defaults with those in 'static_col'
    alt_col[names(static_col)] <- static_col
  }

  if ('' %in% elevs) {
    alt_col[''] <- 'grey50'
  }

  return(alt_col)
}


default_p_modify = function(x, use_UMAP = TRUE) {
  if (x$theme$legend.position != 'none') {
    x <- x +
      theme(
        legend.position = 'bottom',
        legend.direction = 'horizontal'
      )
  }

  x <- x +
    guides(
      size = F, shape = F, text = F,
      colour = ggplot2::guide_legend(
        ncol = 2,
        override.aes = list(alpha = 1)
      ),
      fill = ggplot2::guide_legend(
        override.aes = list(alpha = 1)
      )
    )

  if (use_UMAP) {
    x <- x + labs(x = 'UMAP 1', y = 'UMAP 2')
  } else {
    x <- x + labs(x = 'scVI 1', y = 'scVI 2')
  }

  return(x)
}


plot_embedding_panel <- function(
  p_dat,
  cleanup = TRUE,
  p_modify = default_p_modify) {

  ## This code ideally shouldn't be needed as data should have been
  ## cleaned up elsewhere
  if (cleanup) {
    ## All SC digested samples are 24 h 100 ng/ml IFNy 10 ng/ml TNFa
    p_dat$in_vitro_sc_digest <-
      ifelse(p_dat$sample_origin == 'in_vitro' &
             p_dat$sample_type == 'sc' &
             p_dat$stim_group == '100 ng/ml IFNy 10 ng/ml TNFa' &
             p_dat$duration == 24,
             c('No SC digestion', 'SC digestion')[as.integer(p_dat$sc_digestion) + 1],
             'background') %>%
      factor(levels = c('background', 'No SC digestion', 'SC digestion'))
    # table(p_dat$in_vitro_sc_digest, p_dat$stim_group)
    # table(p_dat$in_vitro_sc_digest, p_dat$duration)

    exp_tally <- table(p_dat$exp)
    p_dat$exp <- factor(p_dat$exp,
      levels = names(exp_tally)[order(exp_tally)])

    p_dat$sample_origin[p_dat$exp == '6071'] <- 'in_vivo'
    p_dat$sample_origin[p_dat$exp == '6072'] <- 'in_vivo'
    p_dat$sample_origin[p_dat$exp == '6073'] <- 'in_vivo'
    p_dat$stim_group <- gsub(' 0 ng/ml TNFa', '', p_dat$stim_group)
    p_dat$stim_group <- gsub(' tnf', ' TNFa', p_dat$stim_group)

    unexpect_stim_group <- setdiff(p_dat$stim_group, stim_group_levels)
    if (length(unexpect_stim_group) != 0) browser()

    p_dat$stim_group <-
      factor(p_dat$stim_group, levels = stim_group_levels)
    # p_dat %<>% arrange(stim_group)
    p_dat$stim_group %<>% as.character
    if ('Phase' %in% colnames(p_dat)) {
      table(p_dat$Phase)
      p_dat$Phase %<>% as.factor
    }
  }

  message('Number of rows with NA stim_group: ',
    nrow(p_dat[is.na(p_dat$stim_group), ]))
  stopifnot(all(p_dat[is.na(p_dat$stim_group), 'hash.ID'] ==
      'Negative'))
  p_dat <- p_dat[!is.na(p_dat$stim_group), ]

  p_dat_bulk <- dplyr::filter(p_dat, sample_type == 'bulk')
  p_dat_sc <- dplyr::filter(p_dat, sample_type == 'sc')

  if (any(class(p_dat_bulk) == 'integer')) browser()

  plots <- list(
    plot_UMAP.data.frame(
      dtf = p_dat_bulk,
      pt_alpha = 1, pt_size = 2,
      stim_hl_group = NULL,
      colour_var = 'stim_group'),
    plot_UMAP.data.frame(
      dtf = p_dat,
      stim_hl_group = NULL,
      colour_var = 'stim_group'),
    plot_UMAP.data.frame(p_dat, stim_hl_group = 'vivo'),
    plot_UMAP.data.frame(p_dat, stim_hl_group = 'injection_vivo'),
    # plot_UMAP.data.frame(dplyr::arrange(p_dat, in_vitro_sc_digest),
    #                      colour_var = 'in_vitro_sc_digest'),
    plot_UMAP.data.frame(p_dat, colour_var = 'sample_origin'),
    plot_UMAP.data.frame(p_dat_bulk,
      pt_alpha = 1, pt_size = 2,
      stim_hl_group = 'tnf'),
    plot_UMAP.data.frame(p_dat, stim_hl_group = 'tnf'),
    plot_UMAP.data.frame(p_dat_bulk,
      pt_alpha = 1, pt_size = 2,
      stim_hl_group = 'ifn'),
    plot_UMAP.data.frame(p_dat, stim_hl_group = 'ifn'),
    plot_UMAP.data.frame(p_dat_bulk,
      pt_alpha = 1, pt_size = 2,
      stim_hl_group = 'combo'),
    plot_UMAP.data.frame(p_dat, stim_hl_group = 'combo'),
    plot_UMAP.data.frame(p_dat, stim_hl_group = 'unstim'),
    plot_UMAP.data.frame(p_dat, colour_var = 'duration'),
    plot_UMAP.data.frame(dplyr::arrange(p_dat, desc(exp)),
                         colour_var = 'exp'),
    plot_UMAP.data.frame(p_dat, colour_var = 'Phase'),
    plot_UMAP.data.frame(p_dat_sc, colour_var = 'N_UMI')
  )

  if (length(unique(p_dat$exp)) > 1) {
    p <- dplyr::arrange(p_dat, desc(exp)) %>%
      plot_UMAP.data.frame(colour_var = 'exp')
    plots %<>% purrr::prepend(list(p))
  }

  if (!is.null(p_dat$clusters)) {
    p <- plot_UMAP.data.frame(p_dat,
      shape_var = NULL, colour_var = 'clusters') +
      theme(legend.position = 'none')

    ## Stim group as function of clusters (clusters as independent
    ## var)
    col_scale <- gen_cyto_inf_col_scale(p_dat$stim_group)
    p_dat$stim_group %<>% factor(levels = names(col_scale))
    p_clust_bar <- p_dat %>%
      dplyr::group_by(clusters, stim_group) %>%
      dplyr::summarize(N = n(), .groups = 'drop') %>%
      dplyr::group_by(clusters) %>%
      dplyr::mutate(frac = N / sum(N), .groups = 'drop') %>%
      dplyr::mutate(clusters = as.factor(clusters)) %>%
      ggplot(aes_string(x = 'clusters',
                        fill = 'stim_group', y = 'frac')) +
      geom_col() +
      scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
      scale_fill_manual(
        values = gen_cyto_inf_col_scale(p_dat$stim_group)) +
      scale_x_discrete(name = 'Cluster', expand = c(0, 0)) +
        theme_cyto_inf()

    p_clust_bar_2 <- p_dat %>%
      dplyr::group_by(clusters, stim_group) %>%
      dplyr::summarize(N = n(), .groups = 'drop') %>%
      dplyr::group_by(stim_group) %>%
      dplyr::mutate(frac = N / sum(N)) %>%
      dplyr::mutate(clusters = as.factor(clusters)) %>%
      ggplot(aes_string(x = 'stim_group',
                        fill = 'clusters', y = 'frac')) +
      geom_col() +
      scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
      # scale_fill_manual(
      #   values = gen_cyto_inf_col_scale(p_dat$stim_group)) +
      scale_x_discrete(name = 'Cluster', expand = c(0, 0)) +
        theme_cyto_inf() +
      coord_flip()

    plots %<>% purrr::prepend(list(p, p_clust_bar_2, p_clust_bar))
  }

  if (T && !is.null(p_modify)) {
    plots <- plots[!sapply(plots, is.null)]
    plots <- purrr::map(plots, p_modify, use_UMAP = T)
  }

  return(invisible(plots))
}


#' Plot all potential features of interest for a Seurat object
#'
#'
plot_variable_features <- function(seurat_obj, id,
                                   redo = F,
                                   features = c('number_reads', 'percent.mt'),
                                   N_features =
                                     length(VariableFeatures(seurat_obj))) {
  if (is.null(GetAssay(seurat_obj)@var.features)) {
    seurat_obj <- FindVariableFeatures(seurat_obj,
                                       selection.method = "vst",
                                       nfeatures = N_features)
  }
  N_features <- min(N_features, length(VariableFeatures(seurat_obj)))
  o_dir <- file.path(img_dir, id)
  if (!dir.exists(o_dir)) dir.create(o_dir)
  if (N_features > 1) {
    features <- c(features, VariableFeatures(pbmc)[1:N_features])
  }
  # library(furrr)
  # plan(transparant)
  map(features, function(feat) {
    o_path <- file.path(o_dir, glue('{feat}.png'))
    if (!file.exists(o_path) || redo) {
      message(feat)
      p <- FeaturePlot(pbmc, features = feat) +
        ggtitle(glue('{feat} - {classify_gene_using_bulk(feat)}'))
      ggsave(o_path, p)
    }
  })
  return(invisible())
}


#' Overrride this function so that we can lower font size
#'
#'
plotGseaTable <- function(pathways, stats, fgseaRes, 
  gseaParam = 1, 
  fgp = gpar(fontsize = 6),
  lgp = gpar(fontsize = 8),
  tgp = gpar(fontsize = 10),
  colwidths = c(5, 3, 0.8, 1.2, 1.2), 
  analysis_title = NULL,
  col_scale = gen_col_fun(stats),
  # col_scale = NULL,
  render = TRUE) {

  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathways <- lapply(pathways, function(p) {
      unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  pathways <- pathways[sapply(pathways, length) > 0]
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    rank_plot <- ggplot() 
    if (F && !is.null(col_scale)) {
      rank_plot <- rank_plot + 
        geom_segment(aes(
          x = 1:length(statsAdj)-.5, 
          xend = 1:length(statsAdj)+.5, 
          y = -Inf, yend = Inf, size = .1,
          color = col_scale(statsAdj)
        ))
      # rank_plot <- rank_plot + geom_rect(aes(
      #     xmin = 1:length(statsAdj)-.5, xmax = 1:length(statsAdj)+.5, 
      #     ymin = -Inf, ymax = Inf,
      #     fill = col_scale(statsAdj)))
      # for (i in 1:length(statsAdj)) {
      #   rank_plot <- rank_plot + geom_rect(aes(
      #       xmin = i-.5, xmax = i+.5, 
      #       ymin = -Inf, ymax = Inf),
      #       fill = col_scale(statsAdj[i]), show.legend = F)
      # }
    }
    zero_idx <- which.min(abs(statsAdj))
    rank_plot <- rank_plot +
      geom_segment(aes(x = zero_idx, xend = zero_idx, 
          y = -max(statsAdj), 
          yend = max(statsAdj)), size = 0.5, col = 'indianred3') 
    rank_plot <- rank_plot +
      geom_segment(aes(x = p, xend = p, y = 0, yend = statsAdj[p]), size = 0.2) + 
      scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) + 
      scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
      xlab(NULL) + ylab(NULL) +
      theme(
        panel.background = element_blank(), 
        axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), axis.title = element_blank(),
        plot.margin = rep(unit(0, "null"), 4), 
        panel.spacing = rep(unit(0, "null"), 4)
      )
    list(
      textGrob(pn, gp = fgp, just = "right", x = unit(0.95, "npc")),
      rank_plot,
      textGrob(sprintf("%.2f", annotation$NES), gp = fgp),
      textGrob(sprintf("%.1e", annotation$pval), gp = fgp), 
      textGrob(sprintf("%.1e", annotation$padj), gp = fgp)
    )
  })

  rankPlot <- ggplot() + geom_blank() + 
    scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
    xlab(NULL) + ylab(NULL) + 
    theme(
      panel.background = element_blank(),
      axis.line = element_blank(), 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), 
      panel.grid = element_blank(),
      axis.title = element_blank(), 
      plot.margin = unit(c(0, 0, 0.5, 0), "npc"), 
      panel.spacing = unit(c(0, 0, 0, 0), "npc")
    )

  grobs <- c(lapply(c("Pathway", "Gene ranks", "NES", "P-val", 
        "Adj. p-val"), textGrob, gp = lgp), 
    unlist(ps, recursive = FALSE), list(nullGrob(),
      rankPlot, nullGrob(), nullGrob(), nullGrob()))
  if (!is.null(analysis_title)) {
    grobs <- c(list(nullGrob(), textGrob(analysis_title, gp = tgp), 
        nullGrob(), nullGrob(), nullGrob()), grobs)
  }
  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))

  # library(gridExtra)
  p <- gridExtra::arrangeGrob(
    grobs = grobs[grobsToDraw], 
    ncol = sum(colwidths != 0), 
    widths = colwidths[colwidths != 0]
  )

  if (render) {
    grid.draw(p)
  }
  else {
    return(p)
  }
}


gen_HT_identifiability_fn <- function(experiment, id) {
  file.path(exp_plot_dir(experiment),
    glue::glue('{experiment}_expression_variation\\
      {make_flag(id)}.pdf'))
}


#' Create a series of plots to identify whether hash tag information
#' is detectably correlated to gene expression
#'
#'
plot_hashtag_identifiability <- function(so,
  o_fn = gen_HT_identifiability_fn(so[['exp']][1], id = ''),
  raster = T) {
  library(Seurat)
  library(ggplot2)

  # so@meta.data %<>%
  #   update_meta(experiment = unlist(so[['exp']])[1])
  so <- order_condition_name(so)
  so <- order_duration(so)
  if (all(is.na(so@meta.data[['condition_name']]))) browser()
  so@meta.data$seurat_clusters <-
    as.factor(so@meta.data$seurat_clusters)

  cluster_cols <- maartenutils::gen_color_vector(
    as.character(sort(unique(so@meta.data$seurat_clusters))),
    name = 'Moonrise3')
  condition_cols <- maartenutils::gen_color_vector(
    levels(so@meta.data[['condition_name']]),
    name = 'FantasticFox1')

  p1 <- DimPlot(so, group.by = 'condition_name',
      label.size = 2, label = TRUE,
      raster = raster) +
    theme_cyto_inf() +
    # ggplot2::scale_colour_manual(values = condition_cols) +
    Seurat::NoLegend()

  p2 <- DimPlot(so, group.by = 'seurat_clusters',
      label.size = 2, label = TRUE,
      raster = raster) +
    theme_cyto_inf() +
    scale_colour_manual(values = cluster_cols) +
    guides(color = 'none') +
    scale_x_discrete(expand = c(0, 0))

  p3 <- DimPlot(so, group.by = 'dbscan_cluster',
      label.size = 2, label = TRUE,
      raster = raster) +
    theme_cyto_inf() +
    guides(color = 'none') +
    scale_x_discrete(expand = c(0, 0))

  hashtag_by_clusters <-
    so@meta.data %>%
      dplyr::group_by(seurat_clusters, condition_name) %>%
      dplyr::summarize(N = n(), .groups = 'drop') %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::mutate(frac = N / sum(N)) %>%
      dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))

  p_clust_bar <- hashtag_by_clusters %>%
    ggplot(aes_string(x = 'seurat_clusters',
                      fill = 'condition_name', y = 'frac')) +
      geom_col() +
      scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
      scale_x_discrete(name = 'Cluster', expand = c(0, 0)) +
      theme_cyto_inf() +
      guides(fill = guide_legend(name = 'Condition', nrow = 6, byrow = F))

  # my_legend <- ggpubr::get_legend(p_clust_bar)

  p5 <- FeaturePlot(so,
    label.size = 2, label = TRUE,
    features = c('N_UMI'), raster = raster) +
    theme_cyto_inf()

  p6 <- DimPlot(so,
    label.size = 2, label = TRUE,
    group.by = c('Phase'), raster = raster) + theme_cyto_inf()

  p_ribo <- FeaturePlot(so,
    label.size = 2, label = TRUE,
    features = c('percent.ribo'), raster = raster) +
    theme_cyto_inf()

  library(patchwork)
  ggsave(o_fn,
    (p1 | p2 | p3) /
    (p_clust_bar) /
    (p5 | p6 | p_ribo) /
    plot_layout(heights = c(1, 1.6, 1)),
    width = 17.4, height = 25, units = 'cm')

  pt = 0.02
  p_mt <- FeaturePlot(so,
    label = T,
    features = c('percent.mt'), raster = raster) +
    theme_cyto_inf(legend.position = 'none')
  p_violin_mt <- VlnPlot(so,
    pt.size = pt,
    features = 'percent.mt', group.by = 'condition_name') +
    coord_flip() +
    theme_cyto_inf(legend.position = 'none')
  p_violin_cl_mt <- VlnPlot(so,
    pt.size = pt,
    features = 'percent.mt', group.by = 'seurat_clusters') +
    coord_flip() +
    theme_cyto_inf(legend.position = 'none')
  mt_plots_fn <- add_flag(o_fn, '-mito')
  ggsave(mt_plots_fn,
    (p_mt / p_violin_mt / p_violin_cl_mt),
    width = 17.4, height = 25, units = 'cm')

  gene_plots <-
    FeaturePlot(so, head(VariableFeatures(so), 24),
      raster = raster) &
    theme_cyto_inf(legend.position = 'none')
  example_plots_fn <- add_flag(o_fn, '-variable_feature_plots')
  ggsave(example_plots_fn,
    gene_plots,
    width = 17.4, height = 25, units = 'cm')

  if ('mouse' %in% colnames(so@meta.data)) {
    p <- DimPlot(so,
      label.size = 2, label = TRUE,
      group.by = c('mouse'), raster = raster) + theme_cyto_inf()

    condition_by_cluster <-
      so@meta.data %>%
        dplyr::mutate(condition_name = paste(condition_name,
            ' - mouse ', mouse, sep = '')) %>%
        dplyr::group_by(condition_name, seurat_clusters) %>%
        dplyr::summarize(N = n(), .groups = 'drop') %>%
        dplyr::group_by(seurat_clusters) %>%
        dplyr::mutate(frac = N / sum(N)) %>%
        dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))

    p_clust_bar_3 <- condition_by_cluster %>%
      ggplot(aes_string(x = 'seurat_clusters',
                        fill = 'condition_name', y = 'frac')) +
        geom_col() +
        scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
        scale_fill_discrete(name = 'Condition') +
        scale_x_discrete(name = 'Cluster', expand = c(0, 0)) +
        theme_cyto_inf() +
        guides(fill = guide_legend(ncol = 2, byrow = F))

    clusters_by_mouse <-
      so@meta.data %>%
        dplyr::group_by(mouse, seurat_clusters) %>%
        dplyr::summarize(N = n(), .groups = 'drop') %>%
        dplyr::group_by(mouse) %>%
        dplyr::mutate(frac = N / sum(N)) %>%
        dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))

    p_clust_bar_4 <- clusters_by_mouse %>%
      ggplot(aes_string(x = 'mouse',
                        fill = 'seurat_clusters', y = 'frac')) +
        geom_col() +
        scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
        scale_fill_discrete(name = 'Clusters') +
        scale_x_discrete(name = 'mouse', expand = c(0, 0)) +
        theme_cyto_inf() +
        guides(fill = guide_legend(ncol = 2, byrow = F))

    mouse_by_condition <-
      so@meta.data %>%
        dplyr::group_by(mouse, condition_name) %>%
        dplyr::summarize(N = n(), .groups = 'drop') %>%
        dplyr::group_by(mouse) %>%
        dplyr::mutate(frac = N / sum(N))

    p_clust_bar_5 <- mouse_by_condition %>%
      ggplot(aes_string(x = 'condition_name',
                        fill = 'mouse', y = 'frac')) +
        geom_col() +
        scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
        scale_fill_discrete(name = 'Clusters') +
        scale_x_discrete(name = 'mouse', expand = c(0, 0)) +
        theme_cyto_inf() +
        guides(fill = guide_legend(ncol = 2, byrow = F))

    mouse_plots_fn <- add_flag(o_fn, '-mouse')
    ggsave(mouse_plots_fn, p / p_clust_bar_3 /
      (p_clust_bar_4 + p_clust_bar_5),
      width = 17.4, height = 25, units = 'cm')
  }

  condition_frac <-
    so@meta.data %>%
      dplyr::group_by(condition_name) %>%
      dplyr::summarize(N = n(), .groups = 'drop') %>%
      dplyr::mutate(frac = N / sum(N))

  condition_pie <- condition_frac %>%
    ggplot(aes(x = '', fill = condition_name, y = N)) +
      geom_bar(stat = 'identity', width = 1) +
      coord_polar(theta = 'y', start = 0) +
      theme_tabula_rasa +
      guides(guide_legend(nrow = 6, byrow = F))

  ## Make separate barplots, repeat the first one
  cluster_by_condition <-
    so@meta.data %>%
      dplyr::group_by(seurat_clusters, condition_name) %>%
      dplyr::summarize(N = n(), .groups = 'drop') %>%
      dplyr::group_by(condition_name) %>%
      dplyr::mutate(frac = N / sum(N)) %>%
      dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))

  p_condition_bar <- cluster_by_condition %>%
    ggplot(aes_string(x = 'condition_name',
                      fill = 'seurat_clusters', y = 'frac')) +
      geom_col() +
      scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
      # scale_fill_manual(name = 'Condition', values = condition_cols) +
      scale_fill_discrete(name = 'Cluster') +
      rotate_x_labels(90) +
      scale_x_discrete(name = 'Condition', expand = c(0, 0)) +
      theme_cyto_inf() +
      guides(fill = guide_legend(nrow = 6, byrow = F))

  bar_fn <- add_flag(o_fn, '-prop_bar')
  ggsave(bar_fn, p_clust_bar / (condition_pie + p_condition_bar),
    width = 17.4, height = 20, units = 'cm')

  return(o_fn)
}


time_informativeness_panel <- function(gene, o_fn = NULL) {
  p1 <- VlnPlot(so_6489, features = gene,
    group.by = 'condition_name') &
    theme_cyto_inf(legend.position = 'none') &
    coord_flip()

  # p1_a <- VlnPlot(so_6489_m, features = gene,
  #   group.by = 'condition_name') &
  #   theme_cyto_inf(legend.position = 'none') &
  #   coord_flip()

  p2 <- VlnPlot(so_5310, features = gene,
    group.by = 'condition_name') &
    theme_cyto_inf(legend.position = 'none') &
    coord_flip()

  p_title <- selection_crit_table_u %>%
    dplyr::filter(gene == .env[['gene']]) %>%
    dplyr::mutate(
      p_title = glue('{gene} - {time_class} - {max_timepoint}')) %>%
    pull(p_title)

  p3 <- plot_expression_dynamics(fn = gene, redo = F) +
    ggtitle(p_title)

  comb <- p1 / (p2 + p3) + plot_layout(ncol = 1)

  print_plot(comb, fn = o_fn, w = 17.4, h = 20)
}


plot_lib_size <- function(p_dat, fn = 'TMM library size') {
  p_dat$duration <- factor(p_dat$duration,
    levels = naturalsort::naturalsort(unique(p_dat$duration)))
  p_dat$stim_group <- p_dat$stim_group %>%
    factor(levels = stim_group_levels)
  stopifnot(!any(is.na(p_dat$stim_group)))
  p_dat <- dplyr::select(p_dat, gexp, col, duration, stim_group)
  p <- plot_expression_dynamics(p_dat = p_dat,
    caching = F, fn = fn, axis_name = fn)
  return(p)
}


plot_perc_zeroes <- function(so_list) {
  sgn <- rlang::exec(find_shared_genes, !!!so_list)
  Ms <- map(so_list, ~so2M(.x)[sgn, ])
  p_dat <- tibble(
    class = flatten_chr(imap(Ms, ~rep(.y, ncol(.x)))),
    perc_nonzero = map(
      Ms, ~apply(.x, 2, function(x) mean(x > 0))) %>%
      flatten_dbl
    )
  ggplot(p_dat,
    aes(x = perc_nonzero, fill = class)) +
    geom_density(alpha = .5) +
    xlab('Fraction of genes non-zero')
}


annotation_par <-
  gpar(fontsize = heatmap_fs, col = 'black',
    lwd = 0, fontfamily = heatmap_font)


legend_params <- list(
  title_gp = gpar(fontsize = heatmap_fs, fontface = 'italic'),
  direction = 'horizontal',
  labels_gp = gpar(fontsize = heatmap_fs),
  grid_height = unit(3, 'mm'),
  grid_width = unit(3, 'mm')
)


if (F) {
  duration_palette <-
    # col <- maartenutils::gen_color_vector(e_levels(sa$duration),
    #   name = maartenutils::col_to_hex('slateblue3'))
  {
    N_cols <- 7
    cols <- viridis::mako(
      N_cols, direction = -1L)[1:(N_cols-2)]
    rlang::set_names(cols, c(0, 2, 6, 12, 24))
  } %>%
    c(c('Unknown' = 'white'))
} else if (F) {
  ## Ugly yellow colors
  duration_palette <-
    setNames(
      c('#FAFCC9', '#FCD929', '#9B8800', '#544E1F',
        '#2E2700', '#1F1A00'),
      c(2, 6, 12, 24, 48, 96)) %>%
    c(c('Unknown' = 'white'))
} else if (T) {
  duration_palette <-
    circlize::colorRamp2(
      # breaks = seq(0, 1, length.out = 6),
      breaks = c(0, 1),
      colors = c('grey90', 'grey10')
    )(seq(0, 1, length.out = 6)) %>%
    setNames(c(2, 6, 12, 24, 48, 96)) %>%
    c(c('Unknown' = 'white')) %>%
    c(c('0' = 'white'))
}

# score_palette <- viridis::mako(3L, direction = -1L)
# score_palette <- pryr::partial(
#   `_f` = viridis::mako,
#   begin = 0, end = 1, direction = -1L)
# score_palette(.7)
# score_palette(seq(0, 1, by = .2))
# viridis::mako(seq(0, 1, by = .2))
# map(seq(0, 1, by = .2), ~viridis::mako(10)(.x))

HM_col_funs <- list(
  'sample_origin' = function(sa) {
     col <- maartenutils::gen_color_vector(
       e_levels(sa$sample_origin, F))
     col
  },
  'name' = function(sa) {
     col <- maartenutils::gen_color_vector(
       e_levels(sa$name, F))
     col
  },
  'exp' = function(sa, cn='exp') {
     col <- maartenutils::gen_color_vector(
       e_levels(sa[[cn]], F))
     col
  },
  'sample_type' = function(sa) {
    col <- maartenutils::gen_color_vector(
      e_levels(sa$sample_type, F),
      name = 'GrandBudapest1')
  },
  'pred_score' = function(sa) {
    col <- maartenutils::gen_color_vector(
      e_levels(sa$pred_score, F),
      name = 'GrandBudapest1')
    return(col)
  },
  'mouse' = function(sa) {
    col <- maartenutils::gen_color_vector(
      e_levels(sa$mouse, F),
      name = 'GrandBudapest1')
    return(col)
  },
  'rep_N' = function(sa) {
    col <- maartenutils::gen_color_vector(
      e_levels(sa$rep_N, F),
      name = 'GrandBudapest1')
    return(col)
  },
  'sn_dilution' = function(sa)
    defensive_stim_col_scale(sa$sn_dilution),
  # 'tnf_conc' = function(sa)
  #   defensive_stim_col_scale(sa$tnf_conc),
  'tnf_conc' = function(sa) {
    elevs <-
      as.character(e_levels(sa$tnf_conc)) %>%
      maartenutils::auto_name() %>%
      purrr::map_chr(function(.x) {
        if (.x != 0) {
          return(glue::glue('{.x} ng/ml TNFa'))
        } else {
          return('')
        }
      })
    out <- defensive_stim_col_scale(elevs)
    names(out) <- names(elevs)[match(names(out), elevs)]
    return(out)
  },
  'ifn_conc' = function(sa) {
    elevs <-
      as.character(e_levels(sa$ifn_conc)) %>%
      maartenutils::auto_name() %>%
      purrr::map_chr(function(.x) {
        if (.x != 0) {
          return(glue::glue('{.x} ng/ml IFNy'))
        } else {
          return('')
        }
      })
    out <- defensive_stim_col_scale(elevs)
    names(out) <- names(elevs)[match(names(out), elevs)]
    return(out)
  },
  'stim_group' = function(sa, cn='stim_group') {
    defensive_stim_col_scale(sa[[cn]])
  },
  'm_sg' = function(sa, cn='m_sg') {
    human_cols <- defensive_stim_col_scale(sa[[cn]])
    new_cols <- setdiff(levels(sa[[cn]]), names(human_cols))
    if (F && length(new_cols) > 0 && all(new_cols == '0.1 ng/ml IFNy')) {
      message('Applying ugly fix')
      
      f <- circlize::colorRamp2(c(1, 0), 
        c(human_cols['100 ng/ml IFNy'], 
          human_cols['0.01 ng/ml IFNy']))
      
      out <- c(human_cols, setNames(f(.5), '0.1 ng/ml IFNy'))
    } else {
      out <- c(human_cols,
        gen_color_vector(new_cols, name = 'FantasticFox1'))
    }
    return(out)
  },
  'partner_cor_nrank' = function(sa) {
    col <- tryCatch(circlize::colorRamp2(
      # colors = c('white', viridis::mako(3)[3L]),
      colors = c(duration_palette['2'], duration_palette['24']),
      breaks = c(.8, 1)
    ), error = function(e) { browser() })
    return(col)
  },
  'duration' = function(sa, cn = 'duration') {
    # levs <- e_levels(sa[[cn]])
    # obs_range <- range(as.numeric(levs), na.rm = T)
    # if (all(obs_range == c(0, 1)))
    #   col <- score_palette(levs)
    # else
    #   col <- duration_palette[levs]
    breaks <- unique(range(sa[[cn]]))
    if (any(!is.finite(breaks)) || length(breaks) <= 1L) {
      return(duration_palette[as.character(breaks)])
    }
    if (max(breaks) == 0) breaks <- c(0, 1)
    col <- tryCatch(circlize::colorRamp2(
      # colors = c('white', viridis::mako(3)[3L]),
      colors = c(duration_palette['2'], duration_palette['24']),
      breaks = breaks
    ), error = function(e) { browser() })
    return(col)
  },
  'duration_discrete' = function(sa, cn = 'duration') {
    levs <- e_levels(sa[[cn]])
    if (any(!levs %in% names(duration_palette))) {
      N_cols <- length(levs)
      col <- viridis::mako(
        N_cols, direction = 1L)[1:(N_cols)] %>%
        rlang::set_names(levs)
    } else {
      col <- duration_palette[levs]
    }
    return(col)
  },
  # 'tnf_duration' = function(sa) {
  #   col <- duration_palette[e_levels(sa$tnf_duration)]
  #   return(col)
  # },
  # 'ifn_duration' = function(sa) {
  #   col <- duration_palette[e_levels(sa$ifn_duration)]
  #   return(col)
  # },
  # 'sn_duration' = function(sa) {
  #   col <- duration_palette[e_levels(sa$sn_duration)]
  #   return(col)
  # },
  'var_ref' = function(sa)
    circlize::colorRamp2(
      colors = c('white', 'salmon3'),
      breaks = c(0, max(sa$var_ref, sa$var_SC))
    )
  ,
  'var_SC' = function(sa)
    circlize::colorRamp2(
      colors = c('white', 'salmon3'),
      breaks = c(0, max(sa$var_ref, sa$var_SC))
    )
  ,
  'geneset' = function(sa) {
    levs <- e_levels(sa$geneset, F)
    if (length(union(levs, c('TNFa', 'IFNy', 'synergy'))) == 3) {
      query <- c('10 ng/ml TNFa', '100 ng/ml IFNy',
            '100 ng/ml IFNy 10 ng/ml TNFa')
      cols <- set_names(defensive_stim_col_scale(query)[query],
        c('TNFa', 'IFNy', 'synergy'))
    } else {
      cols <- maartenutils::gen_color_vector(levs, 'Spectral')
    }
    cols[levs]
  },
  'experiment' = function(sa, cn='experiment')
    maartenutils::gen_color_vector(
      e_levels(sa[[cn]], F), 'GrandBudapest1'),
  # 'synergy_score' = circlize::colorRamp2(
  #   colors = c('royalblue3', 'grey95', 'salmon3'),
  #   breaks = max(abs(ra$synergy_score)) %>% { c(-., 0, .) }
  # ),
  'synergy_score' = function(sa) gen_col_fun(sa$synergy_score),
  'TNFa_bias' = function(sa) circlize::colorRamp2(
    colors = rev(c('royalblue3', 'white', 'salmon3')),
    breaks = c(0, .5, 1)
  ),
  'time_class' = function(sa) maartenutils::gen_color_vector(
    e_levels(sa$time_class),
      name = 'GrandBudapest1'),
  'max_timepoint' = function(sa) maartenutils::gen_color_vector(
    as.character(e_levels(sa$max_timepoint)),
    name = 'GrandBudapest1'),
  'frozen' = function(sa) maartenutils::gen_color_vector(
    as.character(e_levels(sa$rozen)),
    name = 'GrandBudapest1')
)
# ## Append copies with '2' appended to the name
HM_col_funs <- append(HM_col_funs,
  setNames(HM_col_funs, paste0(names(HM_col_funs), '2')))

# stopifnot(sapply(HM_col_funs, is.function))

print_names <- c(
  'sample_origin' = 'Origin',
  'sample_type' = 'Sample type',
  'exp' = 'Experiment',
  'experiment' = 'Experiment',
  'pred_score' = 'Predicted score',
  'ifn_conc' = '[IFNy]',
  'tnf_conc' = '[TNFa]',
  'sn_dilution' = 'SN dilution',
  'pred_ifn_conc' = 'Pred. [IFNy]',
  'pred_tnf_conc' = 'Pred. [TNFa]',
  'pred_sn_dilution' = 'Pred. SN dilution',
  'stim_group' = 'Exposure type',
  'stim_group2' = 'Exposure type2',
  'm_sg' = 'Exposure type',
  'm_sg2' = 'Exposure type2',
  'sg' = 'Exposure type',
  'duration' = 'Duration',
  'duration2' = 'Duration2',
  'partner_cor_nrank' = 'Rank norm. cor.\nwith rep.',
  'tnf_duration' = 'TNFa duration',
  'ifn_duration' = 'IFNy duration',
  'sn_duration' = 'SN duration',
  'pred_tnf_duration' = 'Pred. TNFa duration',
  'pred_ifn_duration' = 'Pred. IFNy duration',
  'pred_sn_duration' = 'Pred. SN duration',
  'var_ref' = 'Reference exp. var.',
  'var_SC' = 'Query exp. var.',
  'geneset' = 'Gene label',
  'frozen' = 'Freezing step',
  'synergy_score' = 'Synergy score',
  'time_class' = 'Time informativeness',
  'max_timepoint' = 'Max expression timepoint',
  'tnf_conc_bin' = 'Any TNFa',
  'ifn_conc_bin' = 'Any IFNy'
)

print_names_nc <- c(
  'sample_origin' = 'origin',
  'sample_type' = 'sample type',
  'exp' = 'experiment',
  'experiment' = 'experiment',
  'pred_score' = 'predicted score',
  'ifn_conc' = '[IFNy]',
  'tnf_conc' = '[TNFa]',
  'sn_dilution' = 'SN dilution',
  'pred_ifn_conc' = 'pred. [IFNy]',
  'pred_tnf_conc' = 'pred. [TNFa]',
  'pred_sn_dilution' = 'pred. SN dilution',
  'stim_group' = 'exposure type',
  'm_sg' = 'exposure type',
  'sg' = 'exposure type',
  'duration' = 'duration',
  'frozen' = 'Freezing step',
  'partner_cor_nrank' = 'rank norm. cor.\nwith rep.',
  'tnf_duration' = 'TNFa duration',
  'ifn_duration' = 'IFNy duration',
  'sn_duration' = 'SN duration',
  'pred_tnf_duration' = 'pred. TNFa duration',
  'pred_ifn_duration' = 'pred. IFNy duration',
  'pred_sn_duration' = 'pred. SN duration',
  'var_ref' = 'reference exp. var.',
  'var_SC' = 'query exp. var.',
  'geneset' = 'gene label',
  'synergy_score' = 'synergy score',
  'time_class' = 'time informativeness',
  'max_timepoint' = 'max expression timepoint'
)


change_colnames <- function(orig) {
  if (is.null(orig)) return(NULL)
  if (length(orig) == 0) return(orig)

  purrr::map_chr(orig, function(cn) {
    if (cn %in% names(print_names))
      return(print_names[cn])
    else
      return(cn)
  })
}


get_HM_cols <- function(sa_sel) {
  cols <-
    auto_name(colnames(sa_sel)) %>%
    purrr::map(function(.x) {
      # if (.x == 'm_sg2') browser()
      fun <- HM_col_funs[[.x]]

      if (is.null(fun)) {
        rlang::warn(glue::glue('No color scheme found for {.x}'))
        lcols <- maartenutils::gen_color_vector(
            arg = as.character(unique(sa_sel[[.x]])),
            name = 'FantasticFox1')
        return(lcols)
      }

      # if (stringr::str_detect(.x, 'duration|conc|dilution')) {
      if (stringr::str_detect(.x, 'duration')) {
        if (is.factor(sa_sel[[.x]]) || is.integer(sa_sel[[.x]])) {
          sa_sel[[.x]] <- order_factor(sa_sel[[.x]])
          lcols <- HM_col_funs[['duration_discrete']](sa_sel, .x)
        } else if (is.numeric(sa_sel[[.x]])) {
          lcols <- HM_col_funs[['duration']](sa_sel, .x)
        }
        return(lcols)
      }

      if (!is.null(fun)) {
        lcols <- fun(sa_sel, cn=.x)
        if (any(is.na(lcols))) {
          rlang::warn('No mapping for: ', lcols[is.na(lcols)])
          lcols <- maartenutils::gen_color_vector(
              unique(sa_sel[[.x]]), 'FantasticFox1')
        } else  {
          return(lcols)
        }
      } else {
        return(NULL)
      }
    }) %>%
    purrr::discard(is.null)

  return(cols)
}


gen_sample_annotation_HM <- function(
  sa,
  annotation_name_side = switch(heatmap_type,
    'column' = 'right', 'row' = 'top'),
  ann_ordering = NULL,
  drop_non_discriminatory = TRUE,
  show_legend = TRUE,
  height =  unit(.33, 'cm'),
  # col_regex = 'score|pred',
  col_regex = '.*',
  heatmap_type = 'column',
  annotation_name_rot = NULL,
  ...) {

  if (maartenutils::null_dat(sa)) {
    return(NULL)
  }

  # if (heatmap_type == 'row') browser()

  sa <- tibble::as_tibble(sa)
  if (length(show_legend) > 1)
    names(show_legend) <- colnames(sa)

  if (F) {
    if ('tnf_conc' %in% colnames(sa))
      sa$tnf_conc <- ifelse(!is.na(sa$tnf_conc) & sa$tnf_conc > 0,
        glue::glue('{sa$tnf_conc} ng/ml TNFa'), NA_character_)
    if ('ifn_conc' %in% colnames(sa))
      sa$ifn_conc <- ifelse(!is.na(sa$ifn_conc) & sa$ifn_conc > 0,
        glue::glue('{sa$ifn_conc} ng/ml IFNy'), NA_character_)
    if ('sn_dilution' %in% colnames(sa))
      sa$sn_dilution <- ifelse(!is.na(sa$sn_dilution) &
        sa$sn_dilution > 0,
        glue::glue('1/{sa$sn_dilution} SN'), NA_character_)
    if ('duration' %in% colnames(sa))
      sa$duration <- factor(sa$duration,
        levels = naturalsort::naturalsort(unique(sa$duration)))
    if ('tnf_duration' %in% colnames(sa))
      sa$tnf_duration <- factor(sa$tnf_duration,
        levels = naturalsort::naturalsort(unique(sa$tnf_duration)))
    if ('ifn_duration' %in% colnames(sa))
      sa$ifn_duration <- factor(sa$ifn_duration,
        levels = naturalsort::naturalsort(unique(sa$ifn_duration)))
    if ('sn_duration' %in% colnames(sa))
      sa$sn_duration <- factor(sa$sn_duration,
        levels = naturalsort::naturalsort(unique(sa$sn_duration)))
    if ('pred_score' %in% colnames(sa)) {
      sa$pred_score <- round(sa$pred_score, 2)
      # sa$pred_score <- factor(sa$pred_score,
      #   levels = sort(unique(sa$pred_score)))
    }
  }

  if (F) {
    for (cn in setdiff(names(alp_c), colnames(sa))) {
      sa[[cn]] <- factor(1)
    }
  }


  sa_sel <- dplyr::select(sa,
    any_of(names(print_names)), matches(col_regex))

  ## Drop non-discriminatory annotations
  if (drop_non_discriminatory) {
    sal_sel <- perform_drop_non_discriminatory(sa_sel)
  }

  if (!is.null(ann_ordering))
    sa_sel <- sa_sel[, intersect(ann_ordering, colnames(sa_sel))]

  if (T) {
    cols <- get_HM_cols(sa_sel)
    # names(cols)
    # map(cols, class)
  } else {
    cols = NULL
  }

  colnames(sa_sel) <- change_colnames(orig = colnames(sa_sel))
  names(cols) <- change_colnames(orig = names(cols))
  names(show_legend) <- change_colnames(orig = names(show_legend))
  if (length(show_legend) > 1)
    show_legend <- show_legend[colnames(sa_sel)]

  # RA2 <- HeatmapAnnotation(
  #   `Mean expression` = anno_barplot(unlist(gene_stats_M[, 'mean']),
  #     height = unit(2, 'cm')),
  #   # `Max expression` = anno_barplot(unlist(gene_stats_M[, 'max']),
  #   #   height = unit(2, 'cm')),
  #   annotation_name_side = 'right',
  #   which = 'column')
  # print_plot_eval(
  #   {
  #     draw(RA2)
  #   },
  #   width = 17.4, height = 15,
  #   filename = file.path(Sys.getenv('img_dir'),
  #     glue::glue('test.png')))
  if (T || 'stim_group' %in% colnames(sa_sel)) {
    alp <- list(
      stim_group = list(ncol = 2),
      `Exposure type` = list(ncol = 1)
    )
  } else {
    alp <- NULL
  }
  alp <- alp[intersect(names(alp), colnames(sa_sel))]

  ann <- HeatmapAnnotation(
    df = as.data.frame(sa_sel),
    gap = unit(0, 'mm'),
    annotation_name_side = annotation_name_side,
    annotation_name_rot = annotation_name_rot,
    annotation_legend_param = alp,
    gp = annotation_par,
    show_annotation_name = T,
    annotation_name_gp = annotation_par,
    show_legend = show_legend,
    simple_anno_size = height,
    which = heatmap_type,
    col = cols
  )

  return(ann)
}


gen_gene_labels_HM <- function(ra, heatmap_type = 'row') {
  if (maartenutils::null_dat(ra)) {
    return(NULL)
  }

  ra <- as_tibble(ra)

  ra <- ra[, scn, drop = F]
  ra_sel <-
    ra[, purrr::map_lgl(ra, ~data.table::uniqueN(.x) > 1), drop = F]

  ## HeatmapAnnotation doesn't like tibblee
  row_labels <- HeatmapAnnotation(
    df = as.data.frame(ra_sel),
    gap = unit(0, 'mm'),
    gp = annotation_par,
    show_annotation_name = T,
    annotation_name_gp = annotation_par,
    col = alc_r[intersect(names(alc_r), colnames(ra_sel))],
    # annotation_legend_param = alp_r[colnames(ra_sel)],
    which = heatmap_type
  )

  return(row_labels)
}


GE_ecdf_plot <- function(
  object_list, assay = NULL,
  q = .99,
  max_samples = 50L,
  pick_largest_samples = TRUE,
  tx = function(x) log10(x + 1)) {

  p_dat <- purrr::imap_dfr(
    object_list,
    function(so, name) {
      if (is.null(assay) || !assay %in% names(so@assays))
        assay <- DefaultAssay(so)
      M <- so2M(so, assay = assay)
      if (!is.null(max_samples) & max_samples < ncol(M)) {
        nM <- ncol(M)
        if (!pick_largest_samples) {
          M <- M[, base::sample(1:nM, min(max_samples, nM))]
        } else {
          idxs <- which(order(colSums(M)) %in% 1:max_samples)
          M <- M[, idxs]
        }
      }
      M2tibble(M) %>%
        tidyr::pivot_longer(cols = where(is.numeric)) %>%
        dplyr::mutate(type = .env[['name']]) %>%
        dplyr::mutate(experiment = first_non_NA(so@meta.data$exp))
    }
  )

  point_summary_dat <- p_dat %>%
    dplyr::group_by(type, name) %>%
    dplyr::summarize(x = q, y = quantile(value, q))

  mean_point <- p_dat %>%
    dplyr::group_by(name) %>%
    dplyr::summarize(
      type, x = mean(value), y = mean(x >= value)[1]) %>%
    dplyr::distinct(.keep_all = T)

  p <- ggplot(p_dat,
    aes(x = tx(value), colour = type, group = name)) +
    stat_ecdf(alpha = .5) +
    geom_point(data = mean_point,
      aes(y = y, x = tx(x)),
      alpha = .5, shape = 4) +
    geom_point(data = point_summary_dat,
      aes(y = x, x = tx(y)),
      alpha = .5, shape = 2) +
    ylab('Cumulative fraction of genes') +
    xlab('Gene expression [log10(x = 1)]') +
    theme(legend.text = element_text(size = 6))
  return(p)
}


#' Unfinished
#'
#'
gene_batch_correlation_plot <- function(object_list, q = .99) {
  p_dat <- imap_dfr(object_list,
    function(so, type) {
      so2M(so) %>%
        M2tibble() %>%
        pivot_longer(cols = where(is.numeric)) %>%
        dplyr::mutate(type = type)
    })

  summary_dat <- p_dat %>%
    group_by(type, name) %>%
    summarize(x = q, y = quantile(value, q))

  p <- ggplot(p_dat,
    aes(x = log10(value + 1), colour = type, group = name)) +
    stat_ecdf(alpha = .5) +
    geom_point(data = summary_dat, aes(y = x, x = log10(y + 1)),
      alpha = .5) +
    ylab('Cumulative fraction of genes') +
    xlab('Gene expression [log10(x = 1)]') +
    theme(legend.text = element_text(size = 6))
}


theme_tabula_rasa <- gg_tabula_rasa <-
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
  )


gen_sc_cluster_comp_plot <- function(sc_so) {
  cluster_cols <- maartenutils::gen_color_vector(
    as.character(sort(unique(sc_so@meta.data$seurat_clusters))),
    name = 'Moonrise3')

  stimulus_cols <- gen_cyto_inf_col_scale()
  condition_names <- unique(sc_so@meta.data[['condition_name']])
  condition_name_matches <-
    purrr::map_chr(condition_names,
    function(q) {
      recover_garbled_string(q, names(stimulus_cols))
    })
  condition_cols <- stimulus_cols[condition_name_matches] %>%
    rlang::set_names(condition_names)
  condition_order <- condition_name_matches %>%
    { .[order(match(names(stimulus_cols), .))] } %>%
    setdiff(NA)

  # condition_cols <- maartenutils::gen_color_vector(
  #   unique(sc_so@meta.data[['condition_name']]),
  #   name = 'FantasticFox1')

  condition_by_clusters <-
    sc_so@meta.data %>%
      dplyr::group_by(seurat_clusters, condition_name) %>%
      dplyr::summarize(N = n(), .groups = 'drop') %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::mutate(frac = N / sum(N)) %>%
      dplyr::mutate(seurat_clusters = as.factor(seurat_clusters)) %>%
      dplyr::mutate(condition_name = factor(condition_name,
          levels = condition_name))

  clust_obj <- tidyr::pivot_wider(condition_by_clusters,
    id_cols = condition_name,
    names_from = seurat_clusters,
    values_from = frac, values_fill = 0) %>%
    dplyr::select(where(is.numeric)) %>%
    as.matrix() %>%
    t() %>%
    dist() %>%
    hclust()

  condition_by_clusters <-
    condition_by_clusters %>%
    dplyr::mutate(seurat_clusters = factor(seurat_clusters,
        levels = clust_obj$order-1))

  p_clust_bar <- condition_by_clusters %>%
    ggplot(aes_string(x = 'seurat_clusters',
                      fill = 'condition_name', y = 'frac')) +
      geom_col() +
      scale_y_continuous(name = 'Fraction', expand = c(0, 0)) +
      scale_fill_manual(name = 'Single cell condition',
        values = condition_cols) +
      scale_x_discrete(name = 'Cluster', expand = c(0, 0)) +
      # theme_cyto_inf() +
      gg_tabula_rasa +
      theme(panel.border = theme_cyto_inf()$panel.border) +
      guides(fill =
        guide_legend(
          # nrow = min(6, floor(length(condition_cols) / 1.5)),
          ncol = 1L,
          byrow = F))

  return(p_clust_bar)
}


add_umap <- function(
  dtf,
  column_selector = where(is.numeric),
  add_source_type = F,
  umap_obj = NULL,
  verbose = F) {

  if (maartenutils::null_dat(dtf)) return(NULL)

  library(rlang)
  library(umap)
  before_atts <- attributes(dtf) %>%
    { .[!names(.) %in% c('names', 'row.names', 'class')] }
  dtf <- dplyr::select(dtf, -matches('UMAP|umap'))

  source_data <- dtf %>%
    dplyr::select(!! enquo(column_selector))

  if (verbose) {
    message(paste0(colnames(source_data), collapse = ', '))
  }

  if (all(is.na(source_data))) {
    coords <- tibble(
      'UMAP1' = rep(NA_real_, nrow(source_data)),
      'UMAP2' = rep(NA_real_, nrow(source_data))
    )
  } else {
    if (is.null(umap_obj)) {
      umap_obj <- umap::umap(source_data)
      coords <- umap_obj$layout
    } else {
      coords <- predict(umap_obj, source_data)
    }

    coords <- coords %>%
      magrittr::set_colnames(c('UMAP1', 'UMAP2')) %>%
      as_tibble()
    # coords <- umap_obj %>%
    #   purrr::pluck('layout') %>%
    #   magrittr::set_colnames(c('UMAP1', 'UMAP2')) %>%
    #   as_tibble()
  }

  out <- dplyr::bind_cols(coords, dtf)
  attributes(out) <- c(attributes(out), before_atts)
  attr(out, 'umap_obj') <- umap_obj

  if (add_source_type) {
    source_types <- c(
      'precise' = 'precise',
      'scVI' = 'scvi|scVI',
      'Harmony' = 'harmony|HARMONY',
      'PCA' = 'PC_*\\d+',
      'ICA' = 'IC\\d+'
    )
    for (i in seq_along(source_types)) {
      if (any(stringr::str_detect(colnames(dtf), source_types[i]))) {
        if (exists('source_type')) {
          cat(source_type, names(source_types)[i], '\n')
          stop('Two umap types detected')
        }
        source_type <- names(source_types)[i]
      }
    }
    stopifnot(exists('source_type'))
    attr(out, 'umap_type') <- source_type
    attr(out, 'source_type') <- source_type
  }

  return(out)
}


get_stim_group_cols <- function(...) UseMethod('get_stim_group_cols')


get_stim_group_cols.character <- function(v) {
  # cols <- tar_read(stim_cols)[e_levels(order_stim_group(v))]
  levs <- order_stim_group(v) %>% sort() %>% as.character()
  cols <- 
    defensive_stim_col_scale(as.character(levs)) %>%
    { .[match(levs, names(.))] }
  missed_levels <- setdiff(v, c(names(cols), 'Unknown'))
  if (length(missed_levels) > 0) {
    rlang::warn('Unexpected stim col requested: ',
      paste(missed_levels, collapse = ', '))
    cols <- c(cols, gen_color_vector(missed_levels, 'FantasticFox1'))
  }
  return(cols)
}


get_stim_group_cols.factor <- function(v) {
  get_stim_group_cols(e_levels(v))
}


get_stim_group_cols.data.frame <-
get_stim_group_cols.tbl_df <- function(dtf) {
  get_stim_group_cols(dtf$stim_group)
}


scale_colour_stim_group <- function(...)
  UseMethod('scale_colour_stim_group')


scale_colour_stim_group.data.frame <- function(x, alpha = .3,
  UMAP_mode = F) {
  scale_colour_stim_group(x$stim_group, alpha = alpha, 
    UMAP_mode = UMAP_mode)
}


scale_colour_stim_group.factor <- function(x, alpha = .3, 
  UMAP_mode = UMAP_mode) {
  scale_colour_stim_group(e_levels(x), alpha = alpha, 
    UMAP_mode = UMAP_mode)
}


scale_colour_stim_group.character <-
scale_colour_stim_group.vector <-
scale_colour_stim_group.numeric <- function(x, alpha = .3, 
  UMAP_mode = F) {
  # col_values <- switch(extract_col_var(colour_var),
  #   'stim_group' = get_stim_group_cols(x), NULL)
  col_values <- get_stim_group_cols(x)
  if (UMAP_mode) {
    col_values['Unexposed in vivo'] <- 'grey80'
  }
  col_values <- alpha(col_values, alpha)
  ggplot2::scale_colour_manual(
    name = 'Stimulus', 
    values = col_values
  )
}
test_plotting <- function()
  print_plot_eval(print(qplot(1:3, 1:3)),
    file.path(img_dir, 'test.pdf'))


scale_fill_stim_group <- function(...)
  UseMethod('scale_fill_stim_group')


scale_fill_stim_group.data.frame <- function(x) {
  scale_fill_stim_group(x$stim_group)
}


scale_fill_stim_group.factor <- function(x) {
  scale_fill_stim_group(e_levels(x))
}


scale_fill_stim_group.character <-
scale_fill_stim_group.vector <-
scale_fill_stim_group.numeric <- function(x) {
  # col_values <- switch(extract_col_var(fill_var),
  #   'stim_group' = get_stim_group_cols(x), NULL)
  col_values <- get_stim_group_cols(x)[x]
  ggplot2::scale_fill_manual(name = 'Stimulus', values = col_values)
}


scale_fill_condition_name <- function(...)
  UseMethod('scale_fill_condition_name')


scale_fill_condition_name.Seurat <- function(so) {
  scale_fill_condition_name(so@meta.data)
}


scale_fill_condition_name.data.frame <-
scale_fill_condition_name.DFrame <- function(x) {
  scale_fill_condition_name(x$condition_name)
}


scale_fill_condition_name.factor <- function(x) {
  scale_fill_condition_name(e_levels(x))
}


scale_fill_condition_name.character <-
scale_fill_condition_name.vector <-
scale_fill_condition_name.numeric <- function(x) {
  # col_values <- switch(extract_col_var(fill_var),
  #   'condition_name' = get_condition_name_cols(x), NULL)
  x <- e_levels(x)
  x_sg <- recover_stim_group(x)
  col_values <- get_stim_group_cols(x_sg)[x_sg]
  # col_values <- set_names(col_values, x[match(names(col_values), x_sg)])
  col_values <- set_names(col_values, x)
  ggplot2::scale_fill_manual(name = 'Stimulus', values = col_values)
}


scale_fill_duration <- function(...) UseMethod('scale_fill_duration')


scale_fill_duration.data.frame <- function(dtf) {
  col_values <- HM_col_funs$duration_discrete(
    dtf[, 'duration', drop=F],
    cn = 'duration'
  )
  ggplot2::scale_fill_manual(name = 'Duration', values = col_values)
}


# scale_fill_duration.factor <- function(v) {
#   scale_fill_duration(dtf$duration)
# }


test_plotting <- function()
  print_plot_eval(print(qplot(1:3, 1:3)),
    file.path(img_dir, 'test.pdf'))


#' Wrapper function for use in pb_harmony.Rmd only
make_umaps <- function(
  p_dat,
  embedding = p_dat$embedding[1] %||% NULL,
  coord_regex =
    glue::glue('^{get_embedding_feat_regex(embedding)}_(1|2)$'),
  fn = NULL,
  fn_app = '',
  colour_var_1 = 'stim_group',
  colour_var_2 = 'tnf_conc',
  width = 17.4, height = 12) {

  stopifnot(is.data.frame(p_dat))
  stopifnot(is.character(coord_regex))

  p1 <- plot_l_umap(
    dtf = p_dat,
    zoom_name = '',
    colour_var = colour_var_1,
    coord_regex = coord_regex,
    # filter_samples = !is.na(MAE) & sample_type == 'sc',
    fn_app = fn_app,
    # zoom_range = experiment %in% c('5029', '6434') & ifn_conc > 1)
    legend_ncol = 2,
    # width = 8.7, height = 10,
    zoom_range = 1 == 1, print_to_file = F
  )

  p2 <- plot_l_umap(
    dtf = p_dat,
    zoom_name = '',
    colour_var = colour_var_2,
    # filter_samples = as.numeric(tnf_conc) > 0 & as.numeric(ifn_conc) > 0,
    coord_regex = coord_regex,
    fn_app = fn_app,
    # zoom_range = experiment %in% c('5029', '6434') & ifn_conc > 1)
    legend_ncol = 3,
    # width = 8.7, height = 10,
    zoom_range = 1 == 1, print_to_file = F
  )

  library(patchwork)
  if (!is.null(fn)) {
    print_plot_eval(print(p1 + p2), filename = fn,
      width = width, height = height, units = 'cm', bg = 'white')
  } else {
    print_umap(
      p = p1 + p2,
      # p = p1,
      dtf = p_dat,
      coord_regex = coord_regex,
      colour_var = NULL,
      filter_name = NULL, zoom_name = NULL,
      fn_app = fn_app, width = width, height = height
    )
  }
}


print_umap <- function(
  p, dtf, coord_regex, colour_var, filter_name,
  zoom_name, fn_app, width = width, height = height) {

  experiments <- gen_exp_string(dtf$exp)
  # coord_regex = "UMAP_1"
  # coord_regex = "PC_1"
  # coord_regex = 'PC_(1|2)'
  # coord_regex = '^harmony_(1|2)$'
  # gsub('[^a-zA-Z]*([a-zA-Z]+).*', '\\1', coord_regex)
  dim_reduc_type <- gsub('[^a-zA-Z]*([a-zA-Z]+).*', '\\1',
    coord_regex)
  axis_vars <- stringr::str_subset(colnames(dtf), coord_regex) %>%
    paste(collapse = '_')
  source_type <- attr(dtf, 'source_type') %||% ''

  o_fn <- file.path(img_dir,
    glue::glue('{source_type}\\
      {dim_reduc_type}\\
      {make_flag(axis_vars)}\\
      {make_flag(experiments)}\\
      {make_flag(colour_var)}\\
      {make_flag(filter_name)}\\
      {make_flag(zoom_name)}\\
      {fn_app}.png'))

  print_plot_eval(print(p), filename = o_fn,
    width = width, height = height, units = 'cm', bg = 'white')

  return(invisible(o_fn))
}

plot_umap_wrapper <- function(
  dtf = p_dat,
  colour_var = 'stim_group',
  filter_name = 'all',
  filter_samples = !is.na(UMAP1),
  zoom_name = 'all',
  zoom_range = 1 == 1,
  point_alpha = .5,
  legend_ncol = 5,
  width = 17.4, height = 20,
  plot_ref_background = T,
  annotation_code = NULL) {
  library(rlang)
}


condition_count_barplot <- function(
  sos, experiment = '5310') {

  dtf <- purrr::imap_dfr(sos, function(so, on) {
    grp_data <- so@meta.data %>%
      order_condition_name() %>%
      { . }
    if ('mouse' %in% colnames(grp_data)) {
        grp_data <- grp_data %>%
          dplyr::group_by(condition_name, mouse)
    } else {
        grp_data <- grp_data %>% dplyr::group_by(condition_name)
    }
    grp_data %>%
      dplyr::summarize(N = n(), .groups = 'drop') %>%
      dplyr::mutate(frac = N / sum(N)) %>%
      dplyr::mutate(type = on)
  }) %>%
  dplyr::mutate(type = factor(type, levels = names(sos))) %>%
  { . }

  # if ('mouse' %in% colnames(dtf)) {
  #   dtf <- dtf %>% dplyr::mutate(
  #     condition_name = paste(condition_name, ' - mouse ', mouse))
  # }

  p <- condition_count <- dtf %>%
    ggplot(aes(x = condition_name, fill = type, y = N)) +
    geom_col(position = 'dodge', width = .8) +
    maartenutils::rotate_x_labels(45) +
    ylab('Number of cells after filtering')

  if ('mouse' %in% colnames(dtf)) {
    p <- p + facet_wrap(~mouse, ncol = 2)
    height = 20
  } else {
    height = 10
  }

  o_fn <- file.path(exp_plot_dir(experiment), glue('{experiment}-\\
      condition_count_barplot.png'))

  ggsave(o_fn, p,
    width = 8.7 / 7 * length(unique(dtf$condition_name)),
    height = height, units = 'cm')

  return(o_fn)
}
# condition_count_barplot(experiment = '6600')

pad_empty_plots <- function(plots, N_in_unit = 7, N_in_grid = 8) {
  plots <- plots[!sapply(plots, is.null)]
  empty_p <- ggplot() + theme_void()
  empty_plots <-
    purrr::map(seq(max(N_in_grid - N_in_unit, 0)), ~empty_p)
  for (i in seq(ceiling(length(plots) / N_in_unit))) {
    # append(list(4,5,6,7), list(3,3,3), after = 3)
    plots <- append(plots, empty_plots,
      after = (i-1)*N_in_grid + N_in_unit)
  }
  return(plots)
}


scale_shape_duration <- function(...)
  UseMethod('scale_shape_duration')


scale_shape_duration.character <-
scale_shape_duration.vector <-
scale_shape_duration.numeric <- function(x) {
  ggplot2::scale_shape_manual(
    name = 'Duration',
    values = rlang::set_names(
      c(15, 16, 17, 17, 18, 18, 0, 1, 19),
      c('2', '6', '12', '16', '24', '44', '48', '96', 'Unknown')
    )[x]
  )
}


scale_shape_duration.factor <- function(x) {
  scale_shape_duration(e_levels(x))
}


scale_shape_duration.data.frame <- function(dtf) {
  stopifnot('duration' %in% colnames(dtf))
  scale_shape_duration(dtf$duration)
}


gen_exp_colors <- function(dtf) {
  l_in_vitro_sc_experiments <- get_in_vitro_sc_experiments(dtf)
  l_in_vivo_sc_experiments <- get_in_vivo_sc_experiments(dtf)
  l_bulk_experiments <- get_bulk_experiments(dtf)

  names_ramp <- function(palette, classes) {
    colorRamp(palette)(seq(0, 1, length.out = length(classes))) %>%
      apply(1, function(r) maartenutils::rgb_to_hex(r[1], r[2], r[3])) %>%
      rlang::set_names(classes)
  }

  palette <- RColorBrewer::brewer.pal(11,'RdYlBu')
  vctrs::vec_c(
    names_ramp(palette[1:2], l_in_vitro_sc_experiments),
    names_ramp(palette[4:5], l_in_vivo_sc_experiments),
    names_ramp(palette[8:11], l_bulk_experiments)
  )
}




print_upset <- function(l,
  min_size = 100L,
  mode = 'intersect',
  fn_app = '',
  height = nrow(mt) / 3,
  width = 8.7,
  mt_mod_code = NULL,
  top_annotation_code = NULL,
  right_annotation_code = NULL,
  ...) {

  m <- ComplexHeatmap::make_comb_mat(l, mode = mode)
  # map_int(l, length)
  # browser()

  if (F) {
    ## Find the binary code corresponding to the intersection of all
    ## concentration coefficients
    conc_idx <-
      which(stringr::str_detect(colnames(attr(m, 'data')), 'conc|sn'))
    conc_intersect_code <-
      rep('0', dim(attr(m, 'data'))[2]) %>%
      { .[conc_idx] = '1'; . } %>%
      paste(collapse = '')
  }

  mt <-
    m %>%
    # {
    #   .[stringr::str_detect(dimnames(.)[[1]],
    #     '^duration', negate = T), ]
    # } %>%
    t() %>%
    { .[comb_size(m) >= min_size] } %>%
    # { .[comb_size(m) >= comb_size(m)[conc_intersect_code]] } %>%
    { . }

  ht_opt(RESET = T)
  set.seed(41)
  print_plot_eval(
    {
      if (!is.null(mt_mod_code)) {
        eval(mt_mod_code)
      }
      HM <- UpSet(
        mt,
        comb_order = order(comb_size(mt)),
        border_gp = gpar(fontsize = 6),
        top_annotation = eval(top_annotation_code),
        right_annotation = eval(right_annotation_code),
        ...
      )
      draw(HM, annotation_legend_side = 'bottom')
    },
    filename = file.path(Sys.getenv('img_dir'),
      glue::glue('upset{fn_app}.pdf')),
    height = height, width = width
  )
}


compare_gene_dynamics_between_experiments <- function(
  gene = 'PNKD',
  experiments = c('5029', '6434', tar_read(bm_e)) %>%
    setdiff('5310') %>%
    sort()) {

  sample_list <- read_preproc_experiments(
    experiments,
    sc_mode = 'pseudobulk'
  )

  settings <-
    tidyr::expand_grid(
      gene = gene,
      experiment = names(sample_list)
    ) %>%
    dplyr::mutate(p_title = glue::glue('{experiment} - {gene}')) %>%
    # tail(n = 1)
    { . }

  plots <-
    settings %>%
      purrr::pmap(function(gene, experiment, p_title) {
        M <- so2M(sample_list[[experiment]])
        if (experiment == '5029') {
          sa <- tar_read(sample_annotation_exp5029)
          merge_cn = 'sample_name'
        } else if (experiment == '6434') {
          sa <- tar_read(sample_annotation_exp6434)
          merge_cn = 'sample_name'
        } else {
          sa <- sample_list[[experiment]]@meta.data
          merge_cn = 'condition_name'
        }
        plot_expression_dynamics(
          fn = gene,
          title = p_title,
          meta_data = sa,
          caching = F,
          merge_cn = merge_cn,
          redo = T,
          leave_out_sn = F,
          lookup_data = M
        )
      })

  dir.create(file.path(img_dir, 'compare_experiments_gene_level'),
    showWarnings = F)
  o_fn <- file.path(img_dir, 'compare_experiments_gene_level',
    glue::glue('feature_dynamics_across_experiments-{gene}.pdf'))

  maartenutils::plot_panel_layout(
    plots = plots,
    labels = NULL,
    nrow = 2, ncol = 2,
    h = 15,
    filename = o_fn
  )

  return(o_fn)
}


gen_clust_object <- function(
  M,
  dist_f = 'pearson',
  clust_method = 'complete',
  reorder_clust = T,
  min_var = NULL,
  M_orig = NULL) {

  pacman::p_load('dendextend')

  if (!is.null(min_var)) {
    idxs <- which(apply(M, 2, var) < min_var)
    if (!is.null(M_orig)) {
      o_idxs <- which(apply(M_orig, 2, var) < min_var)
      idxs <- union(idxs, o_idxs)
    }
  } else {
    idxs <- c()
  }

  if (length(idxs) > 0) {
    message('Dropping: ',
      paste(colnames(M)[idxs], collapse = ', '))
    M <- M[, -idxs]
    if (!is.null(M_orig)) {
      M_orig <- M_orig[, -idxs]
    }
  }

  if (dist_f %in%
    c('spearman', 'pearson', 'kendall')) {
    corM <- cor(M, method = dist_f)
    distM <- dist(1-corM)
  } else if (dist_f %in% c('kendall_unscaled')) {
    stopifnot(!is.null(M_orig))
    corM <- cor(M_orig, method = 'kendall')
    distM <- dist(1-corM)
  } else {
    distM <- dist(t(M), method = dist_f)
  }

  clust_object <- tryCatch(
    hclust(distM, method = clust_method),
    error = function(e) { print(e); NULL })

  clust_object <- as.dendrogram(clust_object)

  if (reorder_clust) {
    wM <- M
    # wM <- t(scale(t(M)))
    # wM[wM < 0] <- 0
    # browser()
    clust_object <- reorder(
      clust_object, 
      wts = colSums(wM),
      agglo.FUN = mean
      # agglo.FUN = median
    )
    # cutree(ref_clust, k = 7)
  }

  return(clust_object)
}


gen_tree_split <- function(clust_object,
  cluster_h = NULL, cluster_k = NULL) {
  if (!is.null(cluster_h)) {
    tree_split <- dendextend::cutree(clust_object, h = cluster_h)
  } else if (!is.null(cluster_k)) {
    tree_split <- dendextend::cutree(clust_object, k = cluster_k)
  } else {
    tree_split <- NULL
  }
  return(tree_split)
}


GE_vs_stim_dur <- function(...) UseMethod('GE_vs_stim_dur')


GE_vs_stim_dur.Seurat <- function(
  so,
  meta_fields = c('stim_group', 'duration'),
  trans = function(x) log2(x + 1),
  # trans = identify,
  genes = VariableFeatures(so),
  ## filter_samples_code does not work
  filter_samples_code = NULL,
  do_calc_ht_size = F,
  filter_samples_mode = 'all',
  ...) {

  M <- trans(as.matrix(subset_feats(GetAssayData(so), genes)))
  # M <- as.matrix(subset_feats(so, genes))
  meta_fields <- intersect(meta_fields, colnames(so@meta.data))

  sa <- extract_sa(so = so, meta_fields = meta_fields)

  GE_vs_stim_dur(
    M = M,
    sa = sa,
    genes = genes,
    trans = NULL,
    do_calc_ht_size = do_calc_ht_size,
    filter_samples_code = filter_samples_code,
    filter_samples_mode = filter_samples_mode,
    ...
  )
}


GE_vs_stim_dur.matrix <- function(
  M,
  sa = NULL,
  plot_name = 'exp5029_var_genes_HM',
  N_hl_genes = 30,
  highlight_genes = c(),
  scale_data = 'Z',
  show_column_names = F,
  show_column_dend = T,
  genes = detected_genes(M),
  min_var = 0.0,
  row_clust_object = NULL,
  row_split = NULL,
  column_clust_object = NULL,
  column_split = NULL,
  # col_dist_f = 'kendall_unscaled',
  # col_dist_f = 'pearson',
  # col_dist_f = 'kendall',
  column_clust_method = 'complete',
  row_dist_f = 'pearson',
  row_clust_method = 'complete',
  col_dist_f = 'euclidean',
  column_cluster_h = NULL,
  column_cluster_k = NULL,
  row_cluster_h = NULL,
  row_cluster_k = NULL,
  filter_samples_code = NULL,
  filter_samples_mode = 'all',
  cell_width_mm = 2,
  cell_height_mm = 0.09,
  max_val = 5,
  top_HM = NULL,
  include_max_e = F,
  do_calc_ht_size = F,
  geneset_labels = NULL,
  draw_c_dend_boxes = F,
  color_name = NULL,
  ...) {

  # filter_samples_code <- rlang::quo(filter_samples_code)
  # condition_call <- substitute(filter_samples_code)
  # eval(condition_call, envir = sa)

  dots <- list(...)
  if ('row_split' %in% names(dots) &&
      is.vector(dots$row_split) &&
      length(dots$row_split) > 1) {
    M <- subset_feats(M,
      genes = intersect(names(dots$row_split), detected_genes(M)))
    dots$row_split <- dots$row_split[
      intersect(names(dots$row_split), detected_genes(M))]
    stopifnot(detected_genes(M) == names(dots$row_split))
  }

  if (T) {
    allowed_genes <- apply(M, 1, var) > min_var
    M <- M[allowed_genes, ]
    if ('row_split' %in% names(dots)) {
      dots$row_split <- dots$row_split[allowed_genes]
    }
    stopifnot(detected_genes(M) == names(dots$row_split))
    if ('ra' %in% names(dots)) {
      dots$ra <- dots$ra[rownames(M), , drop = F]
      # dots$ra <- dots$ra[allowed_genes, ]
    }
    highlight_genes <- intersect(rownames(M), highlight_genes)
    M_orig <- M
    ## Scale the data BEFORE any sample selection is applied
    if (scale_data == 'Z') {
      M <- t(scale(t(M)))
    }

    print(max_val)
    if (!is.null(max_val)) {
      cap_count <- sum(abs(M) > max_val)
      if (cap_count > 0) {
        message('Capping ', cap_count, ' entries to ', max_val)
        locs <- abs(M) > max_val
        M[locs] <- sign(M[locs]) * max_val
      }
    }
    # stopifnot(all(apply(M, 1, var) > 0))
    stopifnot(!any(is.na(M)))
  }

  if (F && !missing(filter_samples_code) &&
      !is.null(filter_samples_code)) {
    idxs <- with(sa, rlang::eval_tidy(filter_samples_code))
    M <- M[, idxs]
    M_orig <- M_orig[, idxs]
    if (!is.null(sa)) {
      sa <- sa[idxs, ]
    }
  }

  ## Sort by stimulus and then duration
  if (FALSE && !is.null(sa)) {
    sa$ri <- 1:nrow(sa)
    idxs <- match(
      dplyr::arrange_at(sa, 
          intersect(c('m_sg', 'stim_group', 'duration'), 
            colnames(sa)))$ri,
      sa$ri
    )
    M <- M[, idxs]
    M_orig <- M_orig[, idxs]
    sa <- sa[idxs, ]
    sa$ri <- NULL
  }

  if (!is.null(sa)) {
    if (filter_samples_mode == 'non_SN') {
      idxs <- stringr::str_detect(sa$stim_group, 'Unstim|IFNy|TNFa')
    } else if (filter_samples_mode == 'SN') {
      idxs <- stringr::str_detect(sa$stim_group, 'Unstim|SN')
    } else if (filter_samples_mode == 'IFNy') {
      idxs <- stringr::str_detect(sa$stim_group, 'Unstim|(IFNy$)')
    } else if (filter_samples_mode == 'TNFa') {
      idxs <- stringr::str_detect(sa$stim_group, 'Unstim|TNFa')
      non_idxs <- stringr::str_detect(sa$stim_group, 'IFNy')
      idxs <- which(idxs & !non_idxs)
    } else {
      idxs <- 1:nrow(sa)
    }
    sa <- sa[idxs, ]
    M <- M[, idxs]
    M_orig <- M_orig[, idxs]
  }

  ## Reassess gene variance. Remove genes that are no longer
  ## informative
  idxs <- which(apply(M, 1, var) >= min_var)
  M <- M[idxs, ]
  M_orig <- M_orig[idxs, ]

  ## Cluster the rows and columns
  if (is.null(row_clust_object)) {
    if (row_clust_method != 'none') {
      row_clust_object <- gen_clust_object(
        M = t(M),
        dist_f = row_dist_f,
        clust_method = row_clust_method,
        M_orig = t(M_orig)
      )
      if (nrow(M) != length(as.hclust(row_clust_object)$labels)) {
        browser()
      }
    } else {
      row_clust_object <- NULL
    }
  }

  if ((!is.null(row_cluster_h) || !is.null(row_cluster_k)) && 
      !is.null(row_clust_object) &&
      row_clust_method != 'none') {
    row_split <- gen_tree_split(
      row_clust_object,
      cluster_h = row_cluster_h,
      cluster_k = min(row_cluster_k, nrow(M), 30)
    )
    if (F) {
      row_split <- gen_tree_split(
        row_clust_object, cluster_h = NULL, cluster_k = 3)
      pacman::p_load('funtimes')
      classes <- geneset_labels %>%
        dplyr::right_join(tibble(gene = names(row_split))) %>%
        dplyr::pull(gs)
      funtimes::purity(classes, row_split)
    }
  } else {
    row_split <- NULL
  }

  if (!is.null(row_split)) {
    M <- subset_feats(M, names(row_split))
    M_orig <- subset_feats(M_orig, names(row_split))
    row_split <- row_split[intersect(names(row_split),
      detected_genes(M))]
  }

  if (is.null(column_clust_object)) {
    if (column_clust_method != 'none') {
      column_clust_object <- 
        gen_clust_object(
          M = M,
          dist_f = col_dist_f,
          clust_method = column_clust_method,
          M_orig = M_orig
        )
      # pacman::p_load('vegan')
      # wM <- M
      wM <- t(scale(t(M)))
      # wM[wM < 0] <- 0
      column_clust_object <- 
        reorder(
          column_clust_object, 
          wts = colSums(wM),
          agglo.FUN = mean
          # agglo.FUN = median
        )
    } else {
      column_clust_object <- NULL
    }
  }

  if ((!is.null(column_cluster_h) || !is.null(column_cluster_k)) && 
      !is.null(column_clust_object)) {
    column_split <- gen_tree_split(column_clust_object,
      cluster_h = column_cluster_h,
      cluster_k = min(column_cluster_k, 30, ncol(M))
    )
    column_split <- sort(column_split)
  } else {
    column_split <- NULL
  }

  if (!is.null(column_split)) {
    M <- M[, names(column_split)]
    M_orig <- M_orig[, names(column_split)]
    if (!is.null(sa)) {
      sa <- sa[names(column_split), ]
    }
  }

  ## Pick highlight genes
  if (!is.null(N_hl_genes) && N_hl_genes > 0) {
    if (!is.null(row_split)) {
      gene_vars <- apply(M_orig, 1, function(x) var(x, na.rm = T))
      # stopifnot(gene_vars )
      gpc <- ceiling(N_hl_genes / max(row_split))
      highlight_genes <-
        map(unique(row_split), function(cl) {
          lgv <- gene_vars[which(row_split == cl)] %>%
            { rank(., ties.method = 'first') } %>%
            { .[. >= max((length(.) - gpc), 1)] } %>%
            { .[seq(gpc)] }
        }) %>%
        map(~names(.x)) %>%
        unlist() %>%
        c(highlight_genes)
    } else {
      extra_highlight_genes <-
        apply(M_orig, 1, function(x) var(x, na.rm = T)) %>%
        { rank(.) } %>% { .[. > (length(.) - N_hl_genes)] } %>%
        names
      highlight_genes <- c(highlight_genes, extra_highlight_genes)
    }
  }

  if (is.null(color_name)) {
    if (scale_data == 'Z') {
      color_name <- paste0(ifelse(scale_data == 'Z', 'Z-scaled\n', ''),
        'gene expression') %>%
      maartenutils::simple_cap(cap_first_word_only = T)
    } else if (scale_data == 'robust') {
      color_name <- paste0(ifelse(scale_data, 'Robust-scaled\n', ''),
        'gene expression') %>%
      maartenutils::simple_cap(cap_first_word_only = T)
    } else {
      color_name <- 'Gene expression'
    }
  }

  plot_args <- list(
    ca = sa,
    ra = NULL,
    row_ann_f = gen_sample_annotation_HM,
    col_ann_f = gen_sample_annotation_HM,
    row_split = row_split,
    highlight_genes = highlight_genes,
    column_split = column_split,
    name = color_name,
    N_hl_genes = N_hl_genes,
    show_column_names = show_column_names,
    show_column_dend = show_column_dend,
    row_dend_reorder = FALSE,
    column_dend_reorder = FALSE,
    cluster_column_slices = FALSE,
    top_HM = top_HM,
    width = unit(ncol(M)*cell_width_mm, 'mm'),
    height = unit(nrow(M)*cell_height_mm, 'mm')
  )

  if (is.null(column_split)) {
    if (is.null(column_clust_object)) {
      plot_args$cluster_columns <- FALSE
    } else {
      plot_args$cluster_columns <- column_clust_object
    }
  } else {
    plot_args$cluster_columns <- TRUE
  }


  if (is.null(row_split)) {
    if (is.null(row_clust_object)) {
      plot_args$cluster_rows <- FALSE
    } else {
      plot_args$cluster_rows <- row_clust_object
    }
  } else {
    plot_args$cluster_rows <- TRUE
  }

  # if (!is.null(plot_args$column_split)) {
  #   plot_args$cluster_columns <- NULL
  # }
  # if (!is.null(plot_args$row_split)) {
  #   plot_args$cluster_rows <- NULL
  # }
  plot_args <- plot_args[setdiff(names(plot_args), names(dots))]
  # class(plot_args$cluster_rows)
  # plot_args$cluster_columns
  # browser()

  # if (include_max_e) {
  #   left_annotation <- rowAnnotation(
  #     `Max expression` = apply(M_orig, 1, max)
  #   )
  # } else {
  #   left_annotation <- NULL
  # }

  if (!is.null(geneset_labels)) {
    plot_args$left_annotation <- NULL
    dots$left_annotation <- NULL
    left_annotation <- rowAnnotation(
      'Gene set' = geneset_labels$gs[
        match(rownames(M), geneset_labels$gene)]
    )
    plot_args$left_annotation <- left_annotation
    dots$geneset_labels <- NULL
    plot_args$geneset_labels <- NULL
  }

  HM <- purrr::exec(gen_HM, M = M, !!!plot_args, !!!dots)

  draw_expr <- rlang::quo({
    # if (!is.null(top_HM)) {
    #   draw(top_HM %v% HM, merge_legend = T, heatmap_legend_side = 'bottom')
    # } else {
    # }
    draw(HM, merge_legend = T, heatmap_legend_side = 'right')
    if (draw_c_dend_boxes) {
      decorate_column_dend(color_name, {
        ind = row_split[order.dendrogram(row_clust_object)]

        N_clusters <- max(row_split)
        first_index = function(l) which(l)[1]
        last_index = function(l) { x = which(l); x[length(x)] }
        last_index = function(l) last(which(l))

        x1 = sapply(1:N_clusters, function(i) first_index(ind == i)) - 1
        x2 = sapply(1:N_clusters, function(i) last_index(ind == i))

        grid.rect(
          # y = 1-x1/length(ind),
          # height = (x2 - x1)/length(ind),
          # just = 'top',
          x = 1-x1/length(ind),
          width = (x2 - x1)/length(ind),
          just = 'right',
          default.units = 'npc',
          gp = gpar(fill = add_alpha(1:5, 0.6), col = NA)
        )
      })
    }
  })

  if (!is.null(highlight_genes)) {
    attr(HM, 'highlight_genes') <- highlight_genes
  }

  if (do_calc_ht_size) {
    HS <- calc_ht_size(draw(HM))
    return(list(draw_expr, HS, HM))
  } else {
    return(HM)
  }
}


calc_ht_size = function(draw_expr, unit = 'inch') {
  pdf(NULL)
  ht = rlang::eval_tidy(draw_expr)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}


calc_ht_size = function(ht, unit = 'inch') {
  pdf(NULL)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}


add_alpha = function(col, alpha = 0.8) {
  rgb(t(col2rgb(col)/255), alpha=alpha)
}


annotate_npc <- function(label, npcx, npcy, ...) {
  ggplot2::annotation_custom(grid::textGrob(
    x = unit(npcx, "npc"), y = unit(npcy, "npc"), label = label, ...))
}


#' Only used in rmd/fig1.Rmd
#'
#'
l_HM <- function(so, genes, fn, make_line_plots = TRUE, 
  out_dir = Sys.getenv('img_dir'), filter_samples_code = NULL, 
  height = 12, ...) {
  # filter_samples_code <- rlang::expr(filter_samples_code)
  dots <- list(...)
  p_args <- list(
    so = so,
    meta_fields = c('stim_group', 'duration'),
    # highlight_genes = 'WARS',
    height = unit(6, 'cm'),
    'filter_samples_code' = filter_samples_code,
    show_column_names = F,
    # column_title = "C%s",
    column_title = ' ',
    row_title = ' ',
    row_gap = unit(.5, 'mm'),
    show_row_dend = F,
    genes = genes
  )
  p_args <- p_args[setdiff(names(p_args), names(dots))]
  # p_args$filter_samples_code
  # DE <- purrr::exec(GE_vs_stim_dur, !!!p_args, !!!dots)
  HM <- do.call(GE_vs_stim_dur, c(p_args, dots))

  print_plot_eval(draw(HM, merge_legends = TRUE),
    width = 17.4, height = height,
    filename = file.path(out_dir, fn)
  )

  if (make_line_plots) {
    exp_dynamics_panel(
      features = attr(HM, 'highlight_genes'),
      merge_cn = 'sample_name',
      out_dir = out_dir,
      lookup_data = as.matrix(tar_read('kallisto_5029')),
      version = fn %>%
        stringr::str_replace('.pdf', '') %>%
        paste0('_hl_genes'),
      redo = F,
      leave_out_sn = T
    )
  }
}


gen_coord_ranges_fun <- function(
  dtf,
  break_interval = .5,
  column_selector = matches('UMAP')) {

  coord_ranges <- tryCatch({
    dtf %>%
      dplyr::select(!! enquo(column_selector)) %>%
      purrr::map(range)
  }, error = function(e) { print(e) })

  out <- function(x) {
    if (is.null(x)) return(NULL)
    out <- x +
      # scale_x_continuous(name = names(coord_ranges)[1],
      #   breaks = function(y)
      #     seq(round(min(y), 1), round(max(y), 1), by = break_interval)) +
      # scale_y_continuous(name = names(coord_ranges)[2],
      #   breaks = function(y)
      #     seq(round(min(y), 1), round(max(y), 1), by = break_interval)) +
      xlab(names(coord_ranges)[1]) +
      ylab(names(coord_ranges)[2]) +
      # ggplot2::expand_limits(
      #   x = coord_ranges[[1]],
      #   y = coord_ranges[[2]]
      #   # x = -c(-10, 10),
      #   # y = -c(-10, 10)
      # ) +
      ggplot2::coord_cartesian(
        xlim = coord_ranges[[1]],
        ylim = coord_ranges[[2]]
      )
    return(out)
  }

  return(out)
}


gen_color_ranges_fun <- function(
  dtf,
  column_selector = matches('UMAP')) {

  cr <- tryCatch({
    dtf %>%
      dplyr::select(!! enquo(column_selector)) %>%
      purrr::map(range)
  }, error = function(e) { print(e) })

  ggplot2::expand_limits(colour = cr)
}


