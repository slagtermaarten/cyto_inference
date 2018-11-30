all_prep_funcs <- list(
  'vanilla' = function(x) x,
  'max-norm' = function(x) x / max(x, na.rm = T)
)


all_integrate_funcs <- list(
  'sum' = function(x) sum(x, na.rm = T),
  'max' = function(x) max(x, na.rm = T)
)


all_weighting_funcs <- list(
  'unweighted' = function(x, genes, gs_name, tV) x,
  ## Max-diff-weighted
  'mdw' = function(x, genes, gs_name, tV) {
    ## Strip off time informativeness indicators
    gs_name <- gsub('(.*)-\\d+h*', '\\1', gs_name)
    cn <- switch(gs_name,
      'TNFa' = 'max_diff_tnf_from_us',
      'TNFa.6h.max' = 'max_diff_tnf_from_us',
      'TNFa.early' = 'max_diff_tnf_from_us',
      'TNFa.late' = 'max_diff_tnf_from_us',
      'TNFa.tafelberg' = 'max_diff_tnf_from_us',
      'TNFa.plateau' = 'max_diff_tnf_from_us',
      'IFNy' = 'max_diff_ifn_from_us',
      'IFNg.late' = 'max_diff_ifn_from_us',
      'IFNy.late' = 'max_diff_ifn_from_us',
      'IFNg.tafelberg' = 'max_diff_ifn_from_us',
      'IFNy.plateau' = 'max_diff_ifn_from_us',
      'synergy' = 'max_diff_combo_from_us',
      'CXCL10' = 'max_diff_combo_from_us',
      'none'
    )
    # if (gs_name == 'CXCL10') browser()
    gs_data <- targets::tar_read(gs_data_step)
    if (cn != 'none') {
      w <- gs_data %>%
        dplyr::right_join(tibble(gene = genes), by = 'gene') %>%
        pull(.env[['cn']])
      w[is.na(w)] <- 0
      if (length(w) != length(genes)) browser()
    } else {
      rlang::warn(
        paste('Could not find weights for geneset', gs_name))
      w <- rep(1, nrow(x))
    }
    x[is.na(x)] <- 0
    if (length(w) > 1) {
      out <- diag(w) %*% x
    } else {
      out <- w * x
    }
    if (!all(out[w > 0,][x[w > 0,] > 0] > 0)) browser()
    stopifnot(all(out[x < 1e-32] < 1e-32))
    out
  }, 
  ## T-weighted
  'tw' = function(x, genes, gs_name, tV) {
    stopifnot(!is.null(tV))
    w <- 2^tV[genes]
    w[is.na(w)] <- 1
    if (length(w) > 1) {
      out <- diag(w) %*% x
    } else {
      out <- w * x
    }
    if (!all(out[w > 0,][x[w > 0,] > 0] > 0)) browser()
    stopifnot(all(out[x < 1e-32] < 1e-32))
    out
  }
)


ComputeGeneSetScores <- function(...) UseMethod('ComputeGeneSetScores')


ComputeGeneSetScores.matrix <- function(
  M,
  genesets,
  tM = NULL,
  tV = NULL,
  prep_funcs = all_prep_funcs['vanilla'],
  # weighting_funcs = all_weighting_funcs['mdw'],
  weighting_funcs = all_weighting_funcs[c('tw', 'unweighted')],
  integrate_funcs = all_integrate_funcs['sum'], 
  simplify_names = TRUE) {

  if (is.null(tV) && !is.null(tM)) {
    tV <- apply(tM, 1, max, na.rm = T)
  }

  scores_list <- 
    purrr::imap(genesets, function(gene_set, gs_name) {
      purrr::map(prep_funcs, function(prep_f) {
        purrr::map(weighting_funcs, function(weighting_f) {
          purrr::map(integrate_funcs, function(integrate_f) {
            M <- subset_feats(M, gene_set)
            invalidity_test <- tryCatch(
              is.na(M) || is.null(M) || nrow(M) == 0,
              error = function(e) { print(e); browser() })
            # if (is.na(invalidity_test)) browser()
            if (invalidity_test) return(NULL)
            # if (is.null(M) || nrow(M) == 0) browser()
            M <- tryCatch(t(apply(M, 1, prep_f)),
              error = function(e) { print(e); browser() })
            M <- tryCatch(weighting_f(M, genes = rownames(M),
                gs_name = gs_name, tV = tV),
              error = function(e) { print(e); browser() })
            out <- apply(M, 2, integrate_f)
          })
        })
      })
    })

  score_M <- scores_list %>%
    unlist(recursive = F) %>%
    unlist(recursive = F) %>%
    unlist(recursive = F) %>%
    dplyr::bind_cols() %>%
    as.matrix %>%
    t() %>%
    set_colnames(colnames(M)) %>%
    ## Seurat objects do not allow underscores in feature names
    { set_rownames(.,
      stringr::str_replace(rownames(.), '_', '-')) } %>%
    { . }
  
  if (simplify_names && 
      length(prep_funcs) == 1 && 
      length(weighting_funcs) == 1 &&
      length(integrate_funcs) == 1) {
    cp <- glue::glue('{names(prep_funcs)}.\\
      {names(weighting_funcs)}.{names(integrate_funcs)}')
    score_M <- score_M %>%
        { .[grepl(cp, rownames(.)), , drop = F] } %>%
        { set_rownames(., gsub(glue::glue('\\.{cp}'), '', rownames(.))) }
  } else if (simplify_names) {
    cpf <- extract_common_prefix(rownames(score_M))
    rownames(score_M) <- stringr::str_replace_all(
      rownames(score_M), cpf, '')
  }

  return(score_M)
}


ComputeGeneSetScores.Seurat <- function(
  so,
  assay = DefaultAssay(so), 
  genesets = tar_read(all_genesets),
  # genesets = tar_read(default_gene_sets),
  prep_funcs = all_prep_funcs['vanilla'],
  weighting_funcs = all_weighting_funcs[c('tw', 'unweighted')],
  integrate_funcs = all_integrate_funcs['sum'], 
  datatype = 'data',
  tM = NULL,
  tV = NULL,
  simplify_names = TRUE) {

  library(Seurat)

  score_M <- ComputeGeneSetScores(
    M = so2M(so, assay = assay, datatype = datatype), 
    genesets = genesets,
    tM = tM, tV = tV,
    prep_funcs = prep_funcs,
    weighting_funcs = weighting_funcs,
    integrate_funcs = integrate_funcs,
    simplify_names = simplify_names
  )

  so[['GS']] <- CreateAssayObject(data = score_M)
  return(so)
}


min_max_scale_assay <- function(so, assay = 'GS') {
  stopifnot(assay %in% names(so@assays))
  so[[glue('{assay}_MM')]] <-
    t(apply(so[[assay]][,], 1,
        function(x) 1e2 * min_max_scaling(x))) %>%
    { CreateAssayObject(data = .) }
  return(so)
}


subset_stim_extract_gs <- function(...)
  UseMethod('subset_stim_extract_gs')


subset_stim_extract_gs.Seurat <- function(so,
  stim_group = NULL, duration = NULL) {
  library(Seurat)
  if (!is.null(stim_group)) {
    so_subs <- so[, so@meta.data$stim_group == stim_group]
  } else if (!is.null(duration)) {
    so_subs <- so[, so@meta.data$duration == duration]
  } else {
    so_subs <- so
  }
  so_subs <- ComputeGeneSetScores(
    so = so_subs,
    genesets = read_all_genesets())
  gs2dtf(so_subs)
}


##################   DEPRECATED   ##############################


ComputeGSFCs <- function(so, assay = 'SCT') {
  lt <- function(x) log2(x + .1)
  comps <- tibble::tribble(
    ~name, ~x, ~y,
    'IFNy-vs-TNFa', 'IFNy', 'TNFa',
    'IFNy6h-vs-IFNy', 'IFNy-6', 'IFNy',
    'IFNy12h-vs-IFNy', 'IFNy-12', 'IFNy',
    'IFNy24h-vs-IFNy', 'IFNy-24', 'IFNy',
    'IFNy24h-vs-IFNy6h', 'IFNy-24', 'IFNy-6',
    'TNFa24h-vs-TNFa', 'TNFa-24', 'TNFa',
    'TNFa24h-vs-2h', 'TNFa-24', 'TNFa-2',
    'synergy24h-vs-IFNy24h', 'synergy-24', 'IFNy-24') %>%
    dplyr::mutate(scores = purrr::pmap(., function(name, x, y) {
      if (any(c(x, y) %nin% rownames(so[[assay]]))) {
        return(NULL)
      }
      lt(so[[assay]][x,]) - lt(so[['GStime']][y,])
    })) %>%
    {
      comb_R <- as.matrix(bind_rows(map(.$scores, as_tibble)))
      set_rownames(comb_R, .$name[!map_lgl(.$scores, is.null)])
    } %>%
    # { t(apply(., 1, function(x) ifelse(!is.finite(x), max(x), x))) } %>%
    # { t(apply(., 1, function(x) ifelse(is.na(x), 0, x))) } %>%
    { . }
  so[[glue('{assay}comps')]] <-
    Seurat::CreateAssayObject(data = comps)
  return(so)
}


#' Edeprecated?!
#'
#'
filtering_to_character <- function(fs) {
  as.character(rlang::quo_get_expr(filtering_settings)) %>%
    { .[c(2,1,3)] } %>%
    paste(collapse = '')
}


compute_gs_ratios <- function(gs_dtf) {
  feats <- attr(gs_dtf, 'feats')
  if (F) {
    feats_comps <- tidyr::expand_grid(f1 = feats, f2 = feats) %>%
      dplyr::filter(as.integer(f1) > as.integer(f2)) %>%
      { . }
  } else {
    feats_comps <- tribble(
      ~f1, ~f2,
      'TNFa', 'IFNy',
      'TNFa', 'synergy',
      'IFNy', 'synergy',
      'TNFa_early', 'TNFa',
      'TNFa_late', 'TNFa_early',
      'IFNy_late', 'IFNy_plateau')
  }

  lt <- function(x) log2(x + 1)
  ratios <- feats_comps %>%
    purrr::pmap_dfc(function(f1, f2) {
      if (f1 %nin% colnames(gs_dtf)) stop('f1 not found')
      if (f2 %nin% colnames(gs_dtf)) stop('f2 not found')
      lt(gs_dtf[[f1]]) - lt(gs_dtf[[f2]])
      # gs_dtf[[f1]] / gs_dtf[[f2]]
    }) %>%
    suppressMessages() %>%
    set_colnames(purrr::pmap_chr(feats_comps, ~glue('{.x}_vs_{.y}')))

  ann <- subselect_ann(gs_dtf)
  R <- SummarizedExperiment::SummarizedExperiment(
    assays = list(ratio = t(ratios)),
    colData = ann)

  return(R)
}


train_tree <- function(
  so = bulk_so,
  filtering_settings = std_bulk_filtering,
  N_reps = 10L, CV = .1) {

  gs_ratio <- ratio_prep(
    so = so,
    filtering_settings = filtering_settings
  )

  gs_ratio <- gs_ratio %>%
    randomly_duplicate_samples(N_reps = N_reps, CV = CV)

  # pacman::p_load('partykit')
  # pacman::p_load('ggparty')
  # pacman::p_load('treeheatr')
  form <-
    # glue('condition_name ~ {paste(colnames(ratios),
    # glue('stim_group ~ {paste(colnames(ratios),
    glue('stim_group + duration ~ {paste(colnames(ratios),
      collapse = \' + \')}') %>%
    as.formula()
  tree <- partykit::ctree(form, data = gs_ratio,
    control = partykit::ctree_control(minbucket = 1*N_reps))

  return(tree)
}


std_bulk_filtering <-
  rlang::quos(ifn_conc != 0.01, tnf_conc != 0.1, sn_dilution != 20000)


ratio_prep <- function(
  so = bulk_so,
  merge_id = c('condition_name', 'seurat_clusters'),
  filtering_settings = std_bulk_filtering) {
  library(SummarizedExperiment)

  merge_id <- colnames(so@meta.data) %>%
    intersect(merge_id) %>%
    match.arg(choices = merge_id)
  sample_id <- rownames(so@meta.data)
  if (any(stringr::str_length(sample_id) >= 100))
    stop('Suspicioulsy long sample names')

  if (is.null(so@assays[['GS']])) {
    so <- ComputeGeneSetScores(
      so = so,
      genesets = read_all_genesets()
    )
  }

  gs_dtf <- gs2dtf(so, merge_id = merge_id) %>%
    apply_data_filtering(filtering_settings)

  gs_ratio <- compute_gs_ratios(gs_dtf)
  metadata(gs_ratio) <- list(
      'merge_id' = merge_id,
      'experiment' = unique(gs_dtf$exp)
    )
  return(gs_ratio)
}


ratio_prep_pseudobulk <- function(scpb_so) {
  library(SummarizedExperiment)
  pb_query_obj <-
    ratio_prep(
      scpb_so,
      merge_id = 'seurat_clusters',
      filtering_settings = NULL
    )
  metadata(pb_query_obj)$experiment <-
    paste(metadata(pb_query_obj)$experiment,
      '-PB', sep = '')
  return(pb_query_obj)
}


#' Failed experiment, probably due to grouping of different
#' durations/concentrations into one bin, confusing the model
#'
#'
train_stim_id <- function(
  ref_obj, N_reps = 10L, CV = .1, mod_type = 'svm') {
  library(tidymodels)
  library(parsnip)

  if (class(ref_obj) == 'SummarizedExperiment') {
    ref_obj <- se2tibble(ref_obj)
  }

  ref_obj <- ref_obj %>%
    dplyr::mutate(across(
      matches('_conc|dilution|duration'), numerify_factor))

  ## Only keep extreme 'archetypes' for initial clustering
  ref_obj <- apply_data_filtering(ref_obj,
    rlang::quos(
      ifn_conc %in% range(ifn_conc),
      tnf_conc %in% range(tnf_conc),
      sn_dilution %in% range(sn_dilution)
    ), verbose = F)

  ## Strip away concentration information from stim_group
  ref_obj <- add_stim(ref_obj)

  # wflow <- train_mod(dtf = ref_obj, mod_type = 'svm')
  # names(wflow)
  # wflow$fit$actions$model
  # str(wflow$fit$actions$model)
  f_mod <- train_mod(
    ref_obj = ref_obj,
    mod_type = mod_type,
    response_var = 'stim',
    var_selector = contains('_vs_'),
    N_reps = N_reps, CV = CV)

  return(f_mod)
}


train_stim_conc_duration <- function(
  ref_obj, N_reps = 10L, CV = .1, mod_type = 'svm') {
  library(tidymodels)
  library(parsnip)

  f_mod <- train_mod(
    ref_obj = ref_obj,
    mod_type = mod_type,
    response_var = 'condition_name',
    var_selector = contains('_vs_'),
    N_reps = N_reps, CV = CV)

  return(f_mod)
}


score_ratio_sim <- function(query_obj, ref_obj, f_mod) {
  SM <- predict_stim(query_obj = query_obj, ref_obj = ref_obj,
    f_mod = f_mod)
  # print_plot_eval(draw(Heatmap(assays(SM)$probs)),
  #   filename = file.path(img_dir, 'testHM.pdf'))
  # nrow(rowData(SM))

  sim_scores <- score_dist_M(
    distM_w = t(assays(SM)$probs),
    ref_ann = rowData(SM),
    query_ann = colData(SM),
    v_mod_query = function(x) x
  )
}


#' Ingest stim_group or condition_name variables and strip them down
#' the stimulus purely
#'
#'
reduce_to_stim <- function(v) {
  as.character(v) %>%
    stringr::str_replace_all('\\d*\\.*\\d+ ng/ml ', '') %>%
    stringr::str_replace('1\\/\\d+ ', '') %>%
    stringr::str_replace(' - \\d+h$', '') %>%
    { . }
}


add_stim <- function(dtf) {
  if ('stim_group' %nin% colnames(dtf)) {
    return(dtf)
  }
  dtf %>% dplyr::mutate(stim = reduce_to_stim(stim_group))
}


stim_class_prediction_panel <- function(obj, f_mod) {
  pacman::p_load('RcppRoll')
  benchmark_dtf <-
    se2tibble(obj) %>%
    add_stim() %>%
    dplyr::mutate(
      pred_class = as.character(unlist(predict(f_mod, .)))) %>%
    dplyr::filter(stim %in% unique(pred_class)) %>%
    dplyr::mutate(stim = factor(stim)) %>%
    dplyr::mutate(pred_class =
      factor(pred_class, levels = levels(stim))) %>%
    { dplyr::bind_cols(., predict(f_mod, ., type = 'prob')) } %>%
    # debug_pipe() %>%
    dplyr::mutate(max_pred_probability =
      apply(dplyr::select(., matches('\\.pred_')), 1, max)) %>%
    # dplyr::mutate(true_stim_pred_probability =
    #   dplyr::select(., matches('\\.pred_'))[1:n(), glue('.pred_{stim}')]) %>%
    # dplyr::mutate(pred_probability_true_stim =
    #   imap_dbl(unname(stim), ~unlist(pred_probs[.y, .x]))) %>%
    dplyr::mutate(pred_correct = pred_class == stim) %>%
    dplyr::arrange(max_pred_probability) %>%
    dplyr::group_by(stim_group) %>%
    dplyr::mutate(pred_correct_roll_avg =
      RcppRoll::roll_mean(as.integer(pred_correct), n = n() / 50,
        align = 'right', fill = NA)
    ) %>%
    dplyr::ungroup() %>%
    # dplyr::bind_cols(pred_probs)
    { . }

  pred_cols <-
    purrr::map_chr(levels(benchmark_dtf$stim),
      ~glue('.pred_{.x}')) %>%
    # rlang::syms()
    { . }
  p0 <- benchmark_dtf %>%
    # exec(yardstick::roc_curve, quo(stim), pred_cols) %>%
    yardstick::roc_curve(stim, any_of(pred_cols)) %>%
    autoplot() + theme_cyto_inf()
  # p0 <- benchmark_dtf %>%
  #   yardstick::roc_curve(stim, `.pred_IFNy`, `.pred_IFNy TNFa`,
  #   `.pred_TNFa`, `.pred_Unstimulated in vitro`) %>%
  #   autoplot() + theme_cyto_inf()
  p1 <-
    ggplot(pred_dat, aes(y = pred_probability, x = pred_correct,
        colour = stim)) +
    geom_jitter(width = .2, alpha = .5) +
    # geom_smooth() +
    geom_violin(fill = NA)
  p1 <-
    pred_dat %>%
    dplyr::filter(as.character(stim) %in%
      unique(as.character(pred_class))) %>%
    ggplot(aes(x = pred_probability,
        y = pred_correct_roll_avg, colour = stim)) +
    geom_point(alpha = .1) +
    scale_colour_discrete(name = 'Stimulus') +
    # geom_line(alpha = .5)
    geom_smooth(span = 2, method = 'loess') +
    labs(
      x = 'Prediction probability of highest stimulation',
      y = 'Fraction correct predictions\n(rolling average)'
    )
  p2 <-
    pred_dat %>%
      dplyr::filter(as.character(stim) %in%
        unique(as.character(pred_class))) %>%
      # dplyr::mutate(stim = droplevels(stim)) %>%
      # dplyr::mutate(pred_class = droplevels(pred_class)) %>%
      dplyr::mutate(stim = factor(stim)) %>%
      dplyr::mutate(pred_class =
        factor(pred_class, levels = levels(stim))) %>%
      dplyr::ungroup() %>%
      yardstick::conf_mat(stim, pred_class) %>%
      # tidy()
      autoplot(type = 'heatmap') +
      theme_cyto_inf() +
      guides(fill = 'none')
      # { . }
  print_plot_eval(print((p1 + plot_spacer()) / (p0 + p2)),
    filename = file.path(img_dir,
      glue::glue('ratio_pred_performance-{metadata(obj)$experiment}.png')),
    width = 17.4, height = 10, units = 'cm')
}


gs2dtf <- function(so, merge_id = c('stim_group', 'duration')) {
  stopifnot('GS' %in% names(so@assays))
  gs_names <- rownames(so@assays$GS) %>%
    stringr::str_subset('vanilla.mdw.sum')
  coi <- c(merge_id, 'stim_group', 'duration', 'condition_name',
    'tnf_conc', 'ifn_conc', 'sn_dilution', 'exp', 'experiment') %>%
    purrr::map(~c(.x, glue('{.x}_med'), glue('{.x}_certainty'))) %>%
    purrr::flatten_chr()

  t_dat <- suppressWarnings(
    Seurat::FetchData(so,
      c(intersect(coi, colnames(so@meta.data)), gs_names))) %>%
    dplyr::rename_with(~stringr::str_replace(.x, 'gs_', '')) %>%
    dplyr::rename_with(~stringr::str_replace(.x,
        '.vanilla.mdw.sum', '')) %>%
    dplyr::rename_with(~stringr::str_replace_all(.x,
        '\\.|-', '_')) %>%
    # dplyr::mutate(stim_group = as.factor(stim_group)) %>%
    # dplyr::mutate(duration = as.factor(duration)) %>%
    tibble::rownames_to_column('sample_name') %>%
    { . }

  feats <- names(which(purrr::map_lgl(t_dat, is.numeric))) %>%
    setdiff(c('duration', 'tnf_conc', 'ifn_conc'))
  feat_order <- t_dat %>%
    dplyr::select(any_of(feats)) %>%
    apply(2, mean) %>%
    sort %>%
    names
  t_dat[, feats] <- t_dat[, feats] / t_dat[, 'HK_genes']
  feats <- factor(feats, levels = feat_order)
  attr(t_dat, 'feats') <- feats
  return(t_dat)
}


scvi2dtf <- function(so, merge_id = c('stim_group', 'duration')) {
  stopifnot('scvi' %in% names(so@assays))
  coi <- c(merge_id, 'stim_group', 'duration', 'condition_name',
    'tnf_conc', 'ifn_conc', 'sn_dilution', 'exp', 'experiment')

  t_dat <- suppressWarnings(
    Seurat::FetchData(so,
      c(intersect(coi, colnames(so@meta.data)),
        rownames(so[['scvi']])))) %>%
    dplyr::rename_with(
      ~stringr::str_replace_all(.x, '\\.|-', '_')) %>%
    # dplyr::mutate(stim_group = as.factor(stim_group)) %>%
    # dplyr::mutate(duration = as.factor(duration)) %>%
    tibble::rownames_to_column('sample_name') %>%
    { . }

  if ('experiment' %nin% colnames(t_dat)) {
    t_dat$experiment <- t_dat$exp
  }

  return(t_dat)
}


#' Get a percentile for a subset of rows
#'
#'
retrieve_anchor_thresh <- function(
  dtf, anchor_gs = 'IFNy', perc = .9, 
  ref_condition = tibble('stim_group' = 'Unexposed in vivo')) {

  stopifnot(anchor_gs %in% colnames(dtf))
  l_group_vars <- c('stim_group', 'duration') %>%
    intersect(colnames(dtf))
  q9s <- 
    dtf %>%
    dplyr::mutate(across(any_of(c('stim_group')), as.character)) %>%
    dplyr::mutate(across(any_of(c('duration')), factor_to_numeric)) %>%
    dplyr::group_by_at(l_group_vars) %>%
    dplyr::summarize(q9 = 
      quantile(.data[[anchor_gs]], perc, na.rm = T))

  if (!is.null(ref_condition)) {
    g1 <- q9s %>% 
      dplyr::right_join(
        ref_condition, 
        by = intersect(l_group_vars, colnames(ref_condition))
      ) %>% 
      { . }
    return(g1)
  } else {
    return(q9s)
  }
}
# retrieve_anchor_thresh(
#   dtf = GS_p_dat, anchor_gs = 'IFNy', perc = .9)


filter_anchor_thresh <- function(
  dtf, anchor_gs = 'IFNy', perc = .9, 
  mode = '>=', 
  ref_condition = tibble('stim_group' = 'Unexposed in vivo')) {

  if (is.na(anchor_gs) || is.null(anchor_gs) || 
      is.na(perc) || is.null(perc)) {
    return(dtf)
  }

  g1 <- retrieve_anchor_thresh(
    dtf = dtf, 
    anchor_gs = anchor_gs, 
    ref_condition = ref_condition, 
    perc = perc
  )

  if (mode == '>=') {
    out <- dplyr::filter(dtf, .data[[anchor_gs]] >= g1$q9)
  } else if (mode == '<') {
    out <- dplyr::filter(dtf, .data[[anchor_gs]] < g1$q9)
  }

  return(out)
}


shuffle_rows <- function(dtf) {
  stopifnot(is.data.frame(dtf))
  si <- sample(1:nrow(dtf), replace = FALSE)
  return(dtf[si, ])
}


normalize_gs <- function(
  dtf, 
  var_cols = colnames(dtf)[1:17],
  ref_condition = tibble('stim_group' = 'Unexposed in vivo')) {

  if (is.null(ref_condition)) {
    return(dtf)
  }

  l_group_vars <- c('stim_group', 'duration') %>%
    intersect(colnames(dtf))

  norm_factors <- 
    dtf %>%
    # group_by_at(colnames(ref_condition)) %>%
    dplyr::right_join(ref_condition) %>%
    # dplyr::group_by_at(l_group_vars) %>%
    # dplyr::select(any_of(var_cols)) %>%
    # dplyr::summarize(q5 = 
    #   quantile(.data[[anchor_gs]], .5, na.rm = T))
    dplyr::summarize(across(any_of(var_cols), 
      list(
        'med' = ~quantile(.x, .5, na.rm = T), 
        'iqr' = ~IQR(.x, na.rm = T)
      ))
    ) %>%
    unlist() %>%
    { . }

  ## Extract components from summarize result, could probably be done
  ## more elegantly
  meds <- 
    norm_factors[seq(1, length(norm_factors), 
      length.out = length(var_cols))] %>%
    { .[!stringr::str_detect(names(.), 'iqr')] }
  iqrs <- norm_factors[seq(2, length(norm_factors), 
    length.out = length(var_cols))] %>%
    # { .[stringr::str_detect(names(.), 'iqr')] } %>%
    { . }

  meta_data <- dtf %>%
    dplyr::select(!any_of(var_cols))
  norm_data <- dtf %>%
    dplyr::select(any_of(var_cols)) %>%
    apply(1, function(x) (x - meds) / iqrs) %>%
    t() %>%
    # rbindlist()
    { . }

  dropped_genesets <- names(which(!iqrs > 0)) %>%
    stringr::str_replace('_iqr', '')
  if (length(dropped_genesets) > 0) {
    rlang::warn(paste0('Dropping: ', 
        paste(dropped_genesets, collapse = ', ')))
    norm_data <- norm_data[, iqrs > 0]
  }

  out <- bind_cols(meta_data, norm_data)
  return(out)
}


compute_perc_labels <- function(dtf, qs, x_off=.05, y_off=.05) {
  dtf %>%
    dplyr::mutate(
      IFNy_hi = IFNy > qs$IFNy,
      TNFa_hi = TNFa > qs$TNFa,
    ) %>%
    dplyr::group_by(cn_simple) %>%
    dplyr::summarize(
      Q1 = mean(!TNFa_hi & !IFNy_hi),
      Q2 = mean(TNFa_hi & !IFNy_hi),
      Q3 = mean(IFNy_hi & !TNFa_hi),
      Q4 = mean(IFNy_hi & TNFa_hi)
    ) %>%
    tidyr::pivot_longer(Q1:Q4) %>%
    # dplyr::mutate(value = round(100 * value, 1)) %>%
    dplyr::right_join(
      tibble(
        name = c('Q1', 'Q2', 'Q3', 'Q4'),
        npc_x = c(x_off, x_off, 1-x_off, 1-x_off),
        npc_y = c(y_off, 1-y_off, y_off, 1-y_off)
      ),
      by = 'name'
    )
}


#' 
#' Q2 Q4
#' Q1 Q3
compute_perc_labels <- function(dtf, qs, 
  x_var = 'IFNy', y_var = 'TNFa', x_off = 0.05, y_off = 0.02) {
  dtf %>%
    dplyr::mutate(
      x_hi = .data[[x_var]] > qs[[x_var]],
      y_hi = .data[[y_var]] > qs[[y_var]]
    ) %>%
    dplyr::summarize(
      Q1 = mean(!y_hi & !x_hi),
      Q2 = mean(y_hi & !x_hi),
      Q3 = mean(x_hi & !y_hi),
      Q4 = mean(x_hi & y_hi)
    ) %>%
    tidyr::pivot_longer(Q1:Q4) %>%
    # dplyr::mutate(value = round(100 * value, 1)) %>%
    dplyr::right_join(
      tibble(
        name = c('Q1', 'Q2', 'Q3', 'Q4'),
        npc_x = c(x_off, x_off, 1-x_off, 1-x_off),
        npc_y = c(y_off, 1-y_off, y_off, 1-y_off)
      ),
      by = 'name'
    )
}
