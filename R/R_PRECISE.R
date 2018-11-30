# if (!stringr::str_detect(reticulate::py_config()$python, 'r4')) {
#   reticulate::use_condaenv(condaenv='r4', required = T)
# }
# rm(precise)

# ls(precise)
# ls(precise$pyobj)


#' Compute PRECISE reconstruction error
#'
#'
compute_RE <- function(M, PVs, weighted_by_gene_loadings = T,
  epsilon = .Machine$double.eps) {
  M_e <- M %*% t(PVs$source_components_)
  M_r <- M_e %*% PVs$source_components_
  D <- (M_r - M) / (M+epsilon)
  # D[M == 0] <- 0
  # D[M == 0] <- NA
  if (weighted_by_gene_loadings) {
    var_e <- var_explained(PVs, M)
    WL <- t(PVs$source_components_) %*% diag(var_e$source)
    gene_weights <- rowSums(abs(WL)) %>% { . / sum(.) }
    CM <- diag(ncol(D) * gene_weights)
    WD <- D %*% CM
    return(WD)
  } else {
    return(D)
  }
}




#' Compute explained variance for a set of PVs and data matrix
#'
#'
var_explained <- function(PVs, X) {
  total_var <- sum(apply(X, 2, var))
  proj_source <- X %*% t(PVs$source_components_)
  proj_target <- X %*% t(PVs$target_components_)
  var_source = apply(proj_source, 2, var) / total_var
  var_target = apply(proj_target, 2, var) / total_var
  return(tibble(
      'pv_idx' = as.character(glue::glue('PV{seq_along(var_source)}')),
      'source' = var_source,
      'target' = var_target))
}


gen_fitter <- function(N_source = 10L, N_target = 10L) {
  precise <- reticulate::import('precise')
  PV_fitter <- precise$PVComputation(
    n_pv = as.integer(min(N_source, N_target)),
    n_factors = list(
      'source' = as.integer(N_source),
      'target' = as.integer(N_target)),
    dim_reduction = 'pca'
  )
}


gen_CF_fitter <- function(X_source, X_target,
  N_source = 10L, N_target = 10L, N_PV = 5L) {
  n_factors <- list(
    # 'source' = max(nrow(X_source)-1L, as.integer(N_source)),
    # 'target' = max(nrow(X_target)-1L, as.integer(N_target))
    'source' = 10L,
    'target' = 10L
  )
  precise <- reticulate::import('precise')
  clf <- precise$ConsensusRepresentation(
    source_data = X_source,
    target_data = X_target,
    n_factors = n_factors,
    n_pv = min(unlist(n_factors), N_PV),
    n_representations = 100L,
    use_data = FALSE,
    mean_center = FALSE,
    std_unit = FALSE
  )
  out <- clf$fit(t(X_source))
  return(out)
}


PRECISE_integrate <- R6::R6Class('PRECISE_integrate',
  lock_objects = FALSE,
  lock_class = FALSE,
  public = list(
  print = function() {
    pf <- print
    pf <- function(x) cat(x, '\n')
    pf(self$sc_so@meta.data$exp[1])
    pf(self$norm_method)
    pf(self$N_CFs)
  },
  initialize = function(
    bulk_so,
    sc_so,
    N_PV_bulk = 10L,
    N_PV_sc = 10) {

    self$N_PV_bulk <- N_PV_bulk
    self$N_PV_sc <- N_PV_sc
    self$N_CFs <- min(self$N_PV_bulk, self$N_PV_sc)
    self$sc_experiment <- sc_so@meta.data$exp[1]
    self$bulk_experiment <- bulk_so@meta.data$exp[1]

    self$bulk_so <- bulk_so
    self$sc_so <- sc_so

    ## Informative features
    self$i_feats <- find_shared_genes(self$sc_so, self$bulk_so)

    prep_data <- precise_preproc(
      sc_so = self$sc_so,
      bulk_so = self$bulk_so,
    )
    self$bulk_so <- prep_data[['bulk_so']]
    self$sc_so <- prep_data[['sc_so']]

    ## Compute mapping based on pseudobulked data
    self$CFs <- tryCatch(gen_CF_fitter(
      t(self$get_bulk_M()),
      t(self$get_sc_M()),
      N_source = self$N_PV_bulk,
      N_target = self$N_PV_sc),
      error = function(e) { browser() }
    )

    if (is.null(self$CFs)) return(NULL)

    self$loadings_M <-
      self$CFs$consensus_representation %>%
      as.matrix %>%
      set_rownames(rownames(self$get_bulk_M())) %>%
      { set_colnames(., paste0('CF', 1:ncol(.))) }
    # vm <- stats::varimax(loadings_M)
    # loadings(vm)[] %>%

    self$bulk_var <- var_explained_CF(self$loadings_M,
      self$get_bulk_M())
    # self$sc_var <- var_explained_CF(self$loadings_M,
    #   self$get_sc_M())
    self$sc_var <- var_explained_CF(self$loadings_M,
      self$get_sc_M())
    self$compute_projections()
    # self$get_optimal_gamma(redo = F)
  },
  project = function(M) {
    M %>% t() %>% { . %*% self$loadings_M }
  },
  reconstruct = function(M) {
    M %>% t() %>% { . %*% self$loadings_M }
  },
  get_bulk_so = function() {
    self$bulk_so[self$i_feats, ]
  },
  get_bulk_M = function() {
    so2M(self$get_bulk_so())
  },
  get_sc_so = function() {
    self$sc_so[self$i_feats, ]
  },
  get_sc_M = function() {
    so2M(self$get_sc_so())
  },
  get_scpb_so = function() {
    self$scpb_so[self$i_feats, ]
  },
  get_scpb_M = function() {
    so2M(self$get_scpb_so())
  },
  compute_projections = function() {

    bsa <- self[['bulk_so']]@meta.data %>%
      dplyr::select(-any_of('sample_name')) %>%
      tibble::rownames_to_column('sample_name') %>%
      { . }
    self$bulk_P_a <-
      self$project(self$get_bulk_M()) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('sample_name') %>%
      maartenutils::controlled_merge(bsa, by_cols = 'sample_name') %>%
      tibble::as_tibble() %>%
      # pivot_longer(cols = matches('_')) %>%
      numerify_regressors() %>%
      # dplyr::mutate(sn_dilution =
      #   as.numeric(as.character(sn_dilution))) %>%
      # dplyr::mutate(sn_dilution =
      #   if_else(sn_dilution == 0, 0, 1/sn_dilution)) %>%
      # dplyr::mutate(
      #   condition_name = glue::glue('{stim_group} - {duration}h')) %>%
      { . }

    ssa <- self[['sc_so']]@meta.data %>%
      dplyr::select(-any_of('sample_name')) %>%
      tibble::rownames_to_column('sample_name') %>%
      { . }
    self$sc_P_a <-
      self$project(self$get_sc_M()) %>%
      as.data.frame() %>%
      tibble::rownames_to_column('sample_name') %>%
      maartenutils::controlled_merge(ssa, by_cols = 'sample_name') %>%
      tibble::as_tibble() %>%
      # pivot_longer(cols = matches('_')) %>%
      numerify_regressors() %>%
      # dplyr::mutate(sn_dilution =
      #   as.numeric(as.character(sn_dilution))) %>%
      # dplyr::mutate(sn_dilution =
      #   if_else(sn_dilution == 0, 0, 1/sn_dilution)) %>%
      # dplyr::mutate(
      #   condition_name = glue::glue('{stim_group} - {duration}h')) %>%
      { . }

    # if (is.null(self$pseudobulk) || self$pseudobulk == 'none') {
    #   if (F) {
    #    sc_P_s <- Seurat::AverageExpression(
    #       sc_so,
    #       assays = 'PRECISE',
    #       group.by = 'condition_name')[[1]] %>%
    #       t()
    #   } else {
    #     sc_P_s <- t(apply(sc_P, 2, function(x)
    #         tapply(x, sc_so@meta.data$condition_name, median)))
    #   }
    # } else {
    #   sc_P_s <- sc_P
    # }
    # ssa <- self[['sc_so']]@meta.data
    # if (!'sample_name' %in% colnames(ssa)) {
    #   ssa <- tibble::rownames_to_column(ssa, 'sample_name')
    # }
    # self$sc_P_a <-
    #   self$project(self$get_sc_M()) %>%
    #   as.data.frame() %>%
    #   tibble::rownames_to_column('sample_name') %>%
    #   maartenutils::controlled_merge(ssa, by_cols = 'sample_name') %>%
    #   numerify_regressors() %>%
    #   tibble::as_tibble() %>%
    #   { . }

    shared_cols <-
      colnames(self$sc_P_a) %>%
      intersect(colnames(self$bulk_P_a)) %>%
      c('duration', 'tnf_conc', 'ifn_conc')
    self$comb_P <- bind_rows(
      self$bulk_P_a[, intersect(colnames(self$bulk_P_a), shared_cols)],
      self$sc_P_a[, intersect(colnames(self$sc_P_a), shared_cols)]) %>%
      dplyr::select(-any_of('condition_name'), everything())
    if (!'duration' %in% colnames(self$comb_P))
      stop('No duration in comb_P')

    # self$UMAP_emb <-
    #   self$comb_P %>%
    #   dplyr::select(matches('CF\\d+')) %>%
    #   as.matrix %>%
    #   umap::umap() %>%
    #   purrr::pluck('layout') %>%
    #   set_colnames(c('UMAP1', 'UMAP2')) %>%
    #   as.data.frame() %>%
    #   { . }
    # self$comb_P <- dplyr::bind_cols(self$comb_P, self$UMAP_emb)
  },
  add_function = function(name, meth) {
    self[[name]] <- meth
    environment(self[[name]]) <- environment(self$add_function)
  },
  compute_feature_dist = function() {
    CFW <- diag(rep(1, self$N_CFs))
    sel_CF <- function(x) {
      # rn <- rownames(x)
      out <- dplyr::select(x, matches('CF\\d*')) %>%
        as.data.frame()
      # rownames(out) <- rn
      rownames(out) <- x$sample_name
      out
    }
    FM <- rbind(
      self[['bulk_P_a']] %>% sel_CF,
      self[['scpb_P_a']] %>% sel_CF
    )
    ## Substances to distances from query (single cell) to reference
    ## (bulk) samples
    self$distM <-
      as.matrix(dist(as.matrix(FM) %*% CFW)) %>% {
        r_idxs <- (nrow(self[['bulk_P_a']])+1):nrow(.)
        c_idxs <- 1:nrow(self[['bulk_P_a']])
        .[r_idxs, c_idxs]
      } %>%
      { 1e3 * . / max(., na.rm = T) }
    return(self$distM)
  },
  compute_sim_scores = function(gamma_v) {
    self$compute_feature_dist()
    distM_w <- exp_similarity(distM = self$distM, gamma_v = gamma_v)
    if (nrow(distM_w) == 0) return(NULL)
    ## Predict scores before dropping unpopular reference samples
    pred_scores <- score_dist_M(
      distM_w,
      self[['bulk_P_a']],
      self[['sc_P_a']],
      grading_vars = get_grading_vars(self$sc_experiments)
    )
    return(invisible(pred_scores))
  },
  gamma_fitness = function(gamma_v = .1) {
    if (!is.finite(abs(gamma_v))) {
      return(list('Score' = 0, 'Pred' = 0))
    }
    scores <- self$compute_sim_scores(gamma_v = gamma_v)
    if (is.null(scores) || nrow(scores) == 0) {
      return(list('Score' = 0, 'Pred' = 0))
    }
    mean_scores <- purrr::map_dbl(auto_name(int_data$grading_vars),
        ~mean(scores[[.x]]*scores[['confidence']], na.rm = T))
    list(
      'Score' = nrow(scores) * (1 - mean(mean_scores)),
      'Pred' = 0
    )
  },
  get_optimal_gamma = function(redo = F) {
    if (!is.null(self$optimal_gamma) && !redo)
      return(self$optimal_gamma)
    self$gamma_grid <- tibble::tibble(
      gamma = 10^seq(-2, 1, by = 1e-2)) %>%
      dplyr::mutate(score = purrr::map_dbl(gamma,
        function(x) {
          tryCatch(self$gamma_fitness(x)[['Score']],
            error = function(e) { 0 })
        }))
    if (all(maartenutils::eps(self$gamma_grid$score, 0))) {
      rlang::warn('All scores are zero')
    }
    self$optimal_gamma <-
      self$gamma_grid$gamma[which.max(self$gamma_grid$score)]
    if (self$optimal_gamma %in%
      unlist(self$gamma_grid[c(1, nrow(self$gamma_grid)), 'gamma'])) {
      rlang::warn('Gamma on edge of parm scan grid')
    }
    return(self$optimal_gamma)
  })
)


update_methods <- function(obj) {
  class_name <- class(obj)
  empty_init_args <-
    map_lgl(formals(obj$initialize), ~.x == '') %>%
    which %>% names

  new_obj <- purrr::exec(get(class_name)$new,
    as.list(obj)[empty_init_args])
  functions <- purrr::map_lgl(as.list(new_obj),
    ~'function' %in% class(.x)) %>%
    which %>% names
  for (f in functions) {
    obj[[f]] <- new_obj[[f]]
  }
  return(obj)
}
# update_methods(prep_data)


AddPreciseEmbedding <- function(
  so,
  ref_M = targets::tar_read('kallisto_5029'),
  genes = rownames(targets::tar_read('kallisto_5029')),
  N_source = 10L, N_target = 20L) {

  genes <- rownames(so@assays$SCT[,]) %>%
    intersect(rownames(ref_M)) %>%
    intersect(genes) %>%
    { .  }

  sc_M <- t(as.matrix(so@assays$SCT[genes,]))
  ref_M <- t(as.matrix(ref_M[genes,]))

  PVs <- gen_fitter()$fit(ref_M, sc_M)

  M <- {sc_M  %*% t(PVs$source_components_) } %>%
    t() %>%
    { set_rownames(., paste0('PV', 1:nrow(.))) }

  # D <- compute_RE(ref_M, PVs, T)
  # D <- compute_RE(sc_M, PVs, T)
  # apply(D, 1, median)
  so[['PRECISE']] <- CreateAssayObject(data = M)

  return(so)
}


plot_var_exp <- function(dtf) {
  library(ggplot2)

  p_dat <- dtf %>%
    dplyr::mutate(source = cumsum(source)) %>%
    dplyr::mutate(target = cumsum(target)) %>%
    dplyr::mutate(
      pv_idx = as.integer(gsub('PV(\\d)', '\\1', pv_idx))
    ) %>%
    tidyr::pivot_longer(cols = c(source, target)) %>%
    { . }

  ggplot(data = p_dat,
    aes(x = pv_idx, y = 100 * value, colour = name)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(name = 'PV number',
      breaks = seq(1, nrow(dtf), by = 2)) +
    scale_y_continuous(name = '% variance explained') +
    scale_colour_discrete(name = '') +
    theme(
      legend.position = c(.05, .95),
      legend.direction = 'vertical',
      legend.justification = c(0,1))
}


# var_explained_CF <- function(CFs, X) {
#   total_var <- sum(apply(X, 2, var))
#   proj_X <- X %*% CFs
#   var_proj_X = apply(proj_X, 2, var) / total_var
#   return(tibble(
#       'cf_idx' = as.character(glue::glue('CF{seq_along(var_proj_X)}')),
#       'var' = var_proj_X))
# }
#' X = [features x samples], Seurat format
var_explained_CF <- function(CFs, X) {
  if (is.null(X) || is.null(CFs)) return(NULL)
  total_var <- sum(apply(X, 1, var))
  proj_X <- t(X) %*% CFs
  var_proj_X = apply(proj_X, 2, var) / total_var
  return(
    tibble::tibble(
      'cf_idx' = as.character(glue::glue('CF{seq_along(var_proj_X)}')),
      'var' = var_proj_X
    )
  )
}


samplewise_exp_var <- function(M, loadings_M) {
  # norm(M, type = 'F')
  ## Projected = L %*% (M.t %*% L).t
  M_s <- t(M) %*% loadings_M
  M_p <- loadings_M %*% t(M_s)
  vec_length <- function(x) { sqrt(sum(x^2)) }
  apply(M_p, 2, vec_length) / apply(M, 2, vec_length) %>%
    tibble::enframe('sample_name', 'exp_var')
}


TPM_M <- function(M, equal_genes = T) {
  if (equal_genes) {
    M <- apply(M, 2, function(x) x / mean(x))
    stopifnot(
      maartenutils::eps(apply(M, 2, function(x) mean(x)), 1, 1e-6))
  }
  M <- t(apply(M, 1, function(x) 1e6 * x / sum(abs(x))))
  return(M)
}


#' Preprocess two Seurat objects in preparation of subsequent
#' integration
#'
#'
precise_preproc <- function(
  bulk_so, sc_so,
  assay = DefaultAssay(bulk_so),
  datatype = 'scale.data') {
  library(Seurat)

  scn <- find_shared_genes(bulk_so, sc_so)
  if (length(scn) == 0) stop('No shared gene names')
  bulk_so <- subset_feats(bulk_so, scn)
  sc_so <- subset_feats(sc_so, scn)

  bulk_M <- so2M(bulk_so, datatype = datatype)
  sc_M <- so2M(sc_so, datatype = datatype)

  bulk_so <- SetAssayData(bulk_so, slot = datatype,
    new.data = bulk_M, assay = assay)
  sc_so <- SetAssayData(sc_so, slot = datatype,
    new.data = sc_M, assay = assay)

  return(list('bulk_so' = bulk_so, 'sc_so' = sc_so))
}


plot_CF_loadings_HM <- function(
  bulk_M, sc_M,
  N_genes = 30, CF_idxs = 1:5,
  loadings_M = CFs$consensus_representation) {

  bulk_var <- var_explained_CF(loadings_M, bulk_M)
  # bulk_var <- var_explained_CF(loadings_M, bulk_M)
  sc_var <- var_explained_CF(loadings_M, sc_M)

  loadings_M <- loadings_M %>%
    # { t(t(.) %*% diag(rowMeans(bulk_M))) } %>%
    # { . %*% diag(sqrt(sc_var$var)) } %>%
    set_rownames(rownames(bulk_M)) %>%
    set_colnames(paste0('CF', 1:ncol(.)))
  loadings_M <- loadings_M[, CF_idxs, drop=F]

  if (length(CF_idxs) > 1) {
    ## In choosing the most weighted genes, give CFs with high amounts
    ## of explained variance priority
    max_loadings_per_gene <-
      loadings_M %>%
        { . %*% diag(sqrt(sc_var$var[CF_idxs])) } %>%
        apply(1, function(x) max(abs(x)))
  } else {
    ## This is pointless now
    max_loadings_per_gene <-
      loadings_M %>%
        { . * sqrt(sc_var$var[CF_idxs]) } %>%
        apply(1, function(x) max(abs(x)))
  }

  if (F && !test_rendering()) {
    pdf(file.path(img_dir, glue::glue('max_loading_ecdf.pdf')),
      height = 5)
    plot(ecdf(max_loadings_per_gene))
    dev.off()
  }

  sel_genes <-
    which(rank(max_loadings_per_gene) >=
    (length(max_loadings_per_gene) - N_genes + 1)) %>%
    names
  loadings_M_sel <- loadings_M[sel_genes, , drop = F]

  gs_data <- targets::tar_read(gs_data_step)
  ra <- tibble(gene = sel_genes) %>%
    left_join(gs_data, by = 'gene') %>%
    dplyr::select(geneset, TNFa_bias, synergy_score, time_class,
      max_timepoint) %>%
    dplyr::mutate(synergy_score = sign(synergy_score) *
      abs(synergy_score)^1/4)
  HM <- gen_HM(loadings_M_sel,
    value_name = 'Gene loading',
    ra = ra,
    sa = tibble(
      'var_ref' = bulk_var$var[CF_idxs],
      'sc_var' = sc_var$var[CF_idxs]),
    cluster_columns = F,
    column_annotation_name_side = 'right',
    trans = NULL,
    show_column_names = TRUE,
    show_row_dend = F,
    column_title = 'Consensus factor\nloadings',
    width = unit(.75 * length(CF_idxs), 'cm'),
    show_row_names = TRUE,
    N_hl_genes = NULL,
    N_genes = NULL)

  M_GE <- bulk_M[sel_genes, ] %>%
    as.matrix()
  lFC_M_GE <- M_GE  %>%
    { log2(. + 1) } %>%
    apply(1, function(x) x - x['2h_unstim']) %>%
    t() %>%
    { .[, -which(abs(colMeans(.)) <= .5)] }
  # M_GE <- 2^M_GE
  # hist(abs(colMeans(M_GE)), breaks = 100)
  # log2(targets::tar_read('kallisto_5029')['IL32', '2h_unstim'] + 1)
  colnames(lFC_M_GE) <-
    tibble(sample_name = colnames(lFC_M_GE)) %>%
      left_join(sample_annotation, by = 'sample_name') %>%
      pull(condition_name)
  HM_GE <- gen_HM(
    M = lFC_M_GE,
    sa = tibble(condition_name = colnames(lFC_M_GE)) %>%
      left_join(sample_annotation, by = 'condition_name'),
    N_genes = NULL,
    N_hl_genes = NULL,
    column_annotation_name_side = 'left',
    show_column_dend = T,
    show_row_names = T,
    row_names_side = 'left',
    show_row_dend = T,
    width = unit(.4 * ncol(lFC_M_GE), 'cm'),
    cluster_rows = F,
    show_column_names = T,
    column_title = 'Gene expression',
    value_name = 'log2FC gene expression\nvs. 2h unstimulated'
  )

  code_block <- rlang::expr(
    draw(HM_GE[row_order(HM)] + HM[row_order(HM)],
      merge_legend = T, newpage = F,
      row_title = 'Genes',
      heatmap_legend_side = 'bottom',
      column_title = glue::glue('CF {paste(CF_idxs, collapse = \', \')}')
    )
  )
  if (!test_rendering()) {
    pdf(file.path(img_dir, glue::glue('loading_heatmap.pdf')), height = 10)
    eval(code_block)
    dev.off()
  } else {
    eval(code_block)
  }
}


numerify_ranks <- add_rank_vars <- function(dtf) {
  if (!is.null(dtf$sn_dilution)) {
    dtf$sn_dilution <- as.numeric(as.character(dtf$sn_dilution))
  }
  if (!is.null(dtf$duration)) {
    dtf$duration <- as.numeric(as.character(dtf$duration))
  }
  if (!is.null(dtf$tnf_rank)) {
    dtf$tnf_rank <- as.numeric(as.character(dtf$tnf_rank))
  }
  if (!is.null(dtf$ifn_rank)) {
    dtf$ifn_rank <- as.numeric(as.character(dtf$ifn_rank))
  }
  if (!is.null(dtf$duration) && !is.numeric(dtf$duration)) {
    browser()
  }
  dtf
}


#' Score a similarity matrix, answering the question how close
#' predicted classes are to expected classes
#'
#'
score_dist_M <- function(
  distM_w,
  ref_ann,
  query_ann,
  query_sample_id = 'sample_name',
  v_mod_ref = function(gv) gv,
  v_mod_query = function(gv) glue::glue('{gv}_med'),
  imp_mod_query = function(gv) glue::glue('{gv}_certainty'),
  # distance_fun = function(x, y) sqrt((x-y)^2),
  distance_fun = function(x, y) abs(x-y)/y,
  grading_vars = c('duration', 'sn_dilution',
    'tnf_conc', 'ifn_conc')) {

  stopifnot(is.character(grading_vars))
  grading_vars <- grading_vars %>%
    intersect(colnames(ref_ann)) %>%
    intersect(colnames(query_ann))
  if (length(grading_vars) == 0) return(NULL)

  if (F) {
    ## Put the annotation in the same order as the rows and columns of
    ## the similarity matrix
    ref_ann <- recover_ann(query = colnames(distM_w),
      lookup_data = ref_ann) %>%
      add_rank_vars()
    query_ann <- recover_ann(query = rownames(distM_w),
      lookup_data = query_ann) %>%
      add_rank_vars()
  }

  # gv <- grading_vars[2]
  # browser()
  stim_cc
  distM_w
  grade_var_scores <- map(grading_vars, function(gv) {
    query_obs <- numerify_factor(query_ann[[v_mod_query(gv)]])
    ref_obs <- numerify_factor(ref_ann[[gv]])
    DM <- outer(query_obs, ref_obs, FUN = distance_fun) %>%
      { (. + 1) / (max(., na.rm = T) + 1) } %>%
      # set_rownames(query_ann[[query_sample_id]]) %>%
      # set_colnames(ref_ann[[query_sample_id]])
      { . }
    { distM_w * DM } %>%
      stats::na.omit() %>%
      apply(1, sum, na.rm = T) %>%
      tibble::enframe(query_sample_id, gv) %>%
      dplyr::mutate(confidence = query_ann[[imp_mod_query(gv)]]) %>%
      # dplyr::mutate({{query_sample_id}} :=
      #   query_ann[[query_sample_id]])
      { . }
  })

  if (F) {
    out <- grade_var_scores %>%
      purrr::reduce(inner_join, by = query_sample_id) %>%
      dplyr::select(-matches('confidence\\.'))
  } else {
    first_occurence <- glue::glue('{query_sample_id}...1')
    out <- suppressMessages(map_dfc(grade_var_scores, ~.x)) %>%
      dplyr::mutate({{query_sample_id}} := .data[[first_occurence]]) %>%
      dplyr::select(-matches(glue::glue('{query_sample_id}\\.\\.\\..*'))) %>%
      { . }
  }

  return(out)
}


plot_scree <- function(...) UseMethod('plot_scree')


#' Do a partial PCA and return the Scree plot
#'
#' Assumes features/variables to be in the rows, samples in the
#' columns
#'
#'
plot_scree.matrix <- function(M, N_pcs = 10L, center = T, scale = T) {
  N_pcs <- min(N_pcs, ncol(M)-1L)
  stopifnot(is.integer(N_pcs))
  library(irlba)

  pc <- irlba::prcomp_irlba(
    M, 
    # n = min(ncol(M)-1), 
    n = min(nrow(M)-1L, ncol(M)-1L, N_pcs), 
    # n = 10L,
    center = center, 
    scale = scale, 
    maxit = 1e6
  )
  var_explained <- pc$sdev^2 / pc$totalvar
  var_explained <- var_explained[1:N_pcs]

  p_dat <- tibble::tibble(
    # idx = factor(seq_along(var_explained)),
    idx = as.integer(seq_along(var_explained)),
    VE = cumsum(var_explained)
  )

  p1 <-
    ggplot2::ggplot(p_dat, aes(x = idx, y = VE)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_continuous(
        name = 'Principal component',
        # breaks = 1:max(p_dat$idx),
        expand = c(0.01, 0.0)) +
      ggplot2::scale_y_continuous(name = 'Variance explained',
        expand = c(0.01, 0.0))

  return(p1)
}


plot_scree.data.frame <- function(dtf, N_pcs = 15L) {
  plot_scree(as.matrix(dtf), N_pcs = N_pcs)
}


plot_scree.Seurat <- function(
  so, title = Project(so), N_pcs = 15L, return_plot = F) {
  M <- so2M(so)
  p <- plot_scree(M, N_pcs = as.integer(min(ncol(M)-1, N_pcs))) +
    ggtitle(glue::glue('{title} scree'))
  if (return_plot) {
    return(p)
  } else {
    magic <- 'MAGIC' %in% names(so@assays)
    o_fn <- file.path(img_dir,
      glue::glue('PCA_scree-{Project(so)}{make_flag(magic)}.png'))
    maartenutils::print_plot(p, w = 8.7, h = 6, fn = o_fn)
    return(o_fn)
  }
}


plot_CF_scores <- function(
  int_data, x_var = 'CF1', y_var = 'CF2', colour_var = 'stim_group') {
  p_clust_bar <- gen_sc_cluster_comp_plot(int_data$sc_so)

  # int_data$bulk_so@meta.data
  bulk_dat <- int_data$comb_P %>%
    dplyr::filter(sample_type == 'bulk')
  if (is.null(bulk_dat$duration)) {
    stop('No duration dat')
    sc_dat$duration
  }

  sc_dat <- int_data$comb_P %>% dplyr::filter(sample_type == 'sc')

  cluster_names <- levels(p_clust_bar$data$seurat_clusters)
  l1 <-
    naturalsort::naturalsort(unique(as.character(sc_dat$sample_name)))
  # l2 <-
  #   unique(as.character(sort(int_data$sc_so@meta.data$seurat_clusters)))
  # if (all(l1 != l2)) {
  #   stop('ldllldf')
  # }
  if (max(as.numeric(l1))-1 == max(as.numeric(cluster_names))) {
    sc_dat$stim_group <-
      as.character(as.numeric(sc_dat$stim_group) - 1)
  }

  shape_values <- levels(droplevels(bulk_dat$duration)) %>%
    { set_names(14+1:length(.), .) }
  N_clusters <- length(cluster_names)
  cluster_cols <- N_clusters %>%
    { set_names(viridis::magma(., begin = .2, end = 1),
      cluster_names) }

  p1 <- bulk_dat %>%
    ggplot(aes_string(x = x_var, y = y_var,
        label = 'sample_name',
        colour = colour_var, shape = 'duration')) +
    geom_point(alpha = .5, size = 3.5) +
    geom_point(size = 1) +
    scale_x_continuous(limits = range(int_data$comb_P[[x_var]])) +
    scale_y_continuous(limits = range(int_data$comb_P[[y_var]])) +
    scale_shape_manual(name = 'Exposure duration',
      values = shape_values) +
    scale_colour_manual(
      name = 'Reference condition',
      values = gen_cyto_inf_col_scale(int_data$bulk_P_a$stim_group)) +
    guides(color = guide_legend(
        # nrow = max(length(cols) / 4, 2)
        ncol = 1L
        )) +
    theme(legend.position = 'bottom') +
    ggtitle('Reference')
    # ggrepel::geom_text_repel(max.overlaps = 100L)

  p2 <- sc_dat %>%
    ggplot(aes_string(x = x_var, y = y_var,
        shape = 'duration',
        label = 'stim_group', colour = colour_var)) +
    geom_point(data = bulk_dat,
      alpha = .5, size = 4, colour = 'grey90', guide = 'none') +
    geom_point(alpha = .5, size = 12, shape = 20) +
    geom_text(guide = 'none', colour = 'black') +
    scale_shape_manual(values = shape_values) +
    scale_x_continuous(limits = range(int_data$comb_P[[x_var]])) +
    scale_y_continuous(limits = range(int_data$comb_P[[y_var]])) +
    scale_colour_manual(
      name = 'Cluster',
      values = cluster_cols) +
    guides(color = 'none', shape = 'none') +
    ggtitle('Pseudobulked single cell')
    # ggrepel::geom_text_repel(max.overlaps = 100L)

  cluster_ann <- tibble(
      # x = seq(.5, N_clusters-.5, length.out = N_clusters), y = 1,
      x = seq(1, N_clusters), y = 1,
      label = cluster_names
    ) %>%
    ggplot(aes(x = x, y = 1, label = label, colour = label,
        fill = label)) +
      geom_point(size = 5, alpha = .5) +
      scale_x_continuous(expand = c(0, 0),
        limits = c(.5, N_clusters+.5)) +
      geom_text(colour = 'black', hjust = .5, size = 3) +
      scale_fill_manual(values = cluster_cols, name = '') +
      scale_colour_manual(values = cluster_cols, name = '') +
      gg_tabula_rasa +
      theme(legend.position = 'none')
  if (F) {
    print_plot(p_clust_bar,
      fn = file.path(img_dir,
        glue::glue('condition_by_clusters{int_data$sc_experiment}.png')),
      w = 8.7, h = 10)
  }
  # p2_dims <- get_dim(p2)
  top0 <- theme_cyto_inf()[['plot.margin']] %>%
    { .[1] <- unit(0, 'cm'); . }
  bottom0 <- theme_cyto_inf()[['plot.margin']] %>%
    { .[3] <- unit(0, 'cm'); . }
  p_clust_bar_ann <-
    # (set_dim(cluster_ann, p2_dims) +
    #   ggplot2::theme(plot.margin = top0)) /
    # (set_dim(p_clust_bar, p2_dims) +
    #   ggplot2::theme(plot.margin = bottom0)) +
    # set_dim(cluster_ann, p2_dims) /
    # set_dim(p_clust_bar, p2_dims) +
    cluster_ann / p_clust_bar +
    plot_layout(heights = c(1, 6.5))
  # p_clust_bar_ann <- gridExtra::arrangeGrob(
  #   (cluster_ann + ggplot2::theme(plot.margin = bottom0)),
  #   (p_clust_bar + ggplot2::theme(plot.margin = top0)),
  #   heights = c(1, 6.5))
  # print_plot(p_clust_bar_ann, file.path(img_dir, 'test.png'),
  #   w = 8.7)

  (p1 + p2) / (plot_spacer() + p_clust_bar_ann) +
    # plot_layout(heights = c(13, 7), guides = 'collect')
    plot_layout(heights = c(13, 7))
}


plot_cosine_sim <- function(M_source, M_target, N_source, N_target) {
  PVs <- gen_fitter(N_source = N_source, N_target = N_target)
  PVs <- PVs$fit(t(M_source), t(M_target))

  permute_bulk_cs <- function(i) {
    M_source_p <- M_source[base::sample(1:nrow(M_source)), ]
    dimnames(M_source_p) <- dimnames(M_source)
    PVs_p <- gen_fitter(
      N_source = N_source,
      N_target = N_target
    )
    PVs_p <- PVs_p$fit(t(M_source_p), t(M_target))
    tibble(
      cs = diag(PVs_p$cosine_similarity_matrix_)) %>%
      dplyr::mutate(idx = 1:n()) %>%
      dplyr::mutate(i = i)
  }
  permuted <- map_dfr(1:N_source, permute_bulk_cs)

  CS_dtf <- permuted %>%
    dplyr::group_by(idx) %>%
    dplyr::summarize(CI_l = quantile(cs, .1),
                     CI_h = quantile(cs, .9)) %>%
    dplyr::right_join(
      tibble(
        idx = 1:ncol(PVs$cosine_similarity_matrix_),
        obs_CS = diag(PVs$cosine_similarity_matrix_)),
      by = 'idx')

  p <- ggplot(CS_dtf,
    aes(x = idx, y = obs_CS, ymin = CI_l, ymax = CI_h)) +
    geom_ribbon(alpha = .1) +
    geom_line() +
    geom_point() +
    scale_x_continuous(name = 'PV number', breaks = CS_dtf$idx,
      expand = c(0, 0)) +
    scale_y_continuous(name = 'Cosine similarity', expand = c(0, 0))
  return(p)

}


plot_cosine_heatmap <- function(
  M_source, M_target, N_source, N_target,
  o_fn = file.path(Sys.getenv('img_dir'),
    glue::glue('cosine_heatmap.png'))) {

  PVs <- gen_fitter(N_source = N_source, N_target = N_target)
  PVs <- PVs$fit(t(M_source), t(M_target))

  HM <- PVs$initial_cosine_similarity_matrix_ %>%
    { set_colnames(., 1:ncol(.)) } %>%
    { set_rownames(., 1:nrow(.)) } %>%
    Heatmap(
      name = 'Cosine similarity',
      cluster_rows = F,
      cluster_columns = F,
      show_row_dend = F,
      show_column_dend = F
    )
  print_plot_eval(
    { draw(HM) },
    width = 17.4, height = 15,
    filename = o_fn
  )

  return(o_fn)
}


plot_var_exp_panel <- function(
  M_source, M_target, N_source, N_target,
  o_fn = file.path(Sys.getenv('img_dir'),
    glue::glue('precise_var_exp.png'))) {

  PVs <- gen_fitter(N_source = N_source, N_target = N_target)
  PVs <- PVs$fit(t(M_source), t(M_target))

  var_source <-
    var_explained(PVs, t(M_source)) %>%
    dplyr::mutate(
      p_title = glue('{pv_idx} - {round(100*source, 1)}%')) %>%
    { . }

  var_target <-
    var_explained(PVs, t(M_target)) %>%
    dplyr::mutate(
      p_title = glue('{pv_idx} - {round(100*source, 1)}%')) %>%
    { . }

  p1 <- plot_var_exp(var_source) +
    ggtitle('Reference/bulk variance')

  p2 <- plot_var_exp(var_target) +
    ggtitle('Target/single cell variance')

  print_plot_eval(print(p1 + p2), width = 12, filename = o_fn)
}


plot_var_exp_CF_panel <- function(
  M_source, M_target, N_source, N_target,
  o_fn = file.path(Sys.getenv('img_dir'),
    glue::glue('precise_var_exp.png'))) {
  browser()

  CFs <- gen_CF_fitter(
    N_source = N_source, N_target = N_target,
    X_source = M_source, X_target = M_target
  )

  var_source <-
    var_explained(CFs, t(M_source)) %>%
    dplyr::mutate(
      p_title = glue('{pv_idx} - {round(100*source, 1)}%')) %>%
    { . }

  var_target <-
    var_explained(CFs, t(M_target)) %>%
    dplyr::mutate(
      p_title = glue('{pv_idx} - {round(100*source, 1)}%')) %>%
    { . }

  p1 <- plot_var_exp(var_source) +
    ggtitle('Reference/bulk variance')

  p2 <- plot_var_exp(var_target) +
    ggtitle('Target/single cell variance')

  print_plot_eval(print(p1 + p2), width = 12, filename = o_fn)
}


gen_PRECISE_wrapper <- function(reference_M, query_M) {
  tmp_dir <- stringr::str_replace(tempfile(), 'file', 'dir')
  dir.create(tmp_dir)

  force(reference_M)
  force(query_M)
  readr::write_csv(reference_M, file.path(tmp_dir, 'reference.csv'))
  readr::write_csv(query_M, file.path(tmp_dir, 'query.csv'))
   
  PRECISE_caller <- function(
    PRECISE_params = list(
      N_source_components = 10L, 
      N_target_components = 10L,
      N_PV = 5L,
      mean_center = FALSE,
      std_unit = FALSE
    ), 
    verbose = F, 
    ncores = 1, 
    perform_validity_tests = TRUE) {

    PRECISE_params$N_source_components <- with(PRECISE_params,
      min(N_source_components, nrow(reference_M)-1))
    PRECISE_params$N_target_components <- with(PRECISE_params, 
      min(N_target_components, nrow(query_M)-1))
    PRECISE_params$N_PV <- with(PRECISE_params, 
      min(N_PV, N_source_components, N_target_components))

    id_string <- PRECISE_params %>%
      modifyList(c(list(ps = '_'))) %>%
      with(paste(
        N_source_components,
        prepend_string(N_target_components, ps),
        prepend_string(N_PV, ps),
        prepend_string(mean_center, ps),
        prepend_string(std_unit, ps),
        sep = '')
      )

    py_file <- file.path(Sys.getenv('python_dir'), 'run_precise.py')
    stopifnot(file.exists(py_file))
    command <- glue::glue('python {py_file} \\
      {tmp_dir} \\
      {PRECISE_params$N_source_components} \\
      {PRECISE_params$N_target_components} \\
      {PRECISE_params$N_PV} \\
      {PRECISE_params$mean_center} \\
      {PRECISE_params$std_unit} \\
      {id_string}')

    system(command, wait = T)
    # print(list.files(tmp_dir))

    tryCatch({
      out <- c('reference_CF', 'query_CF', 'CS_mat', 
        'consensus_representation') %>%
        maartenutils::auto_name() %>%
        purrr::map(function(on) {
          readr::read_csv(
            file.path(tmp_dir, glue::glue('{on}_{id_string}.csv')), 
            show_col_types = FALSE
          )
        })
      out$reference_CF <- as.matrix(out$reference_CF)
      out$query_CF <- as.matrix(out$query_CF)
      out$consensus_representation <-
        as.matrix(out$consensus_representation)
      colnames(out$reference_CF) <- colnames(out$query_CF) <-
        colnames(out$consensus_representation) <-
        glue::glue('CF{1:ncol(out$reference_CF)}')
      rownames(out$reference_CF) <- rownames(reference_M)
      rownames(out$query_CF) <- rownames(query_M)
      rownames(out$consensus_representation) <-
        colnames(reference_M)
      if (F && perform_validity_tests) {
        cs <- svd(out$CS_mat)$d
        out$valid_cs <- all(abs(cs) <= 1)
      }
    }, error = function(e) { print(e); NULL }) 

    return(out)
  }

  return(PRECISE_caller)
}
