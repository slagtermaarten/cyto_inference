exclude_non_distinguishable_stimuli <- function(...)
  UseMethod('exclude_non_distinguishable_stimuli')


exclude_non_distinguishable_stimuli.data.frame <- function(dtf) {
  stopifnot('stim_group' %in% colnames(dtf))

  dtf <- dtf %>%
    dplyr::filter(!stim_group %in% c('1/20000 SN', '0.01 ng/ml IFNy'))

  return(dtf)
}


exclude_non_distinguishable_stimuli.Seurat <- function(so) {
  dtf <- so@meta.data
  dtf$i <- 1:nrow(dtf)

  dtf <- exclude_non_distinguishable_stimuli(dtf)
  return(so[, dtf$i])
}


#' As unstimulated duration cannot be estimated (eyeballing the data),
#' let's not penalize models for getting it wrong. The simplest way of
#' achieving this is by setting 'duration' vars to NA for unstimulated
#' conditions
#'
#' DEPRECATED with the advent of separate TNF and IFNy exposure
#' duration parameters. A zero value for these new parameters can now
#' be interpreted as having had no exposure.
#'
remove_unidentifiable_fields <- function(dtf) {
  if (maartenutils::null_dat(dtf)) return(dtf)
  stopifnot(is.data.frame(dtf))
  stopifnot('stim_group' %in% colnames(dtf))

  dtf <- numerify_regressors(dtf)

  for (vn in c('tnf_duration', 'ifn_duration', 'duration')) {
    if (vn %in% colnames(dtf)) {
      dtf[[vn]][grepl('Unstimulated', dtf$stim_group)] <- NA_real_
    }
  }

  return(dtf)
}


raw_code <- rlang::expr({
  so_P
})


harmony_code <- rlang::expr({
  library(harmony)
  RunHarmony(
    so_P,
    group.by.vars = c('exp'),
    assay = DefaultAssay(so_P),
    max.iter.harmony = N_pcs
  )
})


seurat_rpca_code <- rlang::expr({
  stopifnot(!is.null(reference) && reference %in% experiments)
  var_feats <- VariableFeatures(so_P)
  sample_list <- SplitObject(so_P, split.by = 'exp')
  rm(so_P); gc()

  max_dims <- min(min(purrr::map_int(sample_list, ncol))-1, 30)

  # sample_list <- sample_list[names(sample_list) != '6434']
  sample_list <- PrepSCTIntegration(sample_list)

  # features <- SelectIntegrationFeatures(object.list = sample_list)
  # sample_list <- lapply(sample_list, FUN = function(x) {
  #     x <- ScaleData(x, features = features, verbose = FALSE)
  #     RunPCA(x, features = features, verbose = FALSE,
  #       dims = 1:max_dims)
  # })

  anchors <- Seurat::FindIntegrationAnchors(
    object.list = sample_list,
    # object.list = map(sample_list, ~.x[var_feats, ]),
    # k.anchor = 3,
    # k.filter = 5,
    # k.score = 3,
    # k.anchor = 5,
    # k.filter = 50,
    # k.score = 10,
    # reference = reference,
    reference = which(names(sample_list) == reference),
    dims = 1:max_dims,
    reduction = 'rpca',
    # reduction = 'reference',
    # k.anchor = 2,
    k.anchor = 1,
    k.filter = 1,
    k.score = 1,
    nn.method = 'annoy',
    normalization.method = 'SCT',
    anchor.features = var_feats
  )

  so_I <- Seurat::IntegrateData(
    anchorset = anchors,
    dims = 1:max_dims,
    k.weight = 2
    # normalization.method = 'SCT'
  )

  DefaultAssay(so_I) <- 'integrated'
  so_I[['RNA']] <- NULL
  so_I[['SCT']] <- NULL

  VariableFeatures(so_I) <- var_feats

  so_I <- so_I[, apply(so_I@assays$integrated[,], 2,
    function(x) any(!is.na(x)))]
  so_I <- scale_so(so_I, feature_weights = feature_weights)
  so_I <- scale_so(so_I, feature_weights = 'scale')
  # so_I <- ScaleData(so_I, do.center = T, do.scale = TRUE)
  so_I <- RunPCA(so_I, npcs = N_pcs,
    dims = 1:max_dims, verbose = FALSE)

  so_I
})


seurat_cca_code <- rlang::expr({
  stopifnot(!is.null(reference) && reference %in% experiments)
  var_feats <- VariableFeatures(so_P)
  sample_list <- SplitObject(so_P, split.by = 'exp')
  rm(so_P); gc()

  stop('Incomplete implementation')
  max_dims <- min(min(purrr::map_int(sample_list, ncol))-1, 30)
})



precise_code <- rlang::expr({
  if (length(unique(so_P$exp)) > 2) {
    return(NULL)
  }
  sample_list <- SplitObject(so_P, 'exp')
  # sample_list <- map()
  bulk_idx <- which(names(sample_list) %in% bulk_experiments)
  sc_idx <- which(names(sample_list) %in%
    c(tar_read(e), c('5310_in_vivo', '5310_in_vitro')))

  if (prefilter_ref) {
    exp_rle <- rle(so_P@meta.data$exp) %>%
      { set_names(.$lengths, .$values) }
    corM <- cor(GetAssayData(so_P)) %>%
      { .[1:exp_rle[1], (exp_rle[2]+1):sum(exp_rle)] }
    thresh <- quantile(as.vector(corM), .75)
    # sample_idxs <- which(apply(corM, bulk_idx, max) >= thresh)
    exclude_sample_idxs <- names(which(apply(corM, 2, max) < thresh))
    N_pre <- ncol(so_P)
    so_P <- so_P[, which(!rownames(so_P@meta.data) %in% exclude_sample_idxs)]
    sample_list <- SplitObject(so_P, 'exp')
    so_P <<- so_P
    N_post <- ncol(so_P)
    if (N_post == N_pre) browser()
    # sample_list[[bulk_idx]] <- sample_list[[bulk_idx]][, sample_idxs]
  }

  max_dims <- min(min(purrr::map_int(sample_list, ncol))-1, 30)
  if (!is.null(N_pcs)) {
    max_dims <- min(max_dims, N_pcs)
  }

  if (F) {
    plot_var_exp_panel(
      M_source = so2M(sample_list[[bulk_idx]]),
      M_target = so2M(sample_list[[sc_idx]]),
      # N_source = min(N_PV_bulk, max_dims),
      # N_target = min(N_PV_sc, max_dims)
      N_source = 10L,
      N_target = 10L
    )

    p <- plot_cosine_sim(
      M_source = so2M(sample_list[[bulk_idx]]),
      M_target = so2M(sample_list[[sc_idx]]),
      # N_source = min(N_PV_bulk, max_dims),
      # N_target = min(N_PV_sc, max_dims)
      N_source = 10L,
      N_target = 10L
    )
    print_plot(p, file.path(img_dir, glue('PV_cosine_sim.png')))
    
    # var_dtf <- cbind(p_int$bulk_var$var, p_int$sc_var$var)
    
    source('~/MirjamHoekstra/R/init.R')
    plot_cosine_heatmap(
      M_source = so2M(sample_list[[bulk_idx]]),
      M_target = so2M(sample_list[[sc_idx]]),
      # N_source = min(N_PV_bulk, max_dims),
      # N_target = min(N_PV_sc, max_dims)
      N_source = 10L,
      N_target = 10L
    )
  }

  p_int <- PRECISE_integrate$new(
    bulk_so = sample_list[[bulk_idx]],
    sc_so = sample_list[[sc_idx]],
    N_PV_bulk = min(N_PV_bulk, max_dims),
    N_PV_sc = min(N_PV_sc, max_dims)
  )
  browser()

  # apply(var_dtf, 1, mean)
  precise_Ms <-
    list(p_int$bulk_P_a, p_int$sc_P_a) %>%
    purrr::map(~as.matrix(dplyr::select(.x, matches('CF\\d+')))) %>%
    { .[c(bulk_idx, sc_idx)] } %>%
    # purrr::map(~log2(.x + 1)) %>%
    { . }

  # map(precise_Ms, ~rowSums(abs(.x)))
  if (perform_mnn) {
    gen_sce <- function(M) {
      SingleCellExperiment(list('logcounts' = t(M)))
    }
    # gen_sce(precise_Ms[[1]])
    sces <- purrr::map(precise_Ms, gen_sce)
    # MNN_sce <- rlang::exec(mnnCorrect, !!sces,
    #   k = min(map_dbl(precise_Ms, ~ceiling(nrow(.x) / 10))))
    MNN_sce <- mnnCorrect(sces[[1]], sces[[2]],
      k = k_mult * min(map_dbl(precise_Ms, ~ceiling(nrow(.x) / 10))))
    precise_M <- t(assay(MNN_sce, 'corrected'))
  } else {
    precise_M <- precise_Ms %>%
      { rbind(.[[1]], .[[2]]) } %>%
      { set_colnames(., paste0('precise_', 1:ncol(.))) } %>%
      { . }
  }

  # rowSums(abs(precise_M))
  so_P[['precise']] <-
    CreateDimReducObject(
      embeddings = precise_M,
      key = 'precise_',
      assay = DefaultAssay(so_P)
    )

  so_P
})


embed_experiments <- function(
  experiments = c('5029', '5310'),
  reference = '5029',
  split_5310 = T,
  sc_mode = 'SCT',
  # norm_method = 'CPM',
  norm_method = 'none',
  embedding = 'raw',
  # feature_weights = 'prop',
  # feature_weights = 'scale',
  feature_weights = 'none',
  annotate_umap = TRUE,
  GDR_thresh = NULL,
  filter_gene_reproducibility = NULL,
  genelist = NULL,
  genes = NULL,
  N_PV_bulk = 3L,
  N_PV_sc = 3L,
  prefilter_ref = T,
  k_mult = 2,
  N_pcs = NULL,
  data_mod_code = NULL,
  perform_mnn = FALSE,
  reduce_combine = FALSE) {

  library(Seurat)
  if (is.null(experiments)) return(NULL)

  so_P <-
    read_preproc_experiments(
      experiments = experiments,
      sc_mode = sc_mode
    ) %>%
    normalize_experiments(
      split_5310 = split_5310,
      merge_experiments = T,
      norm_method = norm_method,
      reduce_combine = reduce_combine,
      data_mod_code = data_mod_code,
      GDR_thresh = NULL,
      genelist = genelist,
      genes = genes
    )
  gc()

  if (F) {
    ## After sample scaling/normalization has been done, we can prune
    ## features we don't need/want
    so_P <- subset_feats(
      so = so_P,
      GDR_thresh = GDR_thresh,
      filter_gene_reproducibility = filter_gene_reproducibility,
      genelist = genelist
    )
    gc()
    VariableFeatures(so_P) <- detected_genes(so_P)
  } else if (T) {
    so_P <- set_var_genes(
      so = so_P,
      genelist = genelist,
      genes = genes,
      GDR_thresh = GDR_thresh,
      filter_gene_reproducibility = filter_gene_reproducibility
    )
    # VariableFeatures(so_P)
    # dim(so_P)
    gc()
  }
  if (!nrow(so_P) || !length(VariableFeatures(so_P))) return(NULL)

  so_P <- scale_so(so_P, feature_weights = feature_weights)

  if (T && is.null(N_pcs)) {
    N_pcs <- min(dim(so_P) - 1) %>%
      min(length(VariableFeatures(so_P))) %>%
      min(50L)
    # DefaultAssay(so_P)
    so_P <- suppressWarnings(RunPCA(so_P, npcs = N_pcs, verbose = F))
    # so_P@reductions$pca[,]
    gc()
  }

  ## Avoid having to rerun already run embeddings when a new embedding
  ## is added. Downside is that this method of code lookup precludes
  ## the targets pipelining system to detect and account for changes
  ## in these code blocks
  # so_I <- rlang::tidy_eval(get(glue::glue('{embedding}_code')))
  so_I <- tryCatch(
    eval(get(glue::glue('{embedding}_code'))),
    error = function(e) { print(e); print(traceback()); NULL }
  )
  if (is.null(so_I)) return(NULL)
  rm(so_P); gc()

  reduction <- c(
    'harmony' = 'harmony',
    'raw' = 'pca',
    'seurat_rpca' = 'pca',
    'precise' = 'precise')[embedding]
  if (annotate_umap) {
    so_U <- tryCatch(
      RunUMAP(so_I, reduction = reduction, dims = 1:N_pcs),
      error = function(e) { print(e); so_I })
  } else {
    so_U <- so_I
  }
  rm(so_I)
  map(so_U@reductions, dim)

  dtf <- tryCatch({
    so2dtf(so_U) %>%
      # add_umap(matches('^PC_\\d+$')) %>%
      numerify_regressors() %>%
      # order_duration() %>%
      separate_duration() %>%
      order_condition_name() %>%
      dplyr::mutate(experiment = exp) %>%
      dplyr::mutate(embedding = embedding) %>%
      dplyr::mutate(feature_weights = feature_weights) %>%
      dplyr::mutate(GDR_thresh = GDR_thresh) %>%
      dplyr::mutate(filter_gene_reproducibility =
        filter_gene_reproducibility) %>%
      dplyr::mutate(genelist = genelist)
  }, error = function(e) { print(e); NULL })

  rm(so_I); gc()

  return(dtf)
}


## These do not play well with the targets system it seems
embed_raw_scaled <-
  pryr::partial(embed_experiments,
    embedding = 'raw', feature_weights = 'scale')
# embed_raw_scaled()


embed_harmony_scaled <-
  pryr::partial(embed_experiments,
    embedding = 'harmony', feature_weights = 'scale')


embed_raw_prop <-
  pryr::partial(embed_experiments,
    embedding = 'raw', feature_weights = 'prop')
# embed_raw_prop()


embed_harmony_prop <-
  pryr::partial(embed_experiments,
    embedding = 'harmony', feature_weights = 'prop')


embed_combat <- function(
  experiments = c('5029', '5310'),
  GDR_thresh = NULL,
  filter_gene_reproducibility = 4L,
  genelist = 'informativeV15_monotonic') {

  library(Seurat)
  if (is.null(experiments)) return(NULL)
  so_P <- read_preproc_experiments(
    experiments = experiments
  )

  library(sva)
  batch <- so_P@meta.data$exp
  stopifnot(all(!is.na(batch)))
  M_orig <- GetAssayData(so_P, 'counts', assay = 'SCT')
  M_corr <- sva::ComBat_seq(M_orig, batch = batch)
  # M_corr <- M_orig
  stopifnot(dim(M_orig) == dim(M_corr))
  rm(M_orig)
  # summary(map_dbl(1:nrow(M_orig), ~cor(M_orig[, .x], M_corr[, .x])))

  so_P[['combat']] <- CreateAssayObject(counts = M_corr)

  ## Subset to desired features after having done batch effect
  ## correction
  feats <- read_preproc_experiments(
      experiments = experiments, merge = F,
    ) %>%
    subset_feats(
      genelist = genelist,
      GDR_thresh = GDR_thresh,
      filter_gene_reproducibility = filter_gene_reproducibility
    ) %>%
    rownames()
  so_P <- so_P[feats, ]

  DefaultAssay(so_P) <- 'combat'
  VariableFeatures(so_P) <- rownames(so2M(so_P))
  so_P <- ScaleData(so_P)
  so_P <- RunPCA(so_P, assay = 'combat')

  dtf <- so2dtf(so_P,
    coord_extract = function(so) {
      so@reductions$pca@cell.embeddings
    }) %>%
    add_umap(matches('^PC_\\d+$')) %>%
    numerify_regressors() %>%
    order_duration() %>%
    separate_duration() %>%
    order_condition_name() %>%
    dplyr::mutate(experiment = exp) %>%
    dplyr::mutate(embedding = 'raw')

  rm(so_P)
  return(dtf)
}


embed_quantile <- function(
  experiments = c('5029', '5310'),
  GDR_thresh = NULL,
  sc_mode = 'SCT',
  split_5310 = T,
  genes = NULL,
  score_embedding = F,
  filter_gene_reproducibility = 4L,
  genelist = 'informativeV15') {

  library(Seurat)
  if (is.null(experiments)) return(NULL)
  ## Don't do any gene filtering before quantile norm
  so_P <- read_preproc_experiments(
    GDR_thresh = NULL,
    experiments = experiments,
    sc_mode = sc_mode,
    split_5310 = split_5310,
    merge = T, genelist = NULL)

  M_orig <- as.matrix(GetAssayData(so_P, 'counts', assay = 'SCT'))
  M_corr <- preprocessCore::normalize.quantiles(M_orig)
  dimnames(M_corr) <- dimnames(M_orig)

  stopifnot(dim(M_orig) == dim(M_corr))
  rm(M_orig)

  so_P[['quantile']] <- CreateAssayObject(counts = M_corr)
  DefaultAssay(so_P) <- 'quantile'

  if (!is.null(genelist)) {
    feats <- subset_feats(
      so_P[['quantile']], genelist = genelist) %>% rownames()
    so_P <- so_P[feats, ]
    VariableFeatures(so_P) <- rownames(so2M(so_P))
  } else {
    so_P <- FindVariableFeatures(so_P)
  }

  if (!is.null(genes)) {
    so_P <- subset_feats(so_P, genes)
  }

  so_P <- ScaleData(so_P)
  so_P <- RunPCA(so_P, assay = 'quantile')

  dtf <- so2dtf(so_P,
    coord_extract = function(so) {
      so@reductions$pca@cell.embeddings
    }) %>%
    add_umap(matches('^PC_\\d+$')) %>%
    numerify_regressors() %>%
    order_duration() %>%
    separate_duration() %>%
    order_condition_name() %>%
    dplyr::mutate(experiment = exp) %>%
    dplyr::mutate(embedding = 'raw')

  if (score_embedding) {
    scores <- score_PCA(so_P)
    attr(dtf, 'scores') <- scores
  }

  rm(so_P)
  return(dtf)
}


# embed_raw(experiments = c('5029', '5310'))
embed_seurat_rpca <- function(
  experiments = c('5029', '5310'),
  genelist = 'informativeV15_monotonic',
  GDR_thresh = NULL,
  reference = '5310',
  filter_gene_reproducibility = 4L,
  N_pcs_informative = 15) {

  library(Seurat)
  if (is.null(experiments)) return(NULL)

  sample_list <- read_preproc_experiments(
    experiments,
    genelist = genelist,
    GDR_thresh = GDR_thresh,
    merge_experiments = FALSE
  )

  shared_features <-
    purrr::map(sample_list, ~rownames(.x@assays$SCT)) %>%
    purrr::reduce(intersect) %>%
    intersect(read_geneset(genelist))

  stopifnot(!is.null(reference) && reference %in% names(experiments))
  reference_dataset <- which(names(sample_list) == reference)

  anchors <- Seurat::FindIntegrationAnchors(
    object.list = sample_list,
    # k.anchor = 3,
    # k.filter = 5,
    # k.score = 3,
    # k.anchor = 5,
    # k.filter = 50,
    # k.score = 10,
    reference = reference_dataset,
    # dims = 1:20,
    reduction = 'rpca',
    # reduction = 'reference',
    nn.method = 'annoy',
    normalization.method = 'SCT',
    anchor.features = shared_features
  )

  so_P <- Seurat::IntegrateData(
    anchorset = anchors,
    # k.weight = k_weight,
    normalization.method = 'SCT'
  )

  VariableFeatures(so_P) <- rownames(so2M(so_P))
  so_P <- RunPCA(so_P, verbose = F)

  # so_P <- RunUMAP(so_P, reduction = 'pca',
  #   dims = 1:N_pcs_informative)
  # so_P <- FindNeighbors(so_P,
  #   dims = 1:N_pcs_informative, reduction = 'pca')
  # so_P <- FindClusters(so_P, verbose = FALSE, resolution = 2)
  # so_P <- identify_outlying_clusters(so_P)

  out <- so2dtf(so_P,
    coord_extract = function(so) {
      so@reductions$pca@cell.embeddings
    }) %>%
    add_umap(matches('^PC_\\d+$')) %>%
    numerify_regressors() %>%
    separate_duration() %>%
    order_duration() %>%
    order_stim_group() %>%
    order_condition_name() %>%
    dplyr::mutate(experiment = exp) %>%
    dplyr::mutate(embedding =
      glue::glue('seurat_rpca\\
        {make_flag(N_pcs_informative)}\\
        {make_flag(reference)}')
    )

  return(out)
}


subsample_sc <- function(bm_dtf, N_sc = 1000L) {
  if (maartenutils::null_dat(bm_dtf))
    return(NULL)

  sc_idxs <-
    which(
      bm_dtf$sample_type == 'sc' &
      bm_dtf$exp %in% c('5310', '5310_in_vitro', '6369', '6489')
    ) %>%
    { base::sample(., min(N_sc, length(.)), replace = F) }

  bulk_idxs <- which(bm_dtf$sample_type == 'bulk') %>%
    union(which(bm_dtf$exp %in% bulk_experiments))

  bm_dtf$subsampled <- 1:nrow(bm_dtf) %in% c(bulk_idxs, sc_idxs)

  return(bm_dtf)
}

if (F) {
  bm_dtf <- tar_read(bm_dtf_raw_5029.6369.6434_informativeV15)
  subsample_sc(bm_dtf)
}


if (F) {
  #' DEPRECATED
  #'
  #'
  subset_bm_targets <- function(target_list = scVI_targets) {
    target_names <-
      purrr::map_chr(target_list, extract_target_name) %>%
      stringr::str_detect('bm_dtf')
    target_list[target_names]
  }


  get_embedding_feat_regex <- function(embedding) {
    c('harmony' = 'harmony',
      'raw' = 'PC',
      'precise' = 'precise',
      'seurat_rpca' = 'PC',
      'raw_scaled' = 'PC',
      'raw_prop' = 'PC',
      'harmony_scaled' = 'harmony',
      'harmony_prop' = 'harmony')[as.character(embedding)]
  }


  #' Create heatmap of genes vs. prediction scores
  #'
  #' @param gs Geneset string, geneset will be looked up from targets
  #' database.This argument will be overruled by the \param{genes} if
  #' that one is non-NULL
  #' @param genes Long-formatted tibble/data.frame of two columns:
  #' geneset and column
  #'
  scores_HM <- function(
    experiment = '5310',
    response_var = 'ifn_duration',
    gs = 'informativeV15',
    genes = NULL,
    bm_dtf,
    mod_idx = 1L,
    GDR_thresh = .7,
    hill = NULL) {

    so <- tar_read(filtered_cleaned_so, branch = e2i(experiment))[[1]]

    if (is.null(genes) || is.na(genes) || length(genes) == 0) {
      genes <-
        tibble::enframe(tar_read_raw(gs), 'geneset', 'gene') %>%
        tidyr::unnest(cols = c(gene)) %>%
        dplyr::distinct(gene, .keep_all = T)
    }

    # detected_genes <-
    #   get_GDR_table(so, genes = genes$gene,
    #     thresh = GDR_thresh, format = 'table')
    l_detected_genes <-
      get_GDR_table(so, genes = genes$gene, thresh = GDR_thresh)

    if (length(l_detected_genes) == 0) {
      rlang::warn('No genes')
      return(NULL)
    }

    gene_stats <- compute_gene_stats(so, l_detected_genes)
    l_detected_genes <-
      gene_stats %>%
      dplyr::filter(mean >= 3) %>%
      dplyr::pull(gene)
    # gene_stats_M <- gene_stats %>%
    #   dplyr::left_join(tibble(gene = rownames(scoreM)),
    #     by = 'gene') %>%
    #   dplyr::select(mean, max) %>%
    #   { . }

    if (length(l_detected_genes) == 0) {
      rlang::warn('No genes')
      return(NULL)
    }

    genes <-
      dplyr::right_join(
        genes,
        tibble(gene = l_detected_genes), by = 'gene')

    scoreM <- subset_feats(so2M(so), l_detected_genes)
    if (T) {
      scoreM <- t(apply(scoreM, 1, min_max_scaling))
    } else {
      scoreM <- t(scale(t(scoreM)))
    }
    if (!is.null(hill) && !is.na(hill)) {
      scoreM <- sign(scoreM) * abs(scoreM)^hill / (1 +
        abs(scoreM)^hill)
    }

    if (experiment == '5310') {
      l_experiment <- '5310_in_vitro'
    } else {
      l_experiment <- experiment
    }
    l_dtf <- bm_dtf %>%
      dplyr::filter(.data[['experiment']] == .env[['l_experiment']])
    # dplyr::filter(stringr::str_detect(experiment, {{experiment}}))

    model <- train_knn_reg_mod(
      dtf = bm_dtf,
      filtering_settings = get_filtering_settings(experiment),
      response_var = response_var,
      mod_idx = mod_idx
    )
    l_dtf$pred <- unlist(predict(model, l_dtf)[, 1])

    ## Identially order M and annotation; arrange by predicted score
    shared_ids <- intersect(colnames(scoreM), l_dtf$sample_id)
    l_dtf <- l_dtf %>%
      dplyr::filter(sample_id %in% shared_ids)
    scoreM <- scoreM[, shared_ids]
    l_dtf <- dplyr::arrange(l_dtf, duration, stim_group, pred)
    # l_dtf <- l_dtf[match(shared_ids, l_dtf$sample_id), ]
    scoreM <- scoreM[, match(l_dtf$sample_id, colnames(scoreM))]
    stopifnot(l_dtf$sample_id == colnames(scoreM))

    ra <-
      dplyr::select(l_dtf, stim_group, ifn_duration,
        tnf_duration, pred) %>%
    dplyr::rename('pred_{response_var}' := pred)
  show_ra_legend <- rlang::set_names(c(TRUE, FALSE, FALSE, TRUE),
    colnames(ra))
  ca <-
    dplyr::left_join(genes,
      tibble(gene = rownames(scoreM)),
      by = 'gene') %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(geneset)
HM <- gen_HM(
  M = t(scoreM),
  cluster_columns = T,
  cluster_rows = F,
  use_raster = F,
  raster_by_magick = F,
  show_column_names = T,
  ra = ra,
  show_ra_legend = show_ra_legend,
  row_ann_f = gen_sample_annotation_HM,
  ca = ca,
  col_ann_f = gen_sample_annotation_HM
)
o_fn <- file.path(Sys.getenv('img_dir'),
  glue::glue('score_vs_genes\\
    {make_flag(experiment)}\\
    {make_flag(response_var)}\\
    {make_flag(gs)}.png'))
    names(HM)
    # HM_C <- RA2 %v% HM
    print_plot_eval(
      {
        draw(HM, show_heatmap_legend = T,
          heatmap_legend_side = 'bottom',
          merge_legends = T)
      },
      filename = o_fn)

    return(o_fn)
  }


  score_PCA <- function(so) {
    id_vars = c('tnf_conc', 'ifn_conc', 'duration')

    meta <- fill_in_tnf(so@meta.data) %>%
      dplyr::mutate(sample_idx = 1:n())

    dup_settings <- meta[which(duplicated(meta[, id_vars])), id_vars] %>%
      set_rownames(NULL) %>%
      na.omit() %>%
      dplyr::distinct() %>%
      # head(n = 5) %>%
      { . }

    # str(so@reductions$pca)
    # so@reductions$pca
    # so@reductions$pca@stdev %>% { .^2 / sum(.^2) }

    M <- so@reductions$pca@cell.embeddings
    ## Scores are already weighted; PCs with high variances will have
    ## bigger scores
    # total_var <- so@reductions$pca@stdev %>% { .^2 } %>% sum()
    # total_score <- sum(colSums(abs(M)))
    distM <- dist(M)
    max_dist <- quantile(distM, .99)
    min_dist <- quantile(distM, .01)
    denum <- max_dist - min_dist

    dup_settings_ann <-
      purrr::map_dfr(1:nrow(dup_settings), function(i) {
        l_meta <- dplyr::right_join(meta, dup_settings[i, ],
          by = colnames(dup_settings))
        id_cols <- l_meta[, c('exp', 'condition_name')]
        stopifnot(all(!duplicated(id_cols)))
        return(list(
            'N_samples' = nrow(l_meta),
            'score' = (mean(dist(M[l_meta$sample_idx, ])) - min_dist) / denum
            ))
      })

    out <- bind_cols(dup_settings, dup_settings_ann)

    return(out)
  }


  #' Create the variable part of ML targets from the first row of
  #' data.frame/tibble
  #'
  #'
  settings_to_varpart <- function(settings) {
    var_part <- with(settings, {
      glue::glue('
        _{embedding[1]}\\
        _{params$included_experiments}\\
        _{genelist[1]}\\
        _{feature_weights[1]}\\
        _{norm_method[1]}\\
        _{GDR_thresh[1]}\\
        _{filter_gene_reproducibility[1]}'
      )
      })
  }


  annotate_knn_predictions <- function(
    dtf,
    response_vars = tar_read(response_vars),
    filtering_settings = tar_read(benchmark_sc_filtering)[['no_sn']],
    mod_idx = 1L) {

    for (rv in response_vars) {
      mod <- train_knn_reg_mod(
        dtf = dtf,
        filtering_settings = filtering_settings,
        response_var = rv,
        mod_idx = mod_idx
      )
      dtf[[glue::glue('{rv}_pred')]] <- unlist(predict(mod, dtf))
    }

    actual <- norm_regressors(dtf)[, response_vars]
    error <- actual - dtf[, glue::glue('{response_vars}_pred')]
    colnames(error) <- paste0(colnames(error), '_error')

    dtf <- cbind(dtf, error)

    return(dtf)
  }

  summarize_pred_error <- function(
    dtf,
    experiments = c(tar_read(e), '5310_in_vitro'),
    group_vars = c('experiment', 'stim_group', 'duration'),
    response_vars = tar_read(response_vars),
    output_mode = 'SE') {

    source('~/MirjamHoekstra/R/init.R')
    error_vars <- glue::glue('{response_vars}_error')
    abs_mean <- function(x) mean(abs(x))
    f_list <- list('mean' = mean, 'CI_l' = CI_l)
    f_list <- list('mean' = mean, 'CI_l' = CI_l, 'CI_h' = CI_h)
    f_list <- list('mean' = mean, 'CI_l' = CI_l, 'CI_h' = CI_h, 'CV' = CV)
    # f_list <- list('mean' = mean)
    out <- dtf %>%
      dplyr::filter(experiment %in% experiments) %>%
      dplyr::mutate(across(error_vars, abs)) %>%
      dplyr::group_by_at(group_vars) %>%
      # dplyr::summarize(across(), N = n()) %>%
      dplyr::summarize(across(error_vars, f_list)) %>%
      { . }

    if (output_mode == 'tibble') {
      return(out)
    } else if (output_mode == 'SE') {
      ## Separate the different types of output in the columns of out
      ## into separate assays
      M <- as.matrix(out[, (length(group_vars)+1):ncol(out)])
      types <-
        stringr::str_replace(colnames(M), '.*_error_(.*)$', '\\1')
      idxs <- map(auto_name(unique(types)), ~which(types == .x))

      ## Strip away column information that is now not needed anymore
      strip_cn <- function(x) set_colnames(x,
        stringr::str_replace(colnames(x), '_error.*$', ''))
      l_assays <- imap(idxs, ~strip_cn(M[, .x]))

      n_out <- SummarizedExperiment::SummarizedExperiment(
        assays = l_assays,
        rowData = as.data.frame(out[, 1:length(group_vars)])
      )
      # assay(n_out, 'mean')
      # colData(n_out)
      return(n_out)
    }
  }

  format_tune_res <- function(tune_res) {
    l_id_vars <- c('response_var',
      'experiment', 'experiment_string',
      'embedding', 'genelist', 'GDR_thresh', 'mod_idx',
      'filter_gene_reproducibility', 'feature_weights', 'norm_method') %>%
    intersect(colnames(tune_res)) %>%
    c('mod_idx')

  if ('split_type' %in% colnames(tune_res)) {
    tune_res <- tune_res %>%
      dplyr::filter(split_type == 'assessment') %>%
      dplyr::select(-split_type) %>%
      { . }
  }

  tune_res <-
    tune_res %>%
    dplyr::filter(!is.na(MAE)) %>%
    dplyr::filter(!is.na(GDR_thresh)) %>%
    dplyr::left_join(
      tar_read(model_hyperparam_grid),
      by = 'mod_idx'
      ) %>%
    dplyr::mutate(across(
        any_of(c('embedding', 'feature_weights', 'norm_method',
            'experiment', 'mod_idx', 'split_idx')),
        order_factor)
      ) %>%
    dplyr::mutate(across(
        any_of(c('GDR_thresh', 'filter_gene_reproducibility')),
        numeric2factor)
      ) %>%
    dplyr::mutate(across(
        any_of(c('genelist', 'norm_method', 'genelist',
            'experiment_string', 'weight_func')),
        length_order_factor)
      ) %>%
    dplyr::mutate(across(
        any_of(c('neighbors', 'k', 'dist_power', 'N_spec_levels')),
        numeric2factor)
      ) %>%
    dplyr::mutate(across(
        any_of(c('weight_func', 'response_var')),
        length_order_factor)
      ) %>%
    { . }

  attr(tune_res, 'id_vars') <- l_id_vars
  return(tune_res)
  }
}
