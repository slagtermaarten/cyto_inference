perform_seurat_integration <- result_cacher(
  f = function(
    ds = NULL,
    ## Bulk experiments
    include_5029 = F,
    include_6434 = F,
    ## Single cell experiments
    include_5310 = F,
    include_6369 = F,
    include_6489 = F,
    include_6493 = F,
    hashtag = 'HTO_classification',
    genelist = '',
    regress_only_6369 = F,
    reference = 'exp6434_bulk',
    downsample_6369 = NULL,
    split_5310 = T,
    integration_method = 'reference',
    k_weight = 10,
    add_bulkified_sc = T,
    optimize_silhouette = F,
    filtering_opts = default_filtering_opts,
    source_type = 'filtered_ig',
    vars_to_regress = NULL) {

    sample_list <- call_with(compile_seurat_object_list,
      map(ls(), ~get(.x, envir = parent.frame())))

    elbows <-
      map(sample_list, ElbowPlot,ndims = 20, reduction = 'pca')
    stopifnot(!is.null(reference) && reference %in% names(sample_list))
    reference_dataset <- which(names(sample_list) == reference)

    if (integration_method == 'reference') {
      anchors <- FindIntegrationAnchors(
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
        anchor.features = i_features)

      so <- IntegrateData(
        anchorset = anchors,
        k.weight = k_weight,
        normalization.method = 'SCT'
      )

      if (F) {
        head(anchors@anchors)
        anchors@anchors %>% group_by(dataset1, dataset2) %>%
          dplyr::summarize(
            across(score, c('mean' = mean, 'median' = median)),
            .groups = 'drop')

        anchors@anchors %>%
          group_by(dataset1) %>%
          dplyr::summarize(
            across(score,
              c('mean' = mean, 'median' = median)),
            .groups = 'drop')

        # Seurat:::VlnPlot
        # Seurat:::MetageneBicorPlot
        # Seurat:::CalcAlignmentScore
      }
    } else if (integration_method == 'simple_stack') {
      ## The merge 'y' argument is not working for me, so I do it
      ## iteratively
      so <- purrr::reduce(sample_list, merge, merge.data = T)
      so@assays$integrated <- so@assays$SCT
      VariableFeatures(so) <-
        reduce(map(sample_list, VariableFeatures), intersect)
    } else {
      browser()
    }

    so <- RunPCA(so, verbose = FALSE)

    if (F) {
      p <- ElbowPlot(
        so, ndims = 20, reduction = 'pca')
      print_plot(p, fn = file.path(img_dir, 'integration_elbow.png'))
    }

    N_pcs_informative = 15
    so <- RunUMAP(so, reduction = 'pca', dims = 1:N_pcs_informative)
    so <- FindNeighbors(so,
      dims = 1:N_pcs_informative, reduction = 'pca')
    so <- FindClusters(so, verbose = FALSE, resolution = 2)
    so <- identify_outlying_clusters(so)

    if (optimize_silhouette) {
      ## TODO Finish
      ## Find optimal clustering resolution by optimizing silhouette
      ## score
      silhouette_scores <-
        seq(.8, 2, by = .2) %>%
        auto_name %>%
        purrr::map_dbl(function(res) {
          res <- round(res, 1)
          so <- FindClusters(so, reduction = 'pca', resolution = res)
          clusters <- so@meta.data[[glue('integrated_snn_res.res.{res}')]]
          silhouette()
        })
    }
    return(so)
  },
  filename = function() {
    o_dir <- file.path(rds_dir, 'seurat_integrations')
    dir.create(o_dir, showWarnings = F)
    fn <- file.path(o_dir, glue('
      {make_flag(ds)}\\
      {make_flag(genelist)}\\
      {make_flag(include_6369)}\\
      {make_flag(include_6489)}\\
      {make_flag(include_6493)}\\
      {make_flag(regress_only_6369)}\\
      {make_flag(downsample_6369)}\\
      {make_flag(integration_method)}\\
      {ifelse(integration_method==\'reference\',
        make_flag(reference_dataset), \'\')}\\
      {make_flag(source_type)}.rds'))
    fn
    # saveRDS(so, fn)
  },
  min_mod_time = '2021-09-13 18:05'
)


compile_seurat_object_list <- function() {
  source_fns <- map_chr(auto_name(experiments), function(experiment) {
    compile_so_fns(
      experiment = experiment,
      filtering_opts = filtering_opts,
      lookup_mode = T)[[source_type]]
    }) %>%
    { . }

  sample_list <- list()

  if (include_5029) {
    exp5029 <- readRDS(file.path(rds_dir, 'exp5029_SCT.rds'))
    sample_list <- c(sample_list, list('exp5029_bulk' = exp5029))
  }

  if (include_6434) {
    exp6434 <- readRDS(file.path(rds_dir, 'exp6434_SCT.rds'))
    sample_list <- c(sample_list, list('exp6434_bulk' = exp6434))
  }

  if (include_5310) {
    exp5310_sc <- readRDS(source_fns['5310'])
    exp5310_sc <-
      exp5310_sc[, !is.na(exp5310_sc@meta.data$sample_origin)]
    if (split_5310) {
      es <-
        SplitObject(exp5310_sc, split.by = 'sample_origin') %>%
        set_names(paste0('exp5310_sc_', names(.)))
      sample_list <- c(sample_list, es)
    } else {
      sample_list <- list('exp5310_sc' = exp5310_sc)
    }
  }

  ## Downsample in vivo cells from 5310 to investigate required
  ## depth for follow up experiments
  if (!is.null(ds)) {
    lib_sizes <-
      colSums(as.matrix(sample_list[['in_vivo']]@assays$RNA@counts))
    sampled_counts <- SampleUMI(
      sample_list[['in_vivo']]@assays$RNA@counts,
      max.umi = ds * lib_sizes)
    sample_list[['in_vivo']]@assays$RNA@counts <- sampled_counts
  }

  if (include_6369) {
    sample_list <- c(sample_list,
      list('exp6369_sc' = readRDS(source_fns['6369'])))
    if (!is.null(downsample_6369)) {
      tmp <- sample_list[['exp6369_sc']]
      idxs <- sample(seq(1, ncol(tmp)),
        size = downsample_6369, replace = F)
      tmp <- tmp[, idxs]
      sample_list[['exp6369_sc']] <- tmp
    }
    # sample_list[['exp6369_sc']]@meta.data$G2M.Score <- 1
    # sample_list[['exp6369_sc']]@meta.data$S.Score <- 1
  }

  if (include_6489) {
    sample_list <- c(sample_list,
      list('exp6489_sc' = readRDS(source_fns['6489'])))
  }

  if (include_6493) {
    sample_list <- c(sample_list,
      list('exp6493_sc' = readRDS(source_fns['6493'])))
  }

  if (add_bulkified_sc) {
    ## TODO Finish
    sc_exps <- grep('_sc', names(sample_list), value = T)
    pb_exps <- map(sample_list[sc_exps],
      ~pseudo_bulk(.x, u_var = 'condition_name',
        assay = 'SCT', datatype = 'counts')) %>%
      set_names(gsub('_sc', '_pb', sc_exps))
    sample_list <- c(sample_list, pb_exps)
    # map(sample_list[sc_exps], function(ex) {
    #   sample_levs <-
    #     unique(ex@meta.data[, c('stim_group', 'duration')])
    #   exp_means <- sample_levs %>%
    #     purrr::pmap(function(stim_group, duration) {
    #       cell_idxs <-
    #         find_cells(so = ex, grp = stim_group, duration = duration)
    #       M <- GetAssay(ex)[, cell_idxs]
    #       apply(M, 1, mean)
    #     })
    # ex@meta.data %>%
    #   right_join(sample_levs[1, ]) %>%
    #   dplyr::summarize(ac
    # head(ex@meta.data, n = 1)
    # browser()
    # })
    # Seurat::AggregateExpression
  }

  ## This shouldn't be necessary but alas
  exps_requiring_SCT <-
    which(map_lgl(sample_list, ~'SCT' %nin% names(.x@assays))) %>%
    names()
  for (ex in exps_requiring_SCT) {
    sample_list[[ex]] <- SCTransform(sample_list[[ex]])
  }

  ## Reduce the data to all those detectable in all datasets
  shared_features <-
    purrr::map(sample_list, ~rownames(.x@assays$SCT)) %>%
    purrr::reduce(intersect)
  if (!is.null(genelist) && genelist != '') {
    genelist_genes <- read_geneset(glue('{genelist}_genes'))
    shared_features <- intersect(shared_features, genelist_genes)
  }
  sample_list <- purrr::map(sample_list, ~.x[shared_features, ])

  if (prep_integration) {
    i_features <- SelectIntegrationFeatures(
      object.list = sample_list,
      nfeatures = min(3000, nrow(sample_list[[1]])))
    sample_list <- PrepSCTIntegration(
      object.list = sample_list,
      anchor.features = i_features)
    sample_list <- purrr::map(sample_list, RunPCA,
      features = i_features)
  }
  return(sample_list)
}
formals(compile_seurat_object_list) <-
  formals(perform_seurat_integration)
formals(compile_seurat_object_list)$prep_integration <- F



project_data <- function(
  ref_so = sample_list[['exp6434_bulk']],
  query_so = sample_list[['exp5310_sc_in_vivo']],
  N_dims = 15L, k.filter = 8, k.weight = 4) {

  ref_so@meta.data <- add_condition_name(ref_so@meta.data)

  anchors <- Seurat::FindTransferAnchors(
    reference = ref_so,
    query = query_so,
    dims = 1:N_dims,
    k.filter = k.filter,
    reference.reduction = 'pca')
  N_anchors <- nrow(anchors@anchors)
  message('Found ', N_anchors, ' anchors')

  ## Inventorize all anchor samples in the reference
  idxs <- unique(anchors@anchors[, 'cell1'])
  print(ref_so@meta.data[idxs, ])

  predictions <- tryCatch(Seurat::TransferData(
    anchorset = anchors,
    refdata = ref_so$condition_name,
    k.weight = min(N_anchors, k.weight),
    dims = 1:N_dims),
    error = function(e) { print(e); NULL })
  if (is.null(predictions)) return(NULL)

  ann_preds <- Seurat::AddMetaData(sample_list[[query_name]],
    metadata = predictions)
  return(ann_preds)
}


#' Order an annotation data.frame
#'
#'
order_ann <- function(
  ann,
  sort_vars = c('sn_dilution', 'tnf_conc', 'ifn_conc', 'duration')) {

  sort_vars %<>% intersect(colnames(ann))
  if (length(sort_vars) == 0) return(ann)

  ann$idx <- 1:nrow(ann)
  factor_vars <-
    names(which(map_lgl(t_dat, ~(class(.x) == 'factor')[1]))) %>%
    intersect(sort_vars)
    # setdiff(c('seurat_clusters'))
  for (fac in factor_vars) {
    ann <- dplyr::mutate(ann, 
      {{fac}} := as.numeric(as.character(.data[[fac]])))
  }

  ## remove all NA columns
  ann <- ann[, map_lgl(ann, ~!all(is.na(.x)))]

  ## Determine 'logical' row ordering
  sort_vars <- intersect(sort_vars, colnames(ann))
  ann <- ann %>%
    dplyr::arrange(across(sort_vars)) %>%
    { . }

  return(ann)
}


gen_sim_measure_HM <- function(
  query_obj,
  ref_obj,
  mod,
  # merge_id = attr(query_obj, 'query_obj'),
  min_max_score = .01,
  cell_name = 'Similarity',
  cell_height = NULL,
  cell_width = NULL,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = F,
  ...) {

  # expected_cn <- c('sn_dilution', 'tnf_conc', 'ifn_conc', 'duration')
  # missing_cn <- setdiff(expected_cn, colnames(ref_ann))
  # stopifnot(length(missing_cn) < length(expected_cn))

  library(SummarizedExperiment)
  SM <- predict_stim(query_obj = query_obj, ref_obj = ref_obj, 
    mod = mod)

  # for (cn in missing_cn) {
  #   if (cn %nin% colnames(ref_ann)) {
  #     ref_ann[[cn]] <- 1L
  #   }
  # }

  if (!is.null(min_max_score)) {
    keep_ref <- apply(assays(SM)[[1]], 1, max) >= min_max_score
    print(mean(keep_ref))
    SM <- SM[keep_ref, ]
  }

  # if (is.null(merge_id)) {
  #   ## Infer merge_id
  #   overlap_count <- map_dbl(ref_ann, ~mean(.x %in% colnames(SM)))
  #   if (all(overlap_count == 0)) stop()
  #   merge_id <- names(which.max(overlap_count))
  #   stopifnot(length(setdiff(colnames(SM), ref_ann[[merge_id]])) == 0)
  # }
  # # attr(ref_obj, 'merge_id')
  # ref_ann <- 
  #   recover_ann(query = colnames(SM), lookup_data = ref_obj) %>%
  #   order_ann()
  # # attr(query_obj, 'merge_id')
  # query_ann <- 
  #   recover_ann(query = rownames(SM), lookup_data = query_obj) %>%
  #   order_ann()
  
  if (!is.null(cell_height)) {
    plot_height <- unit(cell_height * ncol(SM), 'cm')
    plot_width <- unit(cell_width * nrow(SM), 'cm')
  } else {
    plot_width <- plot_height <- NULL
  }

  row_ann_f <- 
    if (any(c('stim_group', 'condition_name') %in%
        colnames(colData(SM))))
      gen_sample_annotation_HM else function(...) NULL

  M <- t(assays(SM)[[1]])
  H <- suppressWarnings(
    gen_HM(M,
    sa = as_tibble(rowData(SM)),
    ra = as_tibble(colData(SM)),
    show_row_names = show_row_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_ann_f = row_ann_f,
    col_ann_f = gen_sample_annotation_HM,
    show_row_dend = T,
    show_column_dend = T,
    column_title = glue('Reference samples (exp\\
      {metadata(ref_obj)$experiment})'),
    row_title = glue('Query samples (exp\\
      {metadata(query_obj)$experiment})'),
    height = plot_height,
    width = plot_width,
    N_genes = NULL,
    value_name = cell_name,
    trans = identity, ...))
  attr(H, 'plot_height') <- (plot_height %||% NULL)
  attr(H, 'plot_width') <- (plot_width %||% NULL)

  print_plot_eval(
    draw(H, merge_legend = T, newpage = F,
      heatmap_legend_side = 'bottom'),
    file.path(img_dir,
      glue('RF_scores\\
        {metadata(query_obj)$experiment}.png')),
    height = attr(H, 'plot_height') %||% 12,
    width = attr(H, 'plot_width') %||% 12)

}


#' Process output of project_data
#'
#'
summarize_seurat_preds <- function(ann_preds) {
  ## Rows are query sample names, columns are reference sample names
  summary_dat <- 
    FetchData(ann_preds, c('predicted.id', 'condition_name')) %>%
    dplyr::filter(!grepl('SC digest', condition_name)) %>%
    tally(condition_name, predicted.id) %>%
    dplyr::select(-freq) %>%
    pivot_wider(
      names_from = condition_name, 
      values_from = n,
      values_fill = 0
    ) %>%
    dplyr::mutate(across(is.numeric, ~.x / sum(.x))) %>%
    { 
      cn <- .[['predicted.id']]
      set_rownames(as.matrix(dplyr::select(., where(is.numeric))), cn)
    } %>%
    t() %>%
    { . }
}


annotate_duplication_correlation <- function(so) {
  dup_settings <- 
    so@meta.data[, c('tnf_conc', 'ifn_conc', 'duration')] %>%
    dplyr::group_by(across(everything())) %>%
    dplyr::summarize(across(), rep_N = n()) %>%
    dplyr::mutate(dup = factor(ifelse(rep_N == 2, 'dup', ''))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dup_grp = 1:n())

  tmp <- 
    dplyr::left_join(
      so@meta.data, dup_settings, 
      by = c('tnf_conc', 'ifn_conc', 'duration')
      ) %>%
  set_rownames(rownames(so@meta.data))
  so@meta.data <- tmp
  rm(tmp)

  M <- as.matrix(subset_feats(GetAssayData(so), 
      read_geneset('informativeV15')))
  so@meta.data$partner <- NA_integer_
  # so@meta.data[so@meta.data$dup_grp == 1, ]
  for (i in 1:nrow(so@meta.data)) {
    if (so@meta.data$dup[i] == 'dup') {
      grp <- so@meta.data$dup_grp[i]
      partners <- which(so@meta.data$dup_grp == grp)
      if (length(partners) != 2) browser()
      for (j in partners) {
        so@meta.data$partner[j] <- setdiff(partners, j)
      }
    }
  }

  cor_ranks <- apply(cor(M, method = 'spearman'), 1, frank, ties.method = 'max')
  so@meta.data$partner_cor_nrank <- 
    map_dbl(1:ncol(so), 
      ~cor_ranks[.x, so@meta.data$partner[.x]] / ncol(so))

  return(so)
}


get_pca_embedding <- function(so) {
  M <- so %>%
    GetAssayData() %>%
    subset_feats(read_geneset('informativeV15'))

  pc <- prcomp(M, scale = T, center = T) %>% { . }
  pc_weight <- pc$sdev^2 %>% { . / sum(.) }

  pc_dtf <- pc$rotation %>%
    as.data.frame() %>%
    tibble::rownames_to_column('sample_name') %>%
    # tidyr::pivot_longer(id_cols = 'sample_name', namatches('PC'))
    tidyr::pivot_longer(cols = matches('PC'), 
      names_to = 'PC', values_to = 'PC_score') %>%
    # dplyr::left_join(
    #   tar_read_raw(glue::glue('sample_annotation_exp{experiment}')), 
    #   by = 'sample_name') %>%
    { . }

  sa <- so@meta.data %>%
    ## Force 'sample_names' to be identical to the way 'M' is named
    tibble::rownames_to_column('sn2') %>%
    dplyr::select(-any_of('sample_name')) %>%
    dplyr::rename(sample_name = sn2)

  pc_dtf <- pc_dtf %>%
    dplyr::left_join(sa, by = 'sample_name') %>%
    dplyr::mutate(
      PC = as.numeric(stringr::str_replace(PC, 'PC', ''))) %>%
    dplyr::left_join(
      tibble::enframe(pc_weight, 'PC', 'PC_weight'), 
      by = 'PC')

  return(pc_dtf)
}


identify_homeopathic_stimuli <- function(so) {
  pc_dtf <- get_pca_embedding(so)

  lfc <- function(a, b) log2(a + 1) - log2(b + 1)
  lfc <- function(a, b) log2(a) - log2(b)
  lfc <- function(a, b) log2(a + 1e-3) - log2(b + 1e-3)
  lfc <- function(a, b) log2(a + 1e-6) - log2(b + 1e-6)
  diff_from_us <- pc_dtf %>%
    ## Min-max scaling across all samples/principal components
    dplyr::mutate(PC_score = (PC_score - min(PC_score)) /
      abs(diff(range(PC_score)))) %>%
    dplyr::group_by(experiment, duration) %>%
    dplyr::summarize(
      PC,
      stim_group, 
      fc = lfc(
        PC_score, 
        PC_score[stim_group == 'Unstimulated in vitro']
      ),
      weighted_fc = fc * PC_weight
    ) %>%
    dplyr::filter(stim_group != 'Unstimulated in vitro') %>%
    { . }

  if (F) {
    diff_from_us <- diff_from_us %>%
      dplyr::group_by(duration, stim_group, experiment) %>%
      dplyr::summarize(
        N = n(),
        across(weighted_fc, 
          c('mean' = mean, 'sd'= sd, 'CI_l' = CI_l, 'CI_h' = CI_h))) %>%
      { . }
    
    diff_from_us %>%
      dplyr::ungroup() %>%
      dplyr::distinct(stim_group, experiment, duration) %>%
      dplyr::arrange(stim_group, duration, experiment)
  }

  sgo <- diff_from_us %>%
    dplyr::group_by(experiment, stim_group) %>%
    dplyr::summarize(m = sum(abs(weighted_fc))) %>%
    dplyr::arrange(m) %>%
    dplyr::ungroup() %>%
    dplyr::select(experiment, stim_group)

  diff_from_us <- dplyr::left_join(sgo, diff_from_us)

  allowed_PCs <- 
    diff_from_us %>%
    dplyr::mutate(PC = as.integer(PC)) %>%
    dplyr::group_by(PC) %>%
    dplyr::summarize(mfc = mean(abs(weighted_fc))) %>%
    dplyr::filter(mfc >= 1e-3) %>%
    pull(PC)
  max_val <- max(allowed_PCs) + 1

  out <- diff_from_us %>%
    # dplyr::mutate(PC = cut(PC, breaks = allowed_PCs)) %>%
    dplyr::mutate(PC = as.integer(PC)) %>%
    dplyr::mutate(PC = ifelse(PC >= max_val, max_val, PC)) %>%
    dplyr::mutate(PC = factor(PC, 
      levels = 1:max_val, 
      labels = c(1:max(allowed_PCs), glue::glue('>= {max_val}')))) %>%
    dplyr::group_by(experiment, PC, stim_group) %>%
    dplyr::summarize(
      across(), 
      abs_weighted_fc = mean(abs(weighted_fc))
    ) %>%
    dplyr::distinct(experiment, PC, stim_group, duration, .keep_all = T)

  out %>%
    dplyr::ungroup() %>%
    dplyr::distinct(stim_group, experiment, duration) %>%
    dplyr::arrange(stim_group, duration, experiment)

  return(out)
}
