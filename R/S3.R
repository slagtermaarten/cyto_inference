extract_sc <- function(...) UseMethod('extract_sc')


extract_sc.Seurat <- function(so) {
  dtf <- so@meta.data
  dtf$i <- 1:nrow(dtf)
  dtf_s <- extract_sc(dtf)
  return(so[, dtf_s$i])
}


extract_sc.data.frame <- function(dtf) {
  if ('sample_type' %in% colnames(dtf))
    return(dtf[which(dtf$sample_type == 'sc'), ])
  else if ('exp' %in% colnames(dtf))
    return(dtf[which(!dtf$exp %in% bulk_experiments), ])
  else if ('experiment' %in% colnames(dtf))
    return(dtf[which(!dtf$experiment %in% bulk_experiments), ])
  else
    stop('Do not know how to subset extract SC samples')
}


extract_bulk <- function(...) UseMethod('extract_bulk')


extract_bulk.Seurat <- function(so) {
  dtf <- so@meta.data
  dtf$i <- 1:nrow(dtf)
  dtf_s <- extract_bulk(dtf)
  return(so[, dtf_s$i])
}


extract_bulk.data.frame <- function(dtf) {
  if ('sample_type' %in% colnames(dtf))
    return(dtf[which(dtf$sample_type == 'bulk'), ])
  else if ('exp' %in% colnames(dtf))
    return(dtf[which(dtf$exp %in% bulk_experiments), ])
  else if ('experiment' %in% colnames(dtf))
    return(dtf[which(dtf$experiment %in% bulk_experiments), ])
  else
    stop('Do not know how to subset bulk SC samples')
}


get_sc_experiments <- function(...) UseMethod('get_sc_experiments')


get_sc_experiments.data.frame <- function(dtf) {
  exp_name <- intersect(c('exp', 'experiment'), colnames(dtf))[1]
  setdiff(dtf[[exp_name]], bulk_experiments) %>%
    naturalsort::naturalsort()
}


get_in_vivo_sc_experiments <- function(...) UseMethod('get_in_vivo_sc_experiments')


get_in_vivo_sc_experiments.data.frame <- function(dtf) {
  intersect(get_sc_experiments(dtf), in_vivo_sc_experiments)
}


get_in_vitro_sc_experiments <- function(...) UseMethod('get_in_vitro_sc_experiments')


get_in_vitro_sc_experiments.data.frame <- function(dtf) {
  setdiff(get_sc_experiments(dtf), in_vivo_sc_experiments)
}


get_bulk_experiments <- function(...) UseMethod('get_bulk_experiments')


get_bulk_experiments.data.frame <- function(dtf) {
  exp_name <- intersect(c('exp', 'experiment'), colnames(dtf))[1]
  intersect(dtf[[exp_name]], bulk_experiments) %>%
    naturalsort::naturalsort()
}


strip_dimnames <- function(...) UseMethod('strip_dimnames')


strip_dimnames.matrix <- function(M) {
  rownames(M) <- NULL
  colnames(M) <- NULL
  return(M)
}


compute_gene_stats <- function(...) UseMethod('compute_gene_stats')


compute_gene_stats.Seurat <-
  function(so, genes = detected_genes(so)) {

  M <- subset_feats(so2M(so), genes)
  if (is.null(M)) return(NULL)

  stats <- compute_gene_stats(M)

  tibble::tibble(
    experiment = gen_exp_string(so@meta.data$exp),
    stats
  )
}


compute_gene_stats.matrix <- function(M) {
  if (!all(dim(M) > 0)) return(NULL)

  tibble::tibble(
    gene = rownames(M),
    min = apply(M, 1, min),
    q1 = apply(M, 1, quantile, .1),
    median = apply(M, 1, median),
    mean = apply(M, 1, mean),
    q9 = apply(M, 1, quantile, .9),
    max = apply(M, 1, max),
    q9q1om = (q9 - q1)/mean,
    sd = apply(M, 1, sd),
    CV = sd / mean,
    NULL
  )
}


CombineObject <- function(...) UseMethod('CombineObject')


CombineObject.list <- function(object_list, reduce_combine = F) {
  stopifnot(purrr::map_chr(object_list, class) == 'Seurat')
  
  object_list <- purrr::discard(object_list, is.null)

  if (reduce_combine) {
    purrr::reduce(object_list, merge, merge.data = T)
  } else {
    if (length(object_list) > 1) {
      diagnose_sample_list(object_list)
      out <- merge(
        x = object_list[[1]], 
        y = object_list[2:length(object_list)], 
        merge.data = T
      )
      diagnose_sample_list(list(out))
      rm(object_list)
      gc()
      # table(so@meta.data$exp)
      # attr(out, 'analysis_id') <- gen_exp_string(out@meta.data$exp)
      # out <- order_stim_group(out)
      return(out)
    } else {
      return(object_list[[1]])
    }
  }
}


subset_feats <- function(...) UseMethod('subset_feats')


subset_feats.Seurat <- function(
  so,
  assay = DefaultAssay(so),
  GDR_thresh = NULL,
  filter_gene_reproducibility = NULL,
  genes = NULL,
  genelist = NULL) {

  if (!is.null(genes)) {
    genes <- intersect(genes, detected_genes(so, assay = assay))
    so <- so[genes, ]
  }
  # tracemem(so)

  if (F) {
    ## Convert to list because filering steps were only implemented
    ## for list objects. Very memory intensive however.
    need_split <-
      c(GDR_thresh, filter_gene_reproducibility, genelist) %>%
      purrr::map_lgl(~!is.null(.x)) %>%
      any(na.rm = T)
    if (need_split) {
      object_list <- SplitObject(so, 'exp')
      rm(so)
      object_list <- subset_feats(object_list,
        GDR_thresh = GDR_thresh,
        filter_gene_reproducibility = filter_gene_reproducibility,
        genelist = genelist
      )
      ## Return back to Seurat object
      so <- CombineObject(object_list)
      rm(object_list)
    }
  } else {
    if (!is.null(genelist)) {
      if (stringr::str_detect(genelist, '^vst')) {
        N_genes <- stringr::str_replace(genelist, 'vst', '') %>%
          as.numeric()
        so <- Seurat::FindVariableFeatures(
          so,
          nfeatures = N_genes,
          selection.method = 'vst')
        feats <- VariableFeatures(so)
        # var_feats <- sort(table(feats), decreasing = T)[1:N_genes]
        # var_feats <- names(var_feats)
        so <- so[feats, ]
      } else if (!genelist %in% c('', 'all')) {
        genelist_genes <- read_geneset(glue::glue('{genelist}_genes'))
        so <- subset_feats(so, genes = genelist_genes)
      }
    }

    if (!is.null(filter_gene_reproducibility)) {
      rep_genes <-
        filter_gene_reproducibility_score(filter_gene_reproducibility)
      so <- subset_feats(so, genes = rep_genes)
    }

    ## Apply GDR filtering; require gene to be detected in at least x
    ## % of cells in one or more conditions
    if (!is.null(GDR_thresh) && GDR_thresh > 0) {
      # GDR_thresh = .25
      # source('~/MirjamHoekstra/R/init.R')
      detectable_genes <- get_GDR_table(
        extract_sc(so), thresh = GDR_thresh,
        group_var = 'condition_name')
      so <- subset_feats(so, genes = detectable_genes)
    }
  }
  return(so)
}
if (F) {
  if (!exists('so'))
    so <- tar_read(filtered_cleaned_so, branch = 1)[[1]]
  subset_feats(so,
    genelist = 'vst600',
    filter_gene_reproducibility = 2L)
}


subset_feats.Assay <- function(
  x,
  genelist = NULL,
  genes = NULL,
  filter_gene_reproducibility = NULL,
  GDR_thresh = NULL) {

  shared_features <- rownames(GetAssayData(x))

  if (!is.null(genes) && is.vector(genes) && length(genes) > 0) {
    shared_features <- intersect(shared_features, genes)
  }

  if (!is.null(genelist)) {
    if (stringr::str_detect(genelist, '^vst')) {
      N_genes <- as.numeric(stringr::str_replace(genelist, 'vst', ''))
      x <- Seurat::FindVariableFeatures(x,
        selection.method = 'vst', nfeatures = N_genes)
      shared_features <- intersect(
        shared_features, VariableFeatures(x)
      )
    } else if (!genelist %in% c('', 'all')) {
      genelist_genes <- read_geneset(glue::glue('{genelist}_genes'))
      shared_features <- intersect(shared_features, genelist_genes)
    }
  }

  if (!is.null(filter_gene_reproducibility)) {
    shared_features <- intersect(
      shared_features,
      filter_gene_reproducibility_score(filter_gene_reproducibility)
    )
  }

  ## Apply GDR filtering; require gene to be detected in at least x %
  ## of cells in one or more conditions
  if (!is.null(GDR_thresh) && is.numeric(GDR_thresh) &&
      GDR_thresh > 0) {
    detected_genes <-
      get_GDR_table(
        x, thresh = GDR_thresh,
        genes = shared_features) %>%
      unique()

    shared_features <- intersect(shared_features, detected_genes)
  }

  x <- x[shared_features, ]

  return(x)
}


if (F) {
  subset_feats.list <- function(
  sample_list,
  genes = NULL) {
  # genelist = NULL,
  # GDR_thresh = NULL,
  # filter_gene_reproducibility = NULL) {

  if (is.null(sample_list) || length(sample_list) == 0) return(NULL)

  shared_features <-
    purrr::map(sample_list, ~detected_genes(.x)) %>%
    purrr::reduce(intersect)

  if (!is.null(genes))
    shared_features <- intersect(shared_features, genes)

  # if (!is.null(genelist)) {
  #   if (stringr::str_detect(genelist, '^vst')) {
  #     ## Slightly better could be to concatenate and normalize all
  #     ## data and then do VST on that. But it's too much effort right
  #     ## now for the dubious benefits I expect of it.
  #     # Ms <- purrr::map(sample_list, ~so2M(.x))
  #     # local_shared_features <- map(Ms, ~rownames(.x)) %>%
  #     #   reduce(intersect)
  #     # Ms <- purrr::map(Ms, ~.x[local_shared_features, ])
  #     # M <- purrr::map_dfc(Ms, ~as.data.frame(.x))
  #     sample_list <- purrr::map(sample_list,
  #       ~Seurat::FindVariableFeatures(.x, selection.method = 'vst'))
  #     feats <- purrr::map(sample_list, ~VariableFeatures(.x)) %>%
  #       rlang::flatten_chr() %>%
  #       intersect(shared_features)
  #     N_genes <- as.numeric(stringr::str_replace(genelist, 'vst', ''))
  #     var_feats <- sort(table(feats), decreasing = T)[1:N_genes]
  #     var_feats <- names(var_feats)
  #     shared_features <- intersect(shared_features, var_feats)
  #   } else if (!genelist %in% c('', 'all')) {
  #     genelist_genes <- read_geneset(glue::glue('{genelist}_genes'))
  #     shared_features <- intersect(shared_features, genelist_genes)
  #   }
  # }

  # if (!is.null(filter_gene_reproducibility)) {
  #   shared_features <- intersect(
  #     shared_features,
  #     filter_gene_reproducibility_score(filter_gene_reproducibility)
  #   )
  # }

  # ## Apply GDR filtering; require gene to be detected in at least x %
  # ## of cells in one or more conditions
  # if (!is.null(GDR_thresh)) {
  #   detected_genes <-
  #     purrr::map(sample_list, get_GDR_table, thresh = GDR_thresh,
  #       genes = shared_features) %>%
  #     rlang::flatten_chr() %>%
  #     unique()

  #   shared_features <- intersect(shared_features, detected_genes)
  # }

  sample_list <- purrr::map(sample_list, ~.x[shared_features, ])

  return(sample_list)
}
}


subset_feats.matrix <- subset_feats.dgCMatrix <- function(M, genes) {

  if (is.null(genes) || length(genes) == 0) return(M)
  dim_idx <- match_dimnames(M, genes)

  if (dim_idx == 1) {
    sub_feats <- intersect(genes, rownames(M))
    if (length(sub_feats) == 0) {
      rlang::warn('None of the requested features are available')
    } else {
      M <- M[sub_feats, , drop=F]
    }
  } else if (dim_idx == 2) {
    sub_feats <- intersect(genes, colnames(M))
    if (length(sub_feats) == 0) {
      rlang::warn('None of the requested features are available')
    } else {
      M <- M[, sub_feats, drop=F]
    }
  } else {
    rlang::abort('Unforeseen case')
  }

  return(M)
}


compute_CDR <- function(...) UseMethod('compute_CDR')


compute_CDR.matrix <- compute_CDR.dgCMatrix <- function(M) {
  idx <- setdiff(1:2, 
    match_dimnames(M, c('IDO1', 'CD247', 'TNFAIP9', 'WARS')))
  return(apply(M, idx, function(x) mean(x > 0)))
}


compute_CDR.Seurat <- function(so, assay = DefaultAssay(so)) {
  compute_CDR(GetAssayData(so, assay = assay))
}


compute_GDR <- function(...) UseMethod('compute_GDR')


compute_GDR.matrix <- compute_GDR.dgCMatrix <- function(M) {
  idx <- match_dimnames(M, c('IDO1', 'CD247', 'TNFAIP9', 'WARS'))
  return(apply(M, idx, function(x) mean(x > 0)))
}


compute_GDR.numeric <- function(v) {
  return(mean(v > 0))
}


compute_GDR.Seurat <- function(so, assay = DefaultAssay(so)) {
  compute_GDR(GetAssayData(so, assay = assay))
}


compute_GDR.Assay <- function(x) {
  compute_GDR(x[,])
}


detected_genes <- function(...) UseMethod('detected_genes')


detected_genes.Seurat <- function(so, assay = DefaultAssay(so)) {
  # rownames(so@assays[[assay]])
  rownames(so@assays[[assay]])
}


detected_genes.Assay <- function(x) {
  rownames(x)
}


detected_genes.matrix <- detected_genes.dgCMatrix <- function(x) {
  rownames(x)
}


group_wise_GDR <- function(...) UseMethod('group_wise_GDR')


group_wise_GDR.Seurat <- function(
  so,
  genes = detected_genes(so),
  assay = DefaultAssay(so),
  group_var = 'condition_name') {

  if (!is.null(group_var))
    u_var <- so@meta.data[[group_var]]
  else 
    u_var <- NULL
  stopifnot(length(u_var) == 1L)
  group_wise_GDR(x = GetAssay(so, assay = assay), u_var = u_var)
}


group_wise_GDR.Assay <- function(
  x, u_var = NULL,
  genes = detected_genes(x)) {

  stopifnot(is.null(u_var) || length(u_var) == ncol(x))
  M <- GetAssayData(x)
  group_wise_GDR(M = M, u_var = u_var)
}


group_wise_GDR.matrix <- group_wise_GDR.dgCMatrix <- function(
  M, u_var,
  genes = detected_genes(M)) {

  stopifnot(is.null(u_var) || length(u_var) == ncol(x))
  M <- subset_feats(M, genes)
  if (is.null(u_var)) u_var <- rep('a', ncol(M))
  GDRs <- apply(M, 1,
    function(x) tapply(x, u_var, compute_GDR))
  return(GDRs)
}


if (F) {
  if (!exists(so))
    so <- tar_read(filtered_cleaned_so, branch = 1)[[1]]
  group_wise_GDR(so)
}


get_GDR_table <- function(...) UseMethod('get_GDR_table')


get_GDR_table.Seurat <- function(
  so,
  group_var = 'condition_name',
  thresh = NULL,
  format = 'passing_genes',
  genes = VariableFeatures(so) %|||% detected_genes(so),
  assay = DefaultAssay(so)) {

  if (!is.null(so@meta.data[[group_var]]) && 
      !all(is.na(so@meta.data[[group_var]]))) {
    u_var <- so@meta.data[[group_var]]
  } else {
    u_var <- NULL
  }
  out <- get_GDR_table(
    subset_feats(GetAssay(so, assay = assay), genes = genes),
    u_var = u_var,
    thresh = thresh,
    genes = genes,
    format = format
  )

  return(out)
}


get_GDR_table.Assay <- function(
  x, 
  u_var,
  thresh = NULL,
  format = 'passing_genes',
  genes = detected_genes(x)) {

  get_GDR_table(M = x[,], 
    u_var = u_var, thresh = thresh,
    format = format, genes = genes)
}


get_GDR_table.matrix <- get_GDR_table.dgCMatrix <- function(
  M, u_var, thresh = NULL,
  format = 'passing_genes', 
  genes = detected_genes(M)) {
  
  # if (is.null(M) || any(dim(M) == 0)) return(NULL)

  if ((is.null(thresh) || thresh <= 0) && format == 'passing_genes') {
    return(genes)
  }

  if (nrow(M) == 0 && format == 'passing_genes') return(c())

  gwGDR <- group_wise_GDR(M, u_var = u_var, genes = genes)
  if (is.null(u_var)) {
    max_GDR <- gwGDR
  } else {
    max_GDR <- apply(gwGDR, 2, max, na.rm = T)
  }

  GDR_table <- tibble(
    gene = detected_genes(M),
    global = compute_GDR(M) >= thresh,
    group_wise = max_GDR >= thresh
  )

  if (format == 'table') {
    return(GDR_table)
  } else if (format == 'passing_genes') {
    passing_genes <-
      dplyr::filter(GDR_table, global == T | group_wise == T) %>%
      dplyr::pull(gene)
    return(passing_genes)
  }
}
