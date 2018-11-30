RunTMM <- function(x, ...) {
  UseMethod('RunTMM')
}


RunTMM.default <- function(
  object, 
  integration_feats = rownames(object),
  method = 'TMM',
  logRatioTrim = .3,
  sumTrim = 0.05,
  log = F, 
  verbose = T, 
  return_NF = F,
  ...) {
  library(edgeR)

  if (is.data.frame(x = data)) {
      object <- as.matrix(x = object)
  }
  dgC_input <- F
  if (!inherits(x = object, what = "dgCMatrix")) {
      # object <- as(object = object, Class = "dgCMatrix")
    dgC_input <- T
      object <- as.matrix(object)
  }
  if (verbose) {
      cat(paste0('Performing ', method, '-normalization\n',
          collapse = ''), file = stderr())
  }

  TMM_start_time <- proc.time()

  # Compute cell-specific scaling factors
  DGE_obj <- object[integration_feats, ] %>%
      { DGEList(counts = .) } %>%
      {
        calcNormFactors(.,
          method = method,
          logRatioTrim = logRatioTrim,
          sumTrim = sumTrim, ...)
      }
  # with(DGE_obj$samples, plot(log10(lib.size), norm.factors))
  # names(DGE_obj)
  TMM_norm_factors <- DGE_obj %>%
    { 
      apply(.@.Data[[2]][c('lib.size', 'norm.factors')], 1, 
      function(x) { 1e6 / (x[['lib.size']] * x[['norm.factors']]) }) 
    }

  if (return_NF) {
    # return(DGE_obj@.Data[[2]][['norm.factors']])
    return(TMM_norm_factors)
  }

  if (verbose) {
    cat(paste0('TMM normalization finished: ',
        data.table::timetaken(TMM_start_time)), '\n')
  }

  scale_start_time <- proc.time()
  if (F) {
    norm.data <- as.matrix(object) %*% diag(TMM_norm_factors)
  } else {
    norm.data <- sweep(object, 2, TMM_norm_factors, '*')
  }
  if (verbose) {
    cat(paste0('Scaling finished: ',
        data.table::timetaken(scale_start_time)), '\n')
  }

  if (log) {
    norm.data <- log1p(norm.data)
  }

  colnames(x = norm.data) <- colnames(x = object)
  rownames(x = norm.data) <- rownames(x = object)
  if (dgC_input) {
    norm.data  <- as(object = norm.data, Class = 'dgCMatrix')
  }

  return(norm.data)
}


RunTMM.Assay <- function(object, method = 'TMM', verbose = T, ...) {
  M <- GetAssayData(object = object, slot = 'counts')
  if (any(apply(M, 2, sum) == 0)) {
    rlang::abort('Detected samples with zero counts')
  }
  TMM_norm_data <- RunTMM(
    object = M,
    method = method, verbose = verbose, ...)
  TMM_norm_data <- as.matrix(TMM_norm_data)
  object <- SetAssayData(object = object, slot = 'data',
                         new.data = TMM_norm_data)
  object@key <- 'TMM_'
  return(object)
}


#' S3 TMM method for Seurat objects
#'
#'
RunTMM.Seurat <- function(
  object, assay = DefaultAssay(object), 
  method = 'TMM', verbose = T, ...) {

  assay <- assay %||% DefaultAssay(object = object)
  assay <- assay %||% 'RNA'
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- RunTMM(assay.data, method = method, 
    verbose = verbose, ...)
  object[['TMM']] <- assay.data

  ## Dubious but pragmatic
  object@active.assay <- 'TMM'
  object@assays$TMM@counts <- object@assays$TMM@data

  object <- LogSeuratCommand(object = object)
  return(object)
}


#'
#'
#' dfs: matrices with format [N_features, N_samples]
multiple_TMM <- function(dfs = NULL, integration_feats = NULL,
  logRatioTrim = .3, sumTrim = 0.05, ...) {
  if (is.null(dfs)) { dfs <- as.list(...) }
  NC <- map_dbl(dfs, ~ncol(.x))

  shared_feats <- map(dfs, ~rownames(.x)) %>%
    purrr::reduce(intersect) %>%
    unlist

  if (length(shared_feats) == 0) 
    stop('No shared features')
  dfs <- map(dfs, ~.x[shared_feats, ])
  if (is.null(integration_feats)) {
    integration_feats <- shared_feats
  } else {
    integration_feats <- intersect(shared_feats, integration_feats)
  }

  comb_M <- purrr::reduce(dfs, cbind)
  comb_MT <- RunTMM.default(
    object = comb_M,
    integration_feats = integration_feats,
    logRatioTrim = logRatioTrim,
    sumTrim = sumTrim, ...
  )
  dfs_out <- list()
  for (i in seq_along(NC)) {
    if (i == 1) {
      dfs_out[[i]] <- comb_MT[, 1:NC[i]]
    } else {
      i_start <- sum(NC[seq(1, i-1)]) + 1
      i_end <- sum(NC[seq(1, i)])
      dfs_out[[i]] <- comb_MT[, i_start:i_end]
    }
  }
  names(dfs_out) <- names(dfs)
  return(dfs_out)
}


