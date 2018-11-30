#' Annotate 'condition' (stim_group) for combicoded experiments
#'
#' @param so A Seurat object that has a HTO assay included.
#' @param hash_tag_annotation data.frame object of dimensions [N_cells x
#' (2+N_condition_id_vars)] in which all
#' expected hash tag combinations and their associated conditions are
#' listed with one row for each condition. For example:
#' tribble(
#'   ~HT1, ~HT2, ~condition,
#'   1, 2, 'cond1',
#'   1, 3, 'cond2',
#'   1, 4, 'cond3')
#' , in which the 'winning' combination of HTOs 1 and 2 will be
#' assigned to condition number 1 ('cond1'). This data.frame can have
#' any number of additional relevant variables, which will be copied
#' to the output.
#'
#' @return A seurat object, in which the meta data has been augmented
#' with here computed HTO stats. We advise subsequently filtering
#' cells on the fd_cc and hashtag_evenness variables
#'
annotate_double_barcode_stats <- function(
  so, method = 'raw', ncores = 8, debug = F,
  hash_tag_annotation = tribble(
    ~HT1, ~HT2, ~sample_name,
    1, 2, 'IFNy 2h 0 ng/ml',
    1, 3, 'IFNy 2h 0.01 ng/ml',
    1, 4, 'IFNy 2h 1 ng/ml',
    1, 6, 'IFNy 2h 100 ng/ml',
    2, 4, 'IFNy 6h 0 ng/ml',
    2, 6, 'IFNy 6h 0.01 ng/ml',
    2, 8, 'IFNy 6h 1 ng/ml',
    3, 4, 'IFNy 6h 100 ng/ml',
    3, 5, 'IFNy 12h 0 ng/ml',
    3, 7, 'IFNy 12h 0.01 ng/ml',
    4, 7, 'IFNy 12h 1 ng/ml',
    5, 6, 'IFNy 12h 100 ng/ml',
    5, 7, 'IFNy 24h 0 ng/ml',
    5, 8, 'IFNy 24h 0.01 ng/ml',
    6, 8, 'IFNy 24h 1 ng/ml',
    7, 8, 'IFNy 24h 100 ng/ml')) {
  library(data.table)
  library(Seurat)

  stopifnot(is.data.frame(hash_tag_annotation))
  hash_tag_annotation <- hash_tag_annotation %>%
    data.table::setDT() %>%
    ## This ID indicates the order with which the conditions were
    ## ordered originally. setkey() might change the ordering
    .[, 'condition_i' := 1:.N] %>%
    setkeyv(c('HT1', 'HT2'))

  ## Make sure HT1 < HT2
  hash_tag_annotation[as.integer(HT1) > as.integer(HT2),
    c('HT1', 'HT2') := list(HT2, HT1)]

  if (nrow(unique(hash_tag_annotation, by = c('HT1', 'HT2'))) !=
      nrow(hash_tag_annotation)) {
    stop(
      glue::glue('One or more combinations in hash_tag_annotation\\
        are duplicated'))
  }

  ## This is redundant but I keep it in just to be 100% sure
  hash_tag_annotation <-
    hash_tag_annotation[order(as.integer(HT1), as.integer(HT2))]

  ## Save rownames as they are crucial for Seurat
  rn <- rownames(so@meta.data)
  if (method == 'raw') {
    M <- log(t(as.matrix(GetAssayData(so[['HTO']]))) + 1)
  } else if (method == 'CLR') {
    if (!'NormalizeData.HTO' %in% names(so@commands)) {
      so <- Seurat::NormalizeData(
        so, assay = 'HTO',
        normalization.method = 'CLR')
    }
    so <- HTO_demux_titration(so)
    M <- t(so[['HTO']]@data)
    stopifnot(all(M >= 0))
  }
  stopifnot(all(rownames(M) == rn))

  expected_barcodes <-
    unique(unname(unlist(
          c(hash_tag_annotation[, 'HT1'],
            hash_tag_annotation[, 'HT2'])))) %>%
    as.integer() %>% sort()
  observed_barcodes <-
    stringr::str_replace(colnames(M), '[^\\d]*(\\d+)', '\\1') %>%
    as.integer()
  col_idxs <- setdiff(match(expected_barcodes, observed_barcodes), NA)
  M <- M[, col_idxs]

  if (any(apply(M, 2, function(x) all(is.na(x))))) browser()

  ## Offset the HT names such that the first hashtag also corresponds
  ## to the first column of M
  if (min(expected_barcodes) > 1L) {
    d <- min(expected_barcodes) - 1L
    hash_tag_annotation$HT1 <- hash_tag_annotation$HT1 - d
    hash_tag_annotation$HT2 <- hash_tag_annotation$HT2 - d
  } else {
    d <- 0
  }

  ## Ensure barcodes are sequential, current code cannot deal with
  ## 'gaps' in the HTOs
  stopifnot(all(expected_barcodes == 
      seq(min(expected_barcodes), max(expected_barcodes))))

  if (debug) {
    ## Save HTO so that we can bypass all preceeding steps when
    ## debugging
    saveRDS(M, file.path(Sys.getenv('rds_dir'), 
        glue('norm_HTO_{Project(so)}.rds')))
    # gen_file_overview(Sys.getenv('rds_dir'), pat = '.*norm_HTO.*')
    # readRDS(file.path(Sys.getenv('rds_dir'), glue('norm_HTO_.rds')))
  }

  sc_HT_diagnostics <- 
    find_matching_double_barcodes(
      M = M, 
      hash_tag_annotation = hash_tag_annotation, 
      debug = debug,
      ncores = ifelse(debug, 1, ncores))

  if (any(duplicated(sc_HT_diagnostics$X1))) {
    sc_HT_diagnostics <-
      sc_HT_diagnostics[-which(duplicated(sc_HT_diagnostics$X1)), ]
  }

  ## Correct d again for easier future reference to sample annotation
  sc_HT_diagnostics$HT1 <- as.integer(sc_HT_diagnostics$HT1) + d
  sc_HT_diagnostics$HT2 <- as.integer(sc_HT_diagnostics$HT2) + d

  extra_meta <- sc_HT_diagnostics %>%
    dplyr::select(-matches('HTO')) %>%
    # magrittr::column_to_rownames('X1') %>%
    # dplyr::rename_with(~paste0('HTO_', method, '_', .x))
    { . }

  stopifnot(all(extra_meta$X1 == rownames(so@meta.data)))
  rownames(extra_meta) <- extra_meta$X1
  extra_meta$X1 <- NULL
  so <- AddMetaData(so, extra_meta)

  return(so)
}


find_matching_double_barcodes <- function(M, hash_tag_annotation, 
  ncores = 1, debug = F) {
  library(data.table)

  all_combs <- 
    hash_tag_annotation %>%
    dplyr::select(HT1, HT2) %>%
    unlist %>% 
    unique %>%
    { tidyr::expand_grid(HT1 = ., HT2 = .) } %>%
    setDT %>%
    .[as.integer(HT1) < as.integer(HT2)] %>%
    setkeyv(c('HT1', 'HT2')) %>%
    maartenutils::controlled_merge(
      hash_tag_annotation, by_cols = c('HT1', 'HT2')) %>%
    .[hash_tag_annotation, 'expected' := TRUE] %>%
    .[is.na(expected), 'expected' := FALSE] %>%
    { . }

  stopifnot(sum(all_combs$expected) == nrow(hash_tag_annotation))

  if (ncores > 1) {
    doParallel::registerDoParallel(cores = ncores)
  } else {
    # doParallel::registerDoParallel(cores = 1)
  }

  ## Iterate over all rows of M, corresponding to one cell each, and
  ## determine the likelikehood of each potential hash tag combination
  ## to be true for each cell
  sc_HT_diagnostics <- M %>%
    plyr::adply(1, function(x) {
      if (debug) browser()
      ## Hash tag products
      if (any(is.na(x))) browser()
      if (any(x < 0)) {
        stop('I am not prepared to deal with negative HTO values')
      }
      HTO_scores <- plyr::adply(all_combs, 1, function(i) {
        idxs <- as.integer(unlist(i[, c('HT1', 'HT2'), with = F]))
        vals <- x[idxs]
        c('comb_score' = prod(vals))
      }) %>% setDT
      setkeyv(HTO_scores, c('HT1', 'HT2'))
      HTO_scores_e <- HTO_scores[expected == T]
      ## nrow == 1 if match, else 0
      if (HTO_scores_e[
        comb_score == max(comb_score, na.rm = T), .N] == 1) {
        ## Fold difference closest contenter: the fold-difference of
        ## the winning legal combination and the second-best legal
        ## combination
        fd_cc <- HTO_scores_e[,
          sort(max(comb_score, na.rm = T)/comb_score)[2]]
      } else {
        fd_cc <- NA
      }
      ## Winning HTO combination among ALL potential combinations
      best_match <- HTO_scores[
        comb_score == max(comb_score, na.rm = T)] %>%
        # { .[,3:ncol(.)] }
        { . }
      best_match$bm_exp <- best_match$expected
      best_match$expected <- NULL
      best_match$comb_score <- NULL
      ret_val <- tibble::tibble(
        best_match,
        'HTO_sum' = sum(x),
        'hashtag_evenness' =
          maartenutils::compute_evenness(
            setdiff(HTO_scores$comb_score, NA)),
        'hashtag_evenness_exp' =
          maartenutils::compute_evenness(
            setdiff(HTO_scores_e$comb_score, NA)),
        'hashtag_evenness_unexp' =
          HTO_scores[expected == F, comb_score] %>%
            setdiff(NA) %>%
            maartenutils::compute_evenness(),
        'fd_cc' = fd_cc
      )
      if (length(ret_val) != 6+(ncol(hash_tag_annotation)) && debug)
        browser()
      return(ret_val)
    }, .parallel = (ncores > 1))

  return(sc_HT_diagnostics)
}


#' Deprecated/experimental
#'
#'
remove_doublets <- function(so, N_pcs = 20, resolution = 2.0) {
  # so <- subset(so, subset = nFeature_RNA > 200 & percent.mt < 40)
  # so <- NormalizeData(so, normalization.method = "LogNormalize",
													# scale.factor = 10000)
  so <- FindVariableFeatures(object = so)
  # so <- ScaleData(so)
	# top10 <- head(VariableFeatures(so), 10)
  so <- RunPCA(so, features = VariableFeatures(object = so))
  # so@meta.data$orig.ident <- so@meta.data$dominant_hashtag
	# elbo_plot <- ElbowPlot(so)
	so <- FindNeighbors(so, dims = 1:N_pcs)
	# so <- FindClusters(so, resolution = 1.5)
	# so <- FindClusters(so, resolution = 2.0)
	# so <- FindClusters(so, resolution = 3.0)
	# so <- FindClusters(so, resolution = 4.0)
	so <- FindClusters(so, resolution = resolution, algorithm = 'Leiden')
  message(sprintf('Found %d clusters', max(so@meta.data$seurat_clusters)))
	# so@meta.data$`RNA_snn_res.1.5`
	# so@meta.data$`RNA_snn_res.2`
	# so@meta.data$`RNA_snn_res.3`
	so@meta.data$seurat_clusters <- so@meta.data$seurat_clusters - 1
	# so@meta.data$`RNA_snn_res.5`
  # sc_obj.markers = FindAllMarkers(object = so,
																	# only.pos = TRUE, min.pct = 0.25,
																	# logfc.threshold = 0.25)

  if (!require('DoubletDecon')) {
    devtools::install_github('EDePasquale/DoubletDecon')
  }
  library(DoubletDecon)
  so <- Improved_Seurat_Pre_Process(so, num_genes = 50,
                                        write_files = FALSE)
	# so <- RunUMAP(so, dims = 1:N_pcs)
  return(so)
}


#' Use HTO Demux and use the results of whatever clustering
#' algorithm will not result in an error
#'
HTO_demux_titration <- function(so) {
  library(Seurat)

  param_grid <- tidyr::expand_grid(
    positive.quantile = seq(.99, .75, by = -.01),
    kfunc = c('clara', 'kmeans'))

  for (i in 1:nrow(param_grid)) {
    message(i)
    demux_res <- tryCatch(
      with(param_grid[i, ],
        HTODemux(so, assay = 'HTO',
          kfunc = kfunc, positive.quantile = positive.quantile)),
      error = function(e) { NULL })
    if (!is.null(demux_res)) break
  }

  if (is.null(demux_res)) {
    warning('No reasonable HTO demux setting could be found')
  } else {
    so <- demux_res
    attr(so, 'kfunc') <- param_grid[i, ]$'kfunc'
    attr(so, 'positive.quantile') <-
      param_grid[i, ]$'positive.quantile'
  }
  return(so)
}
