#' DEPRECATED
#'
#'
remove_fibroblast_cluster <- function(so) {
  if (FALSE) {
    p_dat <- bind_cols(so@reductions$pca@cell.embeddings,
                       'sample_id' = rownames(so@meta.data),
                       'clusters' = so@meta.data$seurat_clusters)
    outlier_clusts <- detect_outlying_clusters(p_dat, coord_idxs = 1:50,
                                                p_val_thresh = .002)
  } else {
    p_dat <- bind_cols(so@reductions$umap@cell.embeddings,
                       'sample_id' = rownames(so@meta.data),
                       'clusters' = so@meta.data$seurat_clusters)
    # outlier_clusts <- detect_outlying_clusters(
    #   p_dat, coord_idxs = 1:2, p_val_thresh = .002)
    outlier_clusts <- detect_outlying_clusters(
      p_dat, coord_idxs = 1:2, p_val_thresh = .05)
  }

  for (oc in outlier_clusts) {
    cell_ids <- colnames(so) %>%
      { .[so@meta.data$seurat_clusters %in% oc] }
    markers <- FindMarkers(so, cell_ids)
    # markers_f <- dplyr::filter(markers,
    #   p_val_adj <= 0.25, avg_log2FC > 0)
    markers_f <- dplyr::filter(markers, p_val_adj <= 0.25)
    # marker_genes <- rownames(tail(markers, n = 50))
    # marker_genes <- rownames(head(markers, n = 50))
    detected_col_genes <-
      grep('^COL.*', rownames(markers_f), value = T)
    print(markers_f[detected_col_genes, ])

    collagenase_count <- length(detected_col_genes)
    if (F && collagenase_count > 3) {
      message('Removing cluster: ', oc)
      so <- so[, which(!so@meta.data$seurat_clusters %in% oc)]
    }

    ## If not fibroblasts, then what?
    # source(file.path(r_dir, 'plotting.R'))
    ranks <- sort(set_names(markers_f$avg_log2FC, rownames(markers_f)))
    # ranks <- setNames(seq(1, nrow(markers_f)), rownames(markers_f))
    gtab <- gsea_table(ranks)
    o_fn <-
      file.path(img_dir, glue('outlying_clust_gsea-\\
          exp{so@meta.data$exp[1]}-cluster={oc}.png'))
    ggsave(o_fn, gtab,
      # width = 17.4,
      width = 25,
      unit = 'cm')
  }

  if (F) {
    all_markers <- FindAllMarkers(so)
    all_markers %>% dplyr::filter(cluster == 2)
    all_markers <- all_markers %>%
      arrange(cluster, sign(avg_logFC), p_val)
    all_markers <- all_markers %>%
      arrange(cluster, avg_logFC, p_val)
    gtabs <-
      sort(as.integer(unique(so@meta.data$seurat_clusters))) %>%
      map(function(cl) {
        subs <- all_markers %>%
          dplyr::filter(cluster == cl) %>%
          { . }
        ranks <- setNames(seq(1, nrow(subs)), rownames(subs))
        gtab <- gsea_table(ranks, scoreType = 'pos')
      })
    # all_markers %>% group_by(cluster) %>%
    #   top_n(100) %>%
    #   print(n = 5e3L)
    ggsave(file.path(img_dir, 'tst.png'), gtab)
  }

  so <- FindVariableFeatures(so)
  return(so)
}


#' 
#'
#' @param cluster_indicator seurat_clusters or dbscan_cluster
investigate_cluster <- function(so, clust_id, contrast_cluster_ids, 
  cluster_indicator = 'seurat_clusters') {
  cell_ids <- colnames(so) %>%
    { .[unlist(so@meta.data[[cluster_indicator]]) %in% clust_id] }
  contrast_cell_ids <- colnames(so) %>%
    { .[unlist(so@meta.data[[cluster_indicator]]) %in% contrast_cluster_ids] }
  markers <- FindMarkers(so,
    ident.1 = cell_ids, ident.2 = contrast_cell_ids)
  markers <- dplyr::filter(markers, p_val_adj <= 0.25)
  return(markers)
}


markers_to_gsea <- function(
  markers,
  o_fn =
    file.path(img_dir, glue('outlying_clust_gsea-\\
        exp{so@meta.data$exp[1]}-cluster={clust_id}.png')),
  pathways = reactome_genesets) {

  ranks <- sort(set_names(markers$avg_log2FC, rownames(markers)))
  gtab <- gsea_table(ranks, pathways = pathways)
  ggsave(o_fn, gtab,
    # width = 17.4,
    # width = 25,
    width = 35,
    unit = 'cm')
  return(gtab)
}


#' Read CellMarker
#'
#'
read_CellMarker <- function(organism = 'Human') {
  # http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt
	i_fn <- file.path(data_raw_dir, glue('{organism}_cell_markers.txt'))
	fh <- readr::read_tsv(i_fn, show_col_types = F)

	# fh %>%
	# 	dplyr::select(tissueType, cellType, cellName) %>%
	# 	dplyr::group_by(tissueType, cellType, cellName) %>%
	# 	dplyr::summarize(N = n())

	out <- fh %>%
		dplyr::mutate(gs_name = glue('{tissueType}-{cellName}')) %>%
		dplyr::select(gs_name, geneSymbol) %>%
		dplyr::mutate(geneSymbol = map(geneSymbol,
				~setdiff(strsplit(.x, ', ')[[1]], NA))) %>%
		deframe()

	# outv[1:3]
	out
}


#' Test for enrichment across two variables using Fisher's exact test
#'
#'
test_OE <- function(
  dtf,
  anchor_var = 'seurat_clusters', query_var = 'exp',
  anchor_levs = NULL, query_levs = NULL) {

  if (is.null(anchor_levs)) {
    anchor_levs <- unique(dtf[[anchor_var]])
  }

  if (is.null(query_levs)) {
    query_levs <- unique(dtf[[query_var]])
  }

  ## Establish a crosstable of level combination from the two
  ## variables
  grid <- tidyr::expand_grid(anchor = anchor_levs, query = query_levs)

  fish_tests <- pmap_dfr(grid, function(anchor, query) {
    ## Binarize
    t_dat <- dtf %>%
      dplyr::mutate(
        anchor = .data[[anchor_var]] == .env[['anchor']]) %>%
      dplyr::mutate(
        query = .data[[query_var]] == .env[['query']]) %>%
      { . }
    f_table <- with(t_dat, table(anchor, query))
    if (nrow(f_table) != 2 || ncol(f_table) != 2) {
      return(tibble(estimate = NA, p.value = NA, conf.low = NA,
          conf.high = NA, method = NA, alternative = NA,
          N_dp = NA))
    }
    out <- broom::tidy(fisher.test(f_table))
    ## Number 'double positive'
    out$N_dp <- f_table[2, 2]
    return(out)
  })

  bind_cols(grid, fish_tests) %>%
    dplyr::rename({{anchor_var}} := anchor) %>%
    dplyr::rename({{query_var}} := query)
}


#' Cluster UMAP coordinates using DBSCAN density based clustering
#'
#'
identify_outlying_clusters <- function(x, ...) {
   UseMethod("identify_outlying_clusters", object = x)
   # UseMethod(generic = "identify_outlying_clusters", object = x)
}


#' Annotate a \code{data.frame} with DBSCAN clustering results based
#' on the columns
#'
#' @param min_cluster_size Minimum fold change between DBSCAN cluster
#' size and expected cluster size based on equal cluster sizes to NOT
#' be called an outlying cluster
identify_outlying_clusters.data.frame <- function(
  dtf, 
  column_selector = matches('^UMAP_*\\d+$'),
  eps = .6, min_cluster_size = .2, min_pts = 5,
  query_vars = c('percent.mt', 'N_UMI', 'exp', 
    'stim_group', 'condition_name')) {

  rn <- rownames(dtf)
  stopifnot(is.data.frame(dtf))
  if (is.data.table(dtf)) dtf <- as_tibble(dtf)
  dtf <- dplyr::select(dtf, -matches('dbscan'))
  source_data <- dtf %>%
    dplyr::select(!! enquo(column_selector))
  if (maartenutils::null_dat(source_data) || 
      ncol(source_data) < 2) {
    rlang::abort('Source data not correctly identified')
  }
  dbscan_cl <- fpc::dbscan(source_data,
    eps = eps, MinPts = min_pts)
  dtf$dbscan_cluster <- dbscan_cl$cluster

  rel_cluster_size <- table(dbscan_cl$cluster) %>%
    { . / sum(.) } %>%
    tibble::enframe('dbscan_cluster', 'dbscan_cluster_rel_size') %>%
    dplyr::mutate(dbscan_cluster = as.numeric(dbscan_cluster)) %>%
    dplyr::mutate(dbscan_cluster_rel_size = 
      as.numeric(dbscan_cluster_rel_size)) %>%
    dplyr::mutate(dbscan_cluster_rel_size_oe =
      dbscan_cluster_rel_size / (1/n())) %>%
    dplyr::mutate(outlying_dbscan_cluster = 
      dbscan_cluster_rel_size_oe <= min_cluster_size)

  dtf <-
    dplyr::left_join(dtf, rel_cluster_size, by = 'dbscan_cluster')
    
  ## Early return if no outlying clusters
  if (!any(rel_cluster_size$outlying_dbscan_cluster)) {
    rownames(dtf) <- rn
    return(dtf)
  }
  rownames(dtf) <- rn
  return(invisible(dtf))
}


identify_outlying_clusters.Seurat <- function(so, ...) {
  p_dat <- so2dtf(so)
  dtf <- identify_outlying_clusters(p_dat, ...)
  new_meta <- dtf %>% dplyr::select(matches('dbscan'))
  so <- Seurat::AddMetaData(so, new_meta)
  return(so)
}


#' Pulled out from identify_outlying_clusters.data.frame but not yet
#' used so incomplete
#'
#'
characterize_clusters <- function() {
  ## Characterize outlying clusters
  if (F) {
    regular_clusters <- rel_cluster_size %>%
      dplyr::filter(outlying_dbscan_cluster == F) %>%
      pull(dbscan_cluster)
    outlying_clusters <- rel_cluster_size %>%
      dplyr::filter(outlying_dbscan_cluster == T) %>%
      pull(dbscan_cluster)

    query_types <- dtf %>%
      dplyr::select(any_of(query_vars)) %>%
      map_chr(~class(.x)[1])

    factor_vars <- 
      query_vars[query_types %in% c('character', 'glue', 'factor')]
    numeric_vars <- 
      query_vars[query_types %in% c('numeric', 'integer')]

    ## Test which Seurat clusters correspond to the DBSCAN clusters
    seurat_overlap <- test_OE(dtf, 'dbscan_cluster', 'seurat_clusters',
      anchor_levs = outlying_clusters)
    enriched_seurat <- seurat_overlap %>% dplyr::filter(estimate > 1)

    ## What percentage of this Seurat cluster(s) is outlying?
    dtf %>%
      dplyr::filter(seurat_clusters %in%
        enriched_seurat$seurat_clusters) %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::summarize(mean(outlying_dbscan_cluster)) %>%
      print()

    ## Test for enrichment of discrete variables in outlying clusters
    t_dats <- map(factor_vars, ~test_OE(dtf, .x, 'seurat_clusters',
        query_levs = unique(enriched_seurat$seurat_clusters)))

    for (t_dat in t_dats) {
      t_dat %>%
        dplyr::filter(N_dp > 0) %>%
        dplyr::arrange(p.value) %>%
        print(width = 1000L, n = 1000L)
    }

    ## Numeric vars
    reg_idx <- dtf[['dbscan_cluster']] %in% regular_clusters
    grid <- tidyr::expand_grid(vn = numeric_vars, oc = outlying_clusters)
    out <- grid %>%
      purrr::pmap_dfr(function(vn, oc) {
        out_idx <- dtf[['dbscan_cluster']] %in% oc
        v <- dtf[[vn]]
        x <- v[reg_idx]
        y <- v[out_idx]
        if (length(x) < 20 || length(y) < 20) {
          o <- tibble(
            statistic = NA_real_, 
            p.value = NA_real_, 
            method = NA_character_, 
            alternative = NA_character_,
          )
          o[[glue('med_out')]] <- median(x, na.rm = T)
          o[[glue('med_reg')]] <- median(y, na.rm = T)
          return(o)
        }
        o <- broom::tidy(wilcox.test(x, y))
        o[[glue('med_out')]] <- median(x, na.rm = T)
        o[[glue('med_reg')]] <- median(y, na.rm = T)
        return(o)
      })
    wc_tests <- bind_cols(grid, out)
    print(wc_tests)
  }
}


so2dtf <- function(so, 
  coord_extract = function(so) {
    purrr::map(so@reductions, function(reduc) {
      reduc@cell.embeddings
    }) %>%
    purrr::reduce(cbind)
  }, meta_regex = '') {

  if (!is.null(meta_regex) && !is.na(meta_regex) && 
      meta_regex != '') {
    extra_meta <- so@meta.data %>% dplyr::select(matches(meta_regex))
  } else {
    extra_meta <- NULL
  }

  out <- dplyr::bind_cols(
    coord_extract(so),
    'sample_id' = rownames(so@meta.data),
    so@meta.data %>% 
      dplyr::select(
        matches('seurat_clusters'), 
        matches('duration|conc|rank|^sn'), 
        matches('Score|percent'), 
        matches('exp'), 
        matches('sample'), 
        matches('group'), 
        # where(is.numeric),
        matches('condition_name'), 
        matches('stim_group')),
    extra_meta
  )

  return(out)
}
