## Default filtering settings
default_filtering_opts <- list(
  hashtag = 'HTO_classification',
  clustering_resolution = 2,
  # max_percent_mt = 100,
  max_percent_mt = 30,
  min_UMI = 1000,
  ## Only relevant for doubly hash tagged experiments
  min_fd_cc = 2,
  max_hashtag_evenness = 1
)

old_filtering_opts <- default_filtering_opts %>%
  replace('max_percent_mt', 100)


#' Compute the log fold difference between the winning hash tag and
#' the second best hash tag
#'
compute_fdcc <- function(x, method = 'old') {
  if (is.null(x) || data.table::uniqueN(x) == 2) {
    return(NA)
  }
  if (method == 'new') {
    x_rank <- rank(x, ties.method = 'max')
    max_entry <-
      x[which(x_rank == (sort(x_rank) %>% { .[length(.)] }))[1]]
    sec_best <-
      x[which(x_rank == (sort(x_rank) %>% { .[length(.) - 1] }))[1]]
  } else if (method == 'old') {
    x_order <- order(x)
    max_entry <- x[dplyr::last(x_order)]
    sec_best <- x[x_order[length(x) - 1]]
  }
  return(log2(max_entry + 1) - log2(sec_best + 1))
}
stopifnot(compute_fdcc(base::sample(c(200, 9, 100, 4)), method = 'old') ==
          log2(200+1)-log2(100+1))
stopifnot(compute_fdcc(base::sample(c(100, 9, 100, 4)), method = 'old') ==
          log2(100+1)-log2(100+1))


create_seurat_object <- function(
  data_dir = file.path(data_root, 'raw_exp_5310'),
  gene_column = 2, min_cells = 3,
  obj_id = gsub('raw_exp_', '', basename(data_dir)),
  save_res = F,
  min_hashtag_counts = 0,
  min_features = 1000) {
  library(Seurat)

  data_obj <- Read10X(
    data.dir = data_dir, gene.column = gene_column, strip.suffix = T)
  ht_i <- length(data_obj)
  hashtag_counts <- data_obj[[ht_i]]

  cell_metadat <-
    data.frame(
      fd_cc = apply(data_obj[[ht_i]], 2, compute_fdcc),
      hashtag_evenness = apply(data_obj[[ht_i]], 2, 
        maartenutils::compute_evenness),
      dominant_hashtag = apply(data_obj[[ht_i]], 2, which.max)
    ) %>%
    cbind(t(as.data.frame(data_obj[[ht_i]]))) %>%
    { magrittr::set_rownames(., paste0(obj_id, '_', rownames(.))) }

  if (ht_i > 2) {
    other_metas <- purrr::map_dfc(2:(ht_i-1), ~as.data.frame(data_obj[[.x]]))
    cell_metadat <- cbind(cell_metadat, t(other_metas))
    cell_metadat$CD123 <- NULL
  }
  counts_M <- data_obj[[1]]
  colnames(counts_M) <- paste0(obj_id, '_', colnames(counts_M))
  colnames(hashtag_counts) <- paste0(obj_id, '_', colnames(hashtag_counts))
  so <- CreateSeuratObject(counts = counts_M,
    project = obj_id,
    min.cells = min_cells,
    min.features = min_features,
    meta.data = cell_metadat
  )

  Project(so) <- obj_id
  so@meta.data$exp <- obj_id

  # min_hashtag_counts = 500
  ## Ensure at least one barcode read is available
  valid_cells <- names(which(apply(hashtag_counts, 2, sum) >
      min_hashtag_counts)) %>%
    intersect(colnames(so))
  so <- so[, valid_cells]
  hashtag_counts <- hashtag_counts[, valid_cells]

  if (T) {
    # allowed_HTO <- GetAssayData(so[['HTO']], 'counts') %>%
    allowed_HTO <- hashtag_counts %>%
      { apply(., 1, sum) > ncol(so) } %>%
      which %>% names
    hashtag_counts <- hashtag_counts[allowed_HTO, ]
  }

  so[['HTO']] <-
    suppressWarnings(CreateAssayObject(counts = hashtag_counts))
  so <- NormalizeData(so, assay = 'HTO',
    normalization.method = 'CLR')

  # M <- so@assays[['HTO']][,]
  # all(apply(M, 2, sum) > 0)
  # so <- ScaleData(so, assay = 'HTO', display.progress = TRUE,
  #                     features = rownames(so@assays$HTO@data))

  if (save_res) {
    o_fn <- file.path(rds_dir, glue::glue('exp{obj_id}_sc.rds'))
    message(glue::glue('Saving to {o_fn}'))
    saveRDS(so, o_fn)
    return(invisible(o_fn))
  } else {
    return(invisible(so))
  }
}


seurat_process_HTO <- function(
  so, ncores = 32,
  sa = NULL,
  debug_annotate_double_barcode_stats = F) { 
  source(file.path(r_dir, 'HTO_processing.R'))

  ## Salvage sample annotation from targets, not that clean
  if (is.null(sa)) {
    ## 2022-03-20 14:01 Having tar_read_raw in a function probably
    ## leads it to be permanently marked as out-of-date by targets(),
    ## so I no longer support this 'hackish' behaviour
    # sa <- targets::tar_read_raw(
    #   glue::glue('sc_{experiment}_sample_annotation'))
    if (is.null(sa))
      stop(glue::glue('No sample annotation found for exp {experiment}'))
  }

  ## Meta-data updating
  bak <- so@meta.data
  so@meta.data <- so@meta.data %>% 
    dplyr::select(-any_of(colnames(sa)))
  barcode_mode <- if (all(c('HT1', 'HT2') %in% colnames(sa)))
    'double' else 'single'
  cell_barcodes <- rownames(so@meta.data)

  if (barcode_mode == 'single') {
    so <- HTO_demux_titration(so)
    if (any(grepl('human-Hashtag', so@meta.data$HTO_classification))) {
      so@meta.data$HTO_classification <-
        as.numeric(gsub('.*(\\d+)$', '\\1',
                        so@meta.data$HTO_classification))
    } else if (any(grepl('^HTO', so@meta.data$HTO_classification))) {
      so@meta.data$HTO_classification <-
        as.numeric(gsub('^HTO(\\d+)-.*$', '\\1',
                        so@meta.data$HTO_classification))
        # v <- 'HTO2-CFP-HBSS-2_HTO4-CFP-IFNg-2'
        # gsub('HTO(\\d{1,2})-.*', '\\1', v, perl = TRUE)
    } else {
      browser()
    }

    lead_string <- gsub('\\d+$', '', so@meta.data[['HTO_maxID']][1])
    so@meta.data$dominant_hashtag <-
    so@meta.data$condition_i <-
      as.integer(gsub(lead_string, '', so@meta.data[['HTO_maxID']]))
    stopifnot(mean(!is.na(so@meta.data$dominant_hashtag)) > .75)
    new_meta <- 
      sa[match(so@meta.data$dominant_hashtag, sa$hash_tag), ]
  } else if (barcode_mode == 'double') {
    hash_tag_annotation <- sa %>%
      dplyr::select(HT1, HT2, stim_group, 
        duration, any_of('mouse')) %>%
      { . }

    so_u <- annotate_double_barcode_stats(
      so, hash_tag_annotation = hash_tag_annotation,
      ncores = ncores, debug = debug_annotate_double_barcode_stats,
      method = 'CLR'
    )
    # so_u <- so
    # so_u@meta.data[['condition_i']] <- 1
    # head(so@meta.data)
    # head(so_u@meta.data)
    # so_u@meta.data
    so <- so_u; rm(so_u)
    so@meta.data$dominant_hashtag <- so@meta.data[['condition_i']]
    ## Condition_i recorded the index of the best fitting condition
    ## based on the original ordering, which may have been altered
    ## within annotate_double_barcode_stats()
    new_meta <- sa[so@meta.data$condition_i, ]
  }

  new_meta <- dplyr::mutate(new_meta,
    group_id = sprintf('%s-%s', sample_type, sample_origin))
  stopifnot(nrow(new_meta) == nrow(so@meta.data))
  new_meta <- 
    cbind(so@meta.data, new_meta) %>%
    { .[, !duplicated(colnames(.))] } %>%
    { magrittr::set_rownames(., cell_barcodes) }

  so@meta.data <- new_meta
  head(so@meta.data)
  so@meta.data[['hashtag']] <- so@meta.data[['dominant_hashtag']]
  so@meta.data[['N_UMI']] <- so@meta.data[['nCount_RNA']]

  so <- process_so(so, sct = 'vanilla')

  return(so)
}


#' Deprecated
#'
#'
prep_all_so <- function(
  experiment = '5310',
  filtering_opts = default_filtering_opts,
  prep_filtered = F,
  prep_cleaned = F,
  genelist = '5092_5310_shared_informativeV11',
  prep_mito = F,
  so_rc = maartenutils::gen_time_comparator('2021-09-11 09:30',
    verbose = !parallel),
  redo = F, ncores = 32) {

  fns <- compile_so_fns(
    experiment = experiment,
    filtering_opts = filtering_opts)

  if (so_rc(fns[['vanilla']]) || redo) {
    so <- create_vanilla_so(experiment = experiment)
    saveRDS(so, fns[['vanilla']])
  }

  if (so_rc(fns[['HTO_QC']]) || redo) {
    so <- readRDS(fns[['vanilla']])
    so_f <- filter_HQ_cells(
      so = so,
      hashtag = filtering_opts[['hashtag']],
      max_percent_mt = 100,
      min_UMI = 0,
      min_fd_cc = filtering_opts[['min_fd_cc']],
      max_hashtag_evenness = filtering_opts[['max_hashtag_evenness']]
    )
    so_f <- process_so(so_f, sct = 'vanilla')
    saveRDS(so_f, fns[['HTO_QC']])
  }

  if (prep_filtered || prep_cleaned || prep_mito) {
    if (so_rc(fns[['filtered']]) || redo) {
      so <- readRDS(fns[['HTO_QC']])
      # source('~/MirjamHoekstra/R/init.R')
      # summary(so@meta.data$N_UMI)
      # mean(so@meta.data$N_UMI < filtering_opts[['min_UMI']])
      so_f <- filter_HQ_cells(
          so = so,
          hashtag = filtering_opts[['hashtag']],
          max_percent_mt = filtering_opts[['max_percent_mt']],
          min_UMI = filtering_opts[['min_UMI']])
      so_f <- process_so(so_f, sct = 'vanilla')
      saveRDS(so_f, fns[['filtered']])
    }

    if (so_rc(fns[['filtered_ig']]) || redo) {
      geneset <- read_geneset(glue::glue('{genelist}_genes'))
      so_f <- readRDS(fns[['filtered']])
      so_f_ig <- process_so(so_f[geneset, ], redo_SCT = F)
      saveRDS(so_f_ig, fns[['filtered_ig']])
    }

    if (so_rc(fns[['filtered_ig_mono']]) || redo) {
      geneset <- read_geneset(glue::glue('{genelist}_monoreporter_genes'))
      so_f <- readRDS(fns[['filtered']])
      so_f_ig_mono <- process_so(so_f[geneset, ], redo_SCT = F)
      saveRDS(so_f_ig_mono, fns[['filtered_ig_mono']])
    }
  }

  if (prep_cleaned) {
    if (so_rc(fns[['filtered_cleaned']]) || redo) {
      so_f <- readRDS(fns[['filtered']])
      so_f <- identify_outlying_clusters(so_f)

      # so_f@meta.data %>%
      #   dplyr::select(matches('dbscan')) %>% head
      # so_f %>% dplyr::select(matches('dbscan')) %>% head
      # tab <- so_f %>% dplyr::pull(dbscan_cluster) %>% table
      # tab / sum(tab)

      ## Only save objects if at least one outlying cluster was
      ## detected
      if (mean(so_f[['outlying_dbscan_cluster']] > 0)) {
        outlying_dbscan_clusters <- so_f@meta.data %>%
          dplyr::filter(outlying_dbscan_cluster == T) %>%
          pull(dbscan_cluster) %>%
          unique

        contrast_cluster_ids <-
          setdiff(
            unique(unlist(so_f[['dbscan_cluster']])),
            outlying_dbscan_clusters
          )

        cell_count <- table(so_f@meta.data$dbscan_cluster)
        ## FindMarkers fails without at least 3 cells
        test_clusters <- names(cell_count)[cell_count > 3] %>%
          intersect(outlying_dbscan_clusters)

        if (so_rc(fns[['marker']]) || redo) {
          markers <- map(auto_name(test_clusters),
            function(clust_id) {
            markers <-
              investigate_cluster(so_f,
                clust_id = clust_id,
                contrast_cluster_ids = contrast_cluster_ids,
                cluster_indicator = 'dbscan_cluster'
              )
            })
          saveRDS(markers, fns[['marker']])
        } else {
          markers <- readRDS(fns[['marker']])
        }

        dir.create(file.path(img_dir, glue::glue('exp{experiment}')),
          showWarnings = F)
        for (clust_id in names(markers)) {
          markers_to_gsea(markers = markers[[clust_id]],
            o_fn =
              file.path(img_dir,
                glue::glue('exp{experiment}'),
                glue::glue('outlying_clust_gsea-\\
                exp{experiment}-\\
                cluster={clust_id}-geneset=cellmarker.png')),
            pathways = cellmarker_genesets)
          markers_to_gsea(markers = markers[[clust_id]],
            o_fn =
              file.path(img_dir,
                glue::glue('exp{experiment}'),
                glue::glue('outlying_clust_gsea-\\
                exp{experiment}-\\
                cluster={clust_id}-geneset=reactome.png')),
            pathways = reactome_genesets)
        }

        ## Which hashtags are 'tainted'
        # so_f@meta.data %>%
        #   # dplyr::filter(seurat_clusters %in%
        #   # outlying_dbscan_clusters) %>%
        #   dplyr::mutate(
        #     outlier = seurat_clusters %in% outlying_dbscan_clusters) %>%
        #   dplyr::group_by(dominant_hashtag, outlier) %>%
        #   dplyr::summarize(n = n(), .groups = 'drop')

        # so_f_c <- so_f[,
        #   so_f@meta.data$seurat_clusters %nin%
        #   outlying_dbscan_clusters]
        so_f_c <- so_f[,
          so_f@meta.data$dbscan_cluster %nin% outlying_dbscan_clusters]
        so_f_c <- process_so(so_f_c, redo_SCT = T)
        saveRDS(so_f_c, fns[['filtered_cleaned']])

        if (so_rc(fns[['filtered_ig_cleaned']]) || redo) {
          geneset <- read_geneset(glue::glue('{genelist}_genes'))
          so_f_c <- readRDS(fns[['filtered_cleaned']])
          so_f_ig_cleaned <- process_so(so_f_c[geneset, ],
            redo_SCT = F)
          saveRDS(so_f_ig_cleaned, fns[['filtered_ig_cleaned']])
        }

        if (so_rc(fns[['filtered_ig_mono_cleaned']]) || redo) {
          geneset <- read_geneset(glue::glue('{genelist}_monoreporter_genes'))
          so_f_c <- readRDS(fns[['filtered_cleaned']])
          so_f_ig_mono_cleaned <- process_so(so_f_c[geneset, ],
            redo_SCT = F)
          saveRDS(so_f_ig_mono_cleaned, fns[['filtered_ig_mono_cleaned']])
        }
      }
    }
  }

  if (prep_mito) {
    if (so_rc(fns[['filtered_mito_reg']]) || redo) {
      so <- readRDS(fns[['filtered']])
      so <- process_so(so, redo_SCT = T, sct = 'mito')
      saveRDS(so, fns[['filtered_mito_reg']])
    }

    if (so_rc(fns[['filtered_ig_mito_reg']]) || redo) {
      geneset <- read_geneset(glue::glue('{genelist}_genes'))
      so <- readRDS(fns[['filtered_mito_reg']])
      so <- process_so(so[geneset, ], redo_SCT = F)
      saveRDS(so, fns[['filtered_ig_mito_reg']])
    }

    if (so_rc(fns[['filtered_ig_mono_mito_reg']]) || redo) {
      geneset <- read_geneset(glue::glue('{genelist}_monoreporter_genes'))
      so <- readRDS(fns[['filtered_mito_reg']])
      so <- process_so(so[geneset, ], redo_SCT = F)
      saveRDS(so, fns[['filtered_ig_mono_mito_reg']])
    }
  }

  if (prep_mito && prep_cleaned) {
    ## TODO
  }
}


#' Universal filtering rules for all single cell experiments except
#' for 5906 (combi coded hash tags change things a little)
#'
filter_HQ_cells <- function(so, hashtag = 'raw',
  max_percent_mt = 100, min_UMI = 0,
  min_fd_cc = 1, max_hashtag_evenness = 1) {
  library(Seurat)

  if (is.null(max_percent_mt) || is.na(max_percent_mt))
    max_percent_mt <- 100
  if (is.null(min_UMI) || is.na(min_UMI))
    min_UMI <- 0

  ## Infer hash tagging type for this experiment using sample
  ## annotation
  # experiment <- so@meta.data$exp[1]
  # sa <- get(glue::glue('sc_{experiment}_sample_annotation'))
  sa <- so@meta.data
  csa <- colnames(sa)
  if (all(c('HT1', 'HT2') %in% csa)) {
    so_f <- subset(so,
      bm_exp == T &
      fd_cc >= min_fd_cc &
      hashtag_evenness <= max_hashtag_evenness)
    so <- so_f
  } else if (all(c('HT1', 'HT2', 'hashtag') %in% csa == c(0, 0, 1))) {
    if (hashtag == 'HTO_classification') {
      if ('HTO_classification.global' %in% colnames(so@meta.data)) {
        so <- subset(so, HTO_classification.global == 'Singlet')
      } else if ('HTO_classification' %in% colnames(so@meta.data)) {
        so <- subset(so, HTO_classification == 'Singlet')
      }
    } else if (hashtag == 'dominant_hashtag') {
      so <- subset(so, fd_cc >= min_fd_cc &
        hashtag_evenness <= max_hashtag_evenness)
    }
  } else if (all(c('HT1', 'HT2', 'dominant_hashtag') %in%
      csa == c(0, 0, 1))) {
    so <- subset(so, fd_cc >= min_fd_cc &
      hashtag_evenness <= max_hashtag_evenness)
  } else {
    browser()
  }

  # so@meta.data %>%
  #   dplyr::select(percent.mt, nCount_RNA) %>%
  #   dplyr::summarize(across(everything(),
  #       list(min, median, max)))
  so <- subset(so,
    percent.mt <= max_percent_mt &
    N_UMI >= min_UMI)

  return(so)
}


#' Take a fresh Seurat object and run all the processing steps we
#' invariably want to be done
#'
#'
process_so <- function(so,
  redo_SCT = T, N_dims = 10,
  redo_geneset_scores = F, sct = 'vanilla',
  redo_variable_genes = T,
  clustering_resolution = 2) {

  if (!redo_variable_genes) {
    v_genes <- VariableFeatures(so)
    stopifnot(length(v_genes) > 0)
  }

  if (redo_geneset_scores || is.null(so@meta.data$percent.mt)) {
    so <- PercentageFeatureSet(so, pattern = '^MT-',
                               col.name = 'percent.mt')
    so <- PercentageFeatureSet(so, pattern = '^RP(S|L).*',
                               col.name = 'percent.ribo')
    so <- tryCatch({
      CellCycleScoring(so, s.features = cc.genes$s.genes,
                       g2m.features = cc.genes$g2m.genes,
                       set.ident = FALSE)
    }, error = function(e) { so })
  }

  if (redo_SCT) {
    sct <- switch(sct,
      'cc' = c('S.Score', 'G2M.Score'),
      'mito' = 'percent.mt',
      'vanilla' = NULL,
      sct
    )
    so <- SCTransform(so, vars.to.regress = sct)
  }

  if (redo_variable_genes) {
    so <- FindVariableFeatures(so)
  } else {
    VariableFeatures(so) <- v_genes
  }

  if (is.null(so@assays[['SCT']])) so <- ScaleData(so)
  so <- RunPCA(so, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:N_dims, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:N_dims, verbose = FALSE)
  so <- FindClusters(so, verbose = FALSE,
    resolution = clustering_resolution)
  so <- identify_outlying_clusters(so)
  so <- order_stim_group(so)
  return(so)
}


characterize_outlying_clusters <- function(so) {
  if (is.null(so[['outlying_dbscan_cluster']]))
    stop('DBSCAN has not been run yet for this Seurat object')
  if (!any(so[['outlying_dbscan_cluster']] > 0)) {
    return(NULL)
  }

  outlying_dbscan_clusters <- so@meta.data %>%
    dplyr::filter(outlying_dbscan_cluster == T) %>%
    pull(dbscan_cluster) %>%
    unique

  ## Cluster IDs to compare identified outlier(s) to
  contrast_cluster_ids <-
    setdiff(
      unique(unlist(so[['dbscan_cluster']])),
      outlying_dbscan_clusters
    )

  cell_count <- table(so@meta.data$dbscan_cluster)
  ## FindMarkers fails without at least 3 cells
  test_clusters <- names(cell_count)[cell_count > 3] %>%
    intersect(outlying_dbscan_clusters)

  markers <- purrr::map(maartenutils::auto_name(test_clusters),
    function(clust_id) {
      investigate_cluster(so,
        clust_id = clust_id,
        contrast_cluster_ids = contrast_cluster_ids,
        cluster_indicator = 'dbscan_cluster'
      )
    })
  return(markers)
}


add_6493_strep_counts <- function(so) {
  temp_en <- new.env()
  load(file.path(data_dir, 'raw_exp_6493', 'GEX_FB_6493.Rda'), 
    envir = temp_en)
  so_order <- rownames(so@meta.data) %>%
    stringr::str_replace('6493_', '')
  so <- AddMetaData(so,
    metadata = temp_en$counts[['ADT']]@counts[1, ][so_order],
    col.name = 'strep_counts'
  )
  so <- AddMetaData(so,
    metadata = temp_en$counts[['ADT']]@data[1, ][so_order],
    col.name = 'strep_counts_norm'
  )
  # summary(so@meta.data[['strep_counts']])
  # mean(so@meta.data[['strep_counts']] > 50)
  so@meta.data[['Ag']] <-
    ifelse(so@meta.data[['strep_counts']] > 50, 'Ag+', 'Ag-')
  return(so)
}


