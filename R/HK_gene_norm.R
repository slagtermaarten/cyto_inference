#' Assess the robustness of using non-responsive genes as as estimate
#' of technical library size using different strategies 
#'
#'
HK_gene_norm <- function(M, 
  genes, 
  min_GDR = .75,
  outlier_norm = F,
  int_f = 'sum',
  N_folds = 2,
  N_repeats = 25L,
  strat_sampling = F,
  assay = 'SCT', 
  datatype = 'counts') {

  if (F) {
    # M_raw <- so2M(df, 'RNA', 'counts')
    # summary(colSums(M))
    # summary(colSums(M_raw))
    M <- M[names(which(compute_GDR(M) > .75)), ]
    cov_ <- apply(M, 1, function(x) sd(x) / mean(x))
    genes <- names(which(cov_ < .3))
    selection_crit_table_u %>%
      dplyr::filter(gene %in% genes) %>%
      pull(max_t)
    if (F) {
      exp_dynamics_panel(
        features = genes,
        max_plots_per_doc = 1000L,
        version = 'inferred_HK_in_vitro_5310',
        lookup_data = targets::tar_read('kallisto_5029'),
        redo = F,
        leave_out_sn = F)
    }
    M <- M[, ]
    # summary(cov_)
  } else {
    genes <- names(which(compute_GDR(M) > min_GDR))
    M <- M[genes, ]
  }

  if (F) {
    M_s <- t(apply(M, 1, min_max_scaling))
    tM_GE <- gen_HM(
      M = M_s,
      sa = NULL,
      N_genes = NULL,
      N_hl_genes = NULL,
      cluster_rows = T,
      cluster_columns = T,
      show_row_dend = T,
      show_column_dend = T,
      show_row_names = F,
      show_column_names = F,
      row_names_side = 'right',
      value_name = 'Gene expression'
    )
    pdf(file.path(img_dir, glue('HK_gene_HM.pdf')), height = 10)
    draw(HM_GE)
    dev.off()
  }

  if (outlier_norm) {
    M <- outlier_norm(M)
  }

  if (strat_sampling) {
    ## Stratified sampling, doesn't seem to increase agreement between
    ## differen subsamplings much and shouldn't
    ## come to think of it more deeply
    Amean_bin <- tibble(gene = rownames(M)) %>%
      left_join(selection_crit_table_u, by = 'gene') %>%
      dplyr::mutate(Amean_bin = 
        cut(Amean, breaks = seq(min(Amean), max(Amean), by = .5))
      ) %>%
      pull(Amean_bin) %>%
      as.character()
    dtf <- M2tibble(M) %>% mutate(Amean_bin = Amean_bin)
    M_rs <- rsample::vfold_cv(dtf, v = N_folds, repeats = N_repeats, 
      strata = 'Amean_bin')
  } else {
    M_rs <- rsample::vfold_cv(M, v = N_folds, repeats = N_repeats) 
  }

  compute_gs = switch(int_f[1],
    'median' = function(x) apply(x, 2, median),
    'sum' = function(x) apply(x, 2, sum)
  )

  r.squared <- map_dbl(M_rs$splits, function(.x) {
    if (strat_sampling) {
      M <- .x$data %>% dplyr::select(where(is.numeric))
    } else {
      M <- .x$data
    }
    df <- tibble(
      "in_data" = compute_gs(M[.x$in_id, ]),
      "out_data" = compute_gs(M[-.x$in_id, ])
    )
    mod <- lm(out_data ~ in_data, data = df)
    if (F) {
      p <- qplot(df$in_data, df$out_data, alpha = .2)
      print_plot(p,
        fn = file.path(img_dir, glue('HK_CV_scat.png')))
    }
    summary(mod)$r.squared
  })

  if (N_folds == 2) {
    ## Duplicate values when N_folds == 2
    r.squared <- r.squared[seq(2, N_folds * N_repeats, by = 2)]
  }

  return(
    list(
      'rs.10' = quantile(r.squared, .10),
      'rs.50' = quantile(r.squared, .50),
      'rs.90' = quantile(r.squared, .90),
      'N_genes' = nrow(M)
    )
  )
}


#' Titrate single cell filtering settings and record statistics on the
#' gene expression of selected cells. We aim to put thresholds at a
#' location where cell inclusion is i) maximally lenient but ii)
#' stringent enough to be stable (i.e. gene expression of housekeeping
#' genes not affected by further increasing the leniency of the
#' threshold).
#'
#'
titrate_filtering_settings <- result_cacher(
  f = function(experiment = '6493', filtering_opts) {
    HK_genes <- read_geneset('OVCAR5_HK_genes')

    source_fn <- compile_so_fns(
      experiment = experiment,
      filtering_opts = filtering_opts
    )$HTO_QC

    # object_name <- glue('exp{experiment}_HTO_QC_sc')
    # obj_fn <- file.path(rds_dir, glue('{object_name}.rds'))
    # if (exists(object_name, envir = .GlobalEnv)) {
    #   so <- get(object_name)
    # } else {
    #   so <- readRDS(obj_fn)
    #   assign(object_name, so, envir = .GlobalEnv)
    # }

    so <- readRDS(source_fn)

    settings <- tidyr::expand_grid(
      mito_thresh = seq(.1, 1, by = .1),
      ncount_thresh = c(500, 750, 1000, 1250, 1500, 1750, 2000)
    )

    # settings <- tidyr::expand_grid(
    #   mito_thresh = .1,
    #   ncount_thresh = c(1000)
    # )

    res <- settings %>%
      furrr::future_pmap_dfr(function(mito_thresh, ncount_thresh) {
        so_f <- so@meta.data %>%
          dplyr::filter(
            percent.mt <= mito_thresh & nCount_RNA >= ncount_thresh) %>%
          { so[, rownames(.)] }
        so_f <- tryCatch(SCTransform(so_f),
          error = function(e) { print(e); NULL })
        if (is.null(so_f)) return(NULL)
        out <- tidyr::expand_grid(
          gene = HK_genes,
          ht = unique(so_f@meta.data$dominant_hashtag)) %>%
          pmap(function(gene, ht) {
            if (gene %nin% rownames(so_f)) return(NULL)
            ## Extract gene expression for this gene in this hashtag
            ## group
            so_ff <- so_f@meta.data %>%
              dplyr::filter(dominant_hashtag == ht) %>%
              { so_f[, rownames(.)] }
            v <- as.numeric(so_ff@assays$SCT[gene, ])
            list(
              'gene' = gene,
              'ht' = ht,
              'mito_thresh' = mito_thresh,
              'ncount_thresh' = ncount_thresh,
              'N' = length(v),
              'N_nonzero' = sum(v > 0),
              'var' = var(v, na.rm = T),
              'med' = median(v, na.rm = T),
              'IQR' = IQR(v, na.rm = T))
          })
      })
    return(res)
  }, filename = function(experiment) {
    file.path(rds_dir, 'sc_stat_setting_titration',
      glue('exp{experiment}.rds'))
  }, min_mod_time = '2021-09-09 12:00'
)


