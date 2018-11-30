bulk_targets <- list(
  tar_target(sample_annotation_exp4910, {
    fastq_dir <- file.path(data_root, 'raw_exp_4910', 'fastq')
    fns <- list.files(fastq_dir)
    sn <- gsub('\\d{4}_\\d{1,2}_', '', fns) %>%
      { gsub('_[ATCG]{6}.*.fastq.gz', '', .) } %>%
      { gsub('OVCAR5_Ag-_', '', .) } %>%
      { gsub('IFNg', 'IFNy', .) } %>%
      { gsub('ifng', 'IFNy', .) } %>%
      { gsub('-_', '', .) } %>%
      { . }

    out <- tibble(
      # fastq_fn = fns,
      sample_name = sn,
      'duration' = c(48, 96, 48, 48, 48, 48, 96, 96, 96, 96),
      'stim_group' = c('Unstimulated in vitro', 'SN', 
        '100 ng/ml IFNy', '10 ng/ml TNFa', 'SN', 'SN', 
        'Unstimulated in vitro', '100 ng/ml IFNy', '10 ng/ml TNFa',
        'SN'),
      'IFNy_receptor' = c('+', '-', '+', '+', '+', '-', 
        '+', '+', '+', '+')
    ) %>%
    dplyr::mutate(sample_name = paste0('4910_', sample_name)) %>%
    dplyr::mutate(tnf_conc = ifelse(grepl('TNFa', sample_name), 
        10, 0)) %>%
    dplyr::mutate(ifn_conc = ifelse(grepl('IFNy', sample_name), 
        100, 0)) %>%
    dplyr::mutate(sn_dilution = ifelse(grepl('SN', sample_name), 
        1, 0)) %>%
    add_condition_name() %>%
    dplyr::mutate(ifnr_part = 
      ifelse(IFNy_receptor == '-', glue::glue(' - IFNyR-'), '')) %>%
    dplyr::mutate(condition_name =
      as.character(glue::glue('{condition_name}{ifnr_part}'))) %>%
    dplyr::select(-ifnr_part) %>%
    dplyr::mutate(exp = '4910') %>%
    dplyr::mutate(experiment = '4910') %>%
    dplyr::mutate(group_id = 'bulk-in_vitro')

    return(out)
  }),

  tar_target(kallisto_4910, {
    source(file.path(r_dir, 'bulk_kallisto.R'))
    fn <- file.path(data_root, 'raw_exp_4910', 
      'kallisto_tables', 'counts.tsv')

    out <- read_exp_matrix(fn = fn, mode = 'CPM') %>%
      { magrittr::set_colnames(., 
        cleanup_bulk_sample_names(colnames(.))) } %>%
      { magrittr::set_colnames(., 
        paste('4910_', colnames(.), sep = '')) }
    colnames(out) <- sample_annotation_exp4910$sample_name
    stopifnot(all(!is.na(colnames(out))))
    return(out)
  }),

  tar_target(bulk_4910_so, {
    exp4910 <- CreateSeuratObject(
      counts = kallisto_4910 %>%
        set_colnames(sample_annotation_exp4910$sample_name),
      meta.data = sample_annotation_exp4910 %>%
        as.data.frame() %>%
        set_rownames(NULL) %>%
        tibble::column_to_rownames('sample_name')
    )
    exp4910 <- SCTransform(exp4910)
    exp4910 <- RunPCA(exp4910, assay = 'SCT', npcs = 4, verbose = F)
    return(exp4910)
  }, iteration = 'list'),

  tar_target(kallisto_5029, {
    source(file.path(r_dir, 'bulk_kallisto.R'))
    fn <- file.path(data_root, 'raw_exp_5029', 
      'kallisto_tables', 'counts.tsv')

    out <- read_exp_matrix(fn = fn, mode = 'CPM') %>%
      { magrittr::set_colnames(., 
        cleanup_bulk_sample_names(colnames(.))) } %>%
      { magrittr::set_colnames(., 
        paste('5029_', colnames(.), sep = '')) }
    stopifnot(all(!is.na(colnames(out))))
    return(out)
  }),

  tar_target(sample_annotation_exp5029, {
    fastq_dir <-
      '/shared/gcf/m.slagter/shared/m.hoekstra/request_5029/fastq'
    fns <- list.files(fastq_dir)
    sn <- gsub('\\d{4}_\\d{1,2}_', '', fns) %>%
      { gsub('_[ATCG]{6}.*.fastq.gz', '', .) } %>%
      { . }
    out <- extract_bulk_meta_data_from_sample_names(sn) %>%
      dplyr::mutate(
        sample_name = paste('5029_', sample_name, sep = '')) %>%
      dplyr::mutate(exp = '5029') %>%
      dplyr::mutate(experiment = '5029') %>%
      dplyr::mutate(group_id = 'bulk-in_vitro') %>%
      add_condition_name()
    return(out)
  }),

  tar_target(bulk_5029_so, {
    stopifnot(
      sample_annotation_exp5029$sample_name ==
        colnames(kallisto_5029))
    exp5029 <- CreateSeuratObject(
      counts = kallisto_5029,
      meta.data = sample_annotation_exp5029 %>%
        as.data.frame() %>%
        set_rownames(NULL) %>%
        tibble::column_to_rownames('sample_name')
    )
    exp5029 <- SCTransform(exp5029)
    exp5029 <- RunPCA(exp5029, assay = 'SCT',
      npcs = 20, verbose = F)
    return(exp5029)
  }, iteration = 'list'),

  tar_target(sample_annotation_exp6434, {
    source(file.path(r_dir, 'experiment_meta_data.R'))
    fastq_dir <-
      '/shared/gcf/m.slagter/shared/m.hoekstra/6434/fastq_files'
    fns <- list.files(fastq_dir)
    sn <- gsub('\\d{4}_\\d{1,2}_', '', fns) %>%
      { gsub('_[ATCG]{10}-[ATCG]{10}_S\\d{1,2}_R[12]_001.fastq.gz',
        '', .) } %>%
      { . }
    extract_bulk_meta_data_from_sample_names(sn) %>%
      dplyr::mutate(
        sample_name = paste('6434_', sample_name, sep = '')) %>%
      dplyr::mutate(exp = '6434') %>%
      dplyr::mutate(experiment = '6434') %>%
      { . }
  }),

  # targets::tar_make(callr_function = NULL, names = 'kallisto_6434')
  tar_target(kallisto_6434, {
    source(file.path(r_dir, 'bulk_kallisto.R'))
    fn <- file.path(data_root, 'raw_exp_6434',
      'kallisto_tables', 'counts.tsv')
    out <- read_exp_matrix(fn = fn, mode = 'CPM')
    stopifnot(stringr::str_length(colnames(out)) < 1e2)
    stopifnot(stringr::str_length(rownames(out)) < 1e2)

    colnames(out) <- cleanup_bulk_sample_names(colnames(out))
    colnames(out) <- paste('6434_', colnames(out), sep = '')

    stopifnot(stringr::str_length(colnames(out)) < 1e2)
    stopifnot(stringr::str_length(rownames(out)) < 1e2)

    return(out)
  }),

  tar_target(bulk_6434_so, {
    library(magrittr)
    stopifnot(
      sample_annotation_exp6434$sample_name ==
        colnames(kallisto_6434))
    exp6434 <- CreateSeuratObject(
      counts = kallisto_6434,
      meta.data = sample_annotation_exp6434 %>%
        as.data.frame() %>%
        set_rownames(NULL) %>%
        tibble::column_to_rownames('sample_name'))
    # exp6434 <- Seurat::RenameCells(
    #   exp6434, glue('6434_{rownames(exp6434@meta.data)}'))
    # exp6434 <- Seurat::RenameCells(exp6434, '6434')
    exp6434@meta.data$group_id <- 'bulk-in_vitro'
    exp6434 <- SCTransform(exp6434)
    exp6434 <- RunPCA(exp6434, assay = 'SCT',
      npcs = 20, verbose = F)
    exp6434@meta.data %<>% add_condition_name()
    return(exp6434)
  }, iteration = 'list'),

  tar_target(kallisto_6623, {
    source(file.path(r_dir, 'bulk_kallisto.R'))
    fn <- file.path(data_root, 'raw_exp_6623', 
      'kallisto_tables', 'counts.tsv')
    out <- read_exp_matrix(fn = fn, mode = 'CPM') %>%
      { magrittr::set_colnames(., 
        cleanup_bulk_sample_names(colnames(.))) } %>%
      { magrittr::set_colnames(., 
        paste('6623_', colnames(.), sep = '')) }
    colnames(out) <- 
      stringr::str_replace(colnames(out), 'ifngy', 'ifny')
    stopifnot(all(!is.na(colnames(out))))
    return(out)
  }),

  tar_target(sample_annotation_exp6623, {
    # tar_read(sample_annotation_exp5029)
    fn <- file.path(data_root, 'raw_exp_6623', 
      'experimental_design.txt')
    dtf <- readr::read_tsv(fn)
    dtf$tnf_conc <- 
      as.numeric(stringr::str_replace(dtf$tnf_conc, ' ng/ml .*', ''))
    dtf$ifn_conc <- 
      as.numeric(stringr::str_replace(dtf$ifn_conc, ' ng/ml .*', ''))

    M <- tar_read(kallisto_6623)
    stopifnot(!any(duplicated(colnames(M))))
    ## Translate the meta sample names to the ones in the kallisto
    ## object
    cp <- recover_garbled_string(dtf$sample_name, colnames(M))
    stopifnot(!any(duplicated(cp)))
    ## Ordering of meta in colnames of M
    idxs <- match(colnames(M), unname(cp))
    stopifnot(unname(cp)[idxs] == colnames(M))
    cp <- cp[idxs]
    # cp %>%
    #   { dups = .[duplicated(.)]; .[. %in% dups] }
    ## Find the positions of the translated samples names in the
    ## kallisto object and reorder the meta
    dtf <- dtf[idxs, ]
    stopifnot(all(dtf$sample_name == names(cp)))
    dtf$sample_name <- unname(cp)

    dtf$tnf_part <- with(dtf, 
      ifelse(tnf_conc > 0, glue::glue('{tnf_conc} ng/ml TNFa'), ''))
    dtf$ifn_part <- with(dtf, 
      ifelse(ifn_conc > 0, glue::glue('{ifn_conc} ng/ml IFNy'), ''))
    dtf$stim_group <- with(dtf, 
      paste(ifn_part, tnf_part, sep = ' '))
    dtf$stim_group[dtf$stim_group == ' '] <- 'Unstimulated in vitro'
    duration_part <- with(dtf,
      ifelse(is.na(duration) | duration == 'Unknown', '', 
        glue::glue(' - {as.character(duration)}h')))
    blocker_part <- with(dtf,
        glue::glue(' - {as.character(recovery_duration)}h {as.character(blocker)}'))
    dtf$condition_name <-
      as.character(glue::glue('{dtf$stim_group}{duration_part}{blocker_part}'))
    dtf$condition_name <-
      stringr::str_replace_all(dtf$condition_name, '  ', ' ')
    dtf$condition_name <-
      stringr::str_replace(dtf$condition_name, '^ ', '')
      
    dtf %>%
      dplyr::select(sample_name, condition_name)

    dtf <- dtf
      dplyr::mutate(exp = '6623') %>%
      dplyr::mutate(experiment = '6623') %>%
      { . }

    return(dtf)
  }),

  tar_target(bulk_6623_so, {
    sa <- tar_read(sample_annotation_exp6623)
    M <- tar_read(kallisto_6623)

    stopifnot(all(sa$sample_name == colnames(kallisto_6623)))
    rownames(sa) <- sa$sample_name
    sa <- as.data.frame(sa)
    stopifnot(all(rownames(sa) == colnames(kallisto_6623)))
    exp6623 <- SeuratObject::CreateSeuratObject(
      counts = M, meta.data = sa)
    # exp6623@meta.data
    exp6623 <- SCTransform(exp6623)
    exp6623@meta.data$exp <- '6623'
    exp6623@meta.data$experiment <- '6623'
    exp6623@meta.data$group_id <- 'bulk-in_vitro'
    # exp6623 <- FindVariableFeatures(exp6623)
    exp6623 <- RunPCA(exp6623, 
      # assay = 'SCT',
      npcs = 20, verbose = F)
    return(exp6623)
  }, iteration = 'list'),

  tar_target(bulk_concat, {
    CM1 <- tar_read(kallisto_5029)
    CM2 <- tar_read(kallisto_6434)
    CM3 <- tar_read(kallisto_4910)
    CM4 <- tar_read(kallisto_6623)
    scn <- intersect(rownames(CM1), rownames(CM2)) %>%
      intersect(rownames(CM3)) %>%
      intersect(rownames(CM4))
    cbind(CM1[scn, ], CM2[scn, ], CM3[scn, ], CM4[scn, ])
  }),

  tar_target(bulk_concat_corr, {
    # credentials::set_github_pat()
    # devtools::install_github('zhangyuqing/sva-devel')
    sa1 <- tar_read(sample_annotation_exp5029)
    sa2 <- tar_read(sample_annotation_exp6434)
    sa3 <- tar_read(sample_annotation_exp4910)
    sa4 <- tar_read(sample_annotation_exp6623)
    # stopifnot(
    #   c(sa1$sample_name, sa2$sample_name, sa3$sample_name) == 
    #     colnames(bulk_concat))

    library(sva)
    batch <- c(
      rep('5029', nrow(sa1)), 
      rep('6434', nrow(sa2)),
      rep('4910', nrow(sa3)),
      rep('6623', nrow(sa4))) %>%
      as.factor %>% as.integer()
    stopifnot(ncol(bulk_concat) == length(batch))

    # process_sa <- function(sa) {
    #   dplyr::select(sa, duration, ifn_conc, tnf_conc, sn_dilution) %>%
    #     numerify_regressors()
    # }
    # covar_mat <- dplyr::bind_rows(
    #   process_sa(sa1), 
    #   process_sa(sa2),
    #   process_sa(sa3),
    #   process_sa(sa4)
    # )
    # out <- sva::ComBat_seq(
    #   bulk_concat, batch = batch, covar_mod = covar_mat)
    out <- sva::ComBat_seq(
      bulk_concat, batch = batch)

    return(out)
  }),

  tar_target(comb_bulk_so, {
    b1 <- tar_read(bulk_5029_so)
    b2 <- tar_read(bulk_6434_so)
    b3 <- tar_read(bulk_4910_so)
    # b3@meta.data
    b4 <- tar_read(bulk_6623_so)
    # head(b1@meta.data); head(b2@meta.data)
    sa <- merge(merge(merge(b1, b2), b3), b4)@meta.data
    stopifnot(all(!is.na(sa$tnf_conc)))
    M <- tar_read(bulk_concat_corr)
    # stopifnot(colnames(M) == rownames(sa))
    # colnames(M) <- rownames(sa)

    exp_bulk_comb <- CreateSeuratObject(
      counts = M, 
      meta.data = sa %>% set_rownames(colnames(M))
    )
    # exp_bulk_comb[which(is.na(exp_bulk_comb@meta.data$tnf_conc)), ]
    stopifnot(all(!is.na(exp_bulk_comb@meta.data$tnf_conc)))

    exp_bulk_comb <- SCTransform(exp_bulk_comb)
    exp_bulk_comb <- RunPCA(exp_bulk_comb, assay = 'SCT',
      npcs = 20, verbose = F)

    return(exp_bulk_comb)
  }, iteration = 'list'),

  tar_target(sample_annotation_bulk_comb_no_6623_no_sn, {
    sa1 <- tar_read(sample_annotation_exp5029) %>%
      { .[is.na(.$sn_dilution) | .$sn_dilution == '0', ] }
    sa2 <- tar_read(sample_annotation_exp6434)
    # sa3 <- tar_read(sample_annotation_exp4910) %>%
    #   { .[.$IFNy_receptor == '+' & .$sn_dilution == '0', ] }
    sa <- harmonize_bind_rows(sa1, sa2)
    return(sa)
  }),

  tar_target(covar_bulk_comb_no_6623_no_sn, {
    process_sa <- function(sa) {
      out <- dplyr::select(sa, 
        experiment, duration, ifn_conc, tnf_conc) %>%
        numerify_regressors() %>%
        dplyr::mutate(across(everything(), factor))
      rownames(out) <- sa$sample_name
      return(out)
    }
    process_sa(sample_annotation_bulk_comb_no_6623_no_sn)
  }),

  tar_target(bulk_concat_no_6623_no_sn, {
    CM1 <- tar_read(kallisto_5029)
    CM2 <- tar_read(kallisto_6434)
    # CM3 <- tar_read(kallisto_4910)
    srn <- rownames(CM1) %>%
      intersect(rownames(CM2)) %>%
      # intersect(rownames(CM3)) %>%
      { . }
    filter_samples <- function(M) {  
      allowed_sn <- intersect(
        colnames(M),
        sample_annotation_bulk_comb_no_6623_no_sn$sample_name
      )
      # stopifnot(length(allowed_sn) >= 50L)
      # stopifnot(length(allowed_sn) <= 500L)
      return(M[, allowed_sn])
    }
    cbind(
      filter_samples(CM1[srn, ])
      , filter_samples(CM2[srn, ])
      # , filter_samples(CM3[srn, ])
    )
  }),

  tar_target(bulk_concat_corr_no_6623_no_sn, {
    # pacman::p_load('genefilter')
    library(sva)
    bulk_concat_no_6623_no_sn <- tar_read(bulk_concat_no_6623_no_sn)
    sample_annotation_bulk_comb_no_6623_no_sn <-
      tar_read(sample_annotation_bulk_comb_no_6623_no_sn)
    covar_bulk_comb_no_6623_no_sn <- tar_read(covar_bulk_comb_no_6623_no_sn)

    stopifnot(
      colnames(bulk_concat_no_6623_no_sn) == 
      sample_annotation_bulk_comb_no_6623_no_sn$sample_name
    )
    stopifnot(
      colnames(bulk_concat_no_6623_no_sn) == 
      rownames(covar_bulk_comb_no_6623_no_sn)
    )

    out <- sva::ComBat_seq(
      bulk_concat_no_6623_no_sn, 
      batch = as.integer(covar_bulk_comb_no_6623_no_sn$experiment),
      covar_mod = 
        covar_bulk_comb_no_6623_no_sn %>%
        dplyr::select(-matches('exp')) %>%
        as.data.frame()
    )

    return(out)
  }),

  tar_target(comb_bulk_so_no_6623_no_sn, {
    sa <- as.data.frame(sample_annotation_bulk_comb_no_6623_no_sn)
    rownames(sa) <- sa$sample_name
    exp_bulk_comb <- CreateSeuratObject(
      counts = bulk_concat_corr_no_6623_no_sn, 
      meta.data = sa
    )
    # exp_bulk_comb[which(is.na(exp_bulk_comb@meta.data$tnf_conc)), ]
    stopifnot(all(!is.na(exp_bulk_comb@meta.data$tnf_conc)))
    exp_bulk_comb <- SCTransform(exp_bulk_comb)
    exp_bulk_comb <- RunPCA(exp_bulk_comb, assay = 'SCT',
      npcs = 20, verbose = F)
    return(exp_bulk_comb)
  }, iteration = 'list'),

  tar_target(comb_bulk_so_no_6623_no_sn_uncorrected, {
    sa <- as.data.frame(sample_annotation_bulk_comb_no_6623_no_sn)
    rownames(sa) <- sa$sample_name
    exp_bulk_comb <- CreateSeuratObject(
      counts = bulk_concat_no_6623_no_sn, 
      meta.data = sa
    )
    # exp_bulk_comb[which(is.na(exp_bulk_comb@meta.data$tnf_conc)), ]
    stopifnot(all(!is.na(exp_bulk_comb@meta.data$tnf_conc)))
    exp_bulk_comb <- SCTransform(exp_bulk_comb)
    exp_bulk_comb <- RunPCA(exp_bulk_comb, assay = 'SCT',
      npcs = 20, verbose = F)
    return(exp_bulk_comb)
  }, iteration = 'list'),

  tar_target(comb_bulk_so_no_6623_no_sn_dedup, {
    # so <- tar_read(comb_bulk_so_no_6623_no_sn_uncorrected)
    so <- comb_bulk_so_no_6623_no_sn

    dup_settings <- 
      so@meta.data[, c('tnf_conc', 'ifn_conc', 'duration')] %>%
      distinct()

    m_so <- pseudo_bulk(
      so, 
      u_var = c('tnf_conc', 'ifn_conc', 'duration'),
      datatype = 'data'
    )
    
    m_so@meta.data$experiment <- m_so@meta.data$experiment <- 
      gen_exp_string(so@meta.data$experiment)
    # nrow(m_so@meta.data)

    return(m_so)
  }, iteration = 'list'),

  ## NB different from 'comb_bulk_so'
  tar_target(bulk_comb_so, comb_bulk_so_no_6623_no_sn_dedup),

  tarchetypes::tar_map(
    names = experiment,

    values = tibble::tibble(
      experiment = bulk_e,
      so = rlang::syms(glue::glue('bulk_{experiment}_so'))
    ),

    tar_target(read_count_summary, {
      gen_so_read_count_summary(so)
    })
  ),

  NULL
)
