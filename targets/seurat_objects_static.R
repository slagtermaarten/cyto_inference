sc_metas_targets <- list(
  ## Experiment 5310 (GCF experiment number)
  # Hashtag 1: Ag- cells from In vivo tumor + T cells(in muis)
  # Hashtag 2: Ag- cells from In vivo tumor, no T cells
  # Hashtag 3: in vitro cells, no stim
  # Hashtag 4: in vitro cells, 24h 100ng/ml IFNg stim
  # Hashtag 5: in vitro cells, 24h 10ng/ml TNFa stim
  # Hashtag 6: in vitro cells, 24h 100ng/ml IFNg stim + 10ng/ml TNFa stim
  # Hashtag 7: in vitro cells, 24h 100ng/ml IFNg stim + 10ng/ml TNFa
  # stim “mock” digested and sorted (to test effect digestion and
  # sorting on RNA profile).
  tar_target(sc_5310_sample_annotation, {
    bsa <- sample_annotation_exp5029
    data.frame(
      hash_tag = as.integer(1:7),
      sample_type =
        factor('sc',
          levels = levels(bsa$sample_type)),
      sample_origin = factor(
        c('in_vivo', 'in_vivo', 'in_vitro', 'in_vitro', 'in_vitro',
          'in_vitro', 'in_vitro'),
        levels = levels(bsa$sample_origin)),
      stim_group = factor(
        c('Exposed to T-cells in vivo', 'Unexposed in vivo',
          'Unstimulated in vitro', '100 ng/ml IFNy',
          '10 ng/ml TNFa',
          '100 ng/ml IFNy 10 ng/ml TNFa',
          '100 ng/ml IFNy 10 ng/ml TNFa')),
      duration = factor(c(NA, NA, 24, 24, 24, 24, 24),
        levels = levels(bsa$duration)),
      tau_rank = factor(c(NA, NA, 4, 4, 4, 4, 4),
        levels = c(NA, 1, 2, 3, 4)),
      sn_dilution = factor(c(NA, NA, NA, NA, NA, NA, NA),
        levels = levels(bsa$sn_dilution)),
      tnf_conc = factor(c(NA, NA, 0, 0, 10, 10, 10),
        levels = levels(bsa$tnf_conc)),
      tnf_rank = factor(c(NA, NA, 1, 1, 4, 4, 4),
        levels = c(NA, 1, 2, 3, 4)),
      ifn_conc = factor(c(NA, NA, 0, 100, 0, 100, 100),
        levels = levels(bsa$ifn_conc)),
      ifn_rank = factor(c(NA, NA, 1, 4, 1, 4, 4),
        levels = c(NA, 1, 2, 3, 4)),
      sc_digestion = c(T, T, F, F, F, F, T)
    ) %>% add_sc_condition_name()
  }),

  # sc_5906_sample_annotation <-
  #   matrix(c(1, 2, 'IFNg 2h 0 ng/ml',
  #            1, 3, 'IFNg 2h 0.01 ng/ml',
  #            1, 4, 'IFNg 2h 1 ng/ml',
  #            1, 6, 'IFNg 2h 100 ng/ml',
  #            2, 4, 'IFNg 6h 0 ng/ml',
  #            2, 6, 'IFNg 6h 0.01 ng/ml',
  #            2, 8, 'IFNg 6h 1 ng/ml',
  #            3, 4, 'IFNg 6h 100 ng/ml',
  #            3, 5, 'IFNg 12h 0 ng/ml',
  #            3, 7, 'IFNg 12h 0.01 ng/ml',
  #            4, 7, 'IFNg 12h 1 ng/ml',
  #            5, 6, 'IFNg 12h 100 ng/ml',
  #            5, 7, 'IFNg 24h 0 ng/ml',
  #            5, 8, 'IFNg 24h 0.01 ng/ml',
  #            6, 8, 'IFNg 24h 1 ng/ml',
  #            7, 8, 'IFNg 24h 100 ng/ml'), ncol = 3, byrow = T) %>%
  #   as.data.frame(stringsAsFactors = F) %>%
  #   magrittr::set_colnames(c('HT1', 'HT2', 'sample_group')) %>%
  #   dplyr::mutate(sample_name = copy(sample_group)) %>%
  #   dplyr::mutate(
  #     sample_name = gsub('([^ ]*) (.*)', '\\2 \\1', sample_name)) %>%
  #   dplyr::mutate(
  #     sample_name = maartenutils::variabilize_character(sample_name)) %>%
  #   dplyr::mutate(
  #     sample_name = gsub('ifng', 'ifn', sample_name)) %>%
  #   dplyr::mutate(sample_name = gsub('/', '_', sample_name)) %>%
  #   dplyr::mutate(
  #     sample_name = gsub('(\\d)_(\\d)', '\\1.\\2', sample_name)) %>%
  #   dplyr::mutate(
  #     duration = factor(as.numeric(gsub('.* (\\d+)h.*', '\\1', sample_group)))) %>%
  #   dplyr::mutate(
  #     tnf_conc = factor('0', levels =
#     levels(bsa$tnf_conc))) %>%
  #   dplyr::mutate(
  #     ifn_conc = factor(as.numeric(gsub('.*h ([\\.0-9]+) ng/ml$',
  #           '\\1', sample_group)))) %>%
  #   dplyr::mutate(
  #     sample_type = factor('sc',
  #       levels = levels(bsa$sample_type))) %>%
  #   dplyr::mutate(
  #     sample_origin = factor(c('in_vitro'),
  #       levels = levels(bsa$sample_origin))) %>%
  #   dplyr::mutate(stim_group = sprintf('%s ng/ml IFNy', ifn_conc)) %>%
  #   dplyr::mutate(
  #     sn_dilution = factor(0,
  #       levels = levels(bsa$sn_dilution))) %>%
  #   dplyr::mutate(sc_digestion = T) %>%
  #   dplyr::mutate(stim_group = ifelse(grepl('^0 ng/ml', stim_group),
  #                              'Unstimulated in vitro', stim_group))

  # tar_target(sc_6071_sample_annotation, {
  #   bsa <- sample_annotation_exp5029
  #   data.frame(
  #   hash_tag = 1:10,
  #   sample_type = factor('sc', levels =
  #     levels(bsa$sample_type)),
  #   sample_origin = factor(rep('in_vivo', 10),
  #     levels = levels(bsa$sample_origin)),
  #   stim_group = c('Unexposed in vivo',
  #                  'Unexposed in vivo',
  #                  'In vivo 100 ng/ml IFNy',
  #                  'In vivo 100 ng/ml IFNy',
  #                  'In vivo 10 ng/ml TNFa',
  #                  'In vivo 10 ng/ml TNFa',
  #                  'In vivo 100 ng/ml IFNy 10 ng/ml TNFa',
  #                  'In vivo 100 ng/ml IFNy 10 ng/ml TNFa',
  #                  'In vivo 100 ng/ml IFNy',
  #                  'In vivo 10 ng/ml TNFa'),
  #   duration = factor(rep(NA, 10),
  #     levels = levels(bsa$duration)),
  #   sc_digestion = c(rep(T, 10)),
  #   tumor_type = 'OVCAR5',
  #   aB2M_Ab = c(rep(F, 8), rep(T, 2))
  #   )}),

  # tar_target(sc_6072_sample_annotation, {
  #   bsa <- sample_annotation_exp5029
  #   dtf <- sc_6071_sample_annotation[1:8,]
  #   dtf$tumor_type <- 'NMM'
  #   dtf$stim_group <-
  #     c(rep('Unexposed in vivo', 2),
  #       rep('Exposed to T-cells in vivo', 6))
  #   return(dtf)
  #   }),

  tar_target(sc_6369_sample_annotation, {
    bsa <- sample_annotation_exp5029
    dtf <- tibble::tribble(~HT1, ~HT2, ~sample_group,
        1, 2, 'Unstimulated 2h',
        1, 3, 'IFNg 2h 1 ng/ml',
        1, 4, 'IFNg 2h 10 ng/ml',
        1, 6, 'IFNg 2h 100 ng/ml',
        2, 4, 'Unstimulated 6h',
        2, 6, 'IFNg 6h 1 ng/ml',
        2, 8, 'IFNg 6h 10 ng/ml',
        3, 4, 'IFNg 6h 100 ng/ml',
        3, 5, 'Unstimulated 12h',
        3, 7, 'IFNg 12h 1 ng/ml',
        4, 7, 'IFNg 12h 10 ng/ml',
        5, 6, 'IFNg 12h 100 ng/ml',
        5, 7, 'Unstimulated 24h',
        5, 8, 'IFNg 24h 1 ng/ml',
        6, 8, 'IFNg 24h 10 ng/ml',
        7, 8, 'IFNg 24h 100 ng/ml',
        1, 5, 'Unstimulated 24h frozen',
        2, 5, 'IFNg 24h 100 ng/ml frozen') %>%
    dplyr::mutate(across(c(HT1, HT2), as.integer)) %>%
    dplyr::mutate(sample_name = sample_group) %>%
    dplyr::mutate(
      sample_name = gsub('([^ ]*) (.*)', '\\2 \\1', sample_name)) %>%
    dplyr::mutate(
      sample_name = maartenutils::variabilize_character(sample_name)) %>%
    dplyr::mutate(
      sample_name = gsub('ifng', 'ifn', sample_name)) %>%
    dplyr::mutate(
      sample_name = gsub('/', '_', sample_name)) %>%
    dplyr::mutate(
      sample_name = gsub('(\\d)_(\\d)', '\\1.\\2', sample_name)) %>%
    dplyr::mutate(
      duration = factor(
        as.numeric(gsub('.* (\\d+)h.*', '\\1', sample_group)))) %>%
    dplyr::mutate(ifn_conc = suppressWarnings(as.numeric(
          gsub('.*h ([\\.0-9]+) ng/ml.*', '\\1', sample_group))
      )) %>%
    # debug_pipe %>%
    dplyr::mutate(ifn_conc = ifelse(is.na(ifn_conc), 0, ifn_conc)) %>%
    dplyr::mutate(ifn_conc = factor(ifn_conc)) %>%
    dplyr::mutate(sample_type =
      factor('sc',
        levels = levels(bsa$sample_type))) %>%
    dplyr::mutate(sample_origin =
      factor(c('in_vitro'),
        levels = levels(bsa$sample_origin))) %>%
    dplyr::mutate(stim_group = sprintf('%s ng/ml IFNy', ifn_conc)) %>%
    # dplyr::mutate(
    #   sn_dilution = factor(0, levels = levels(bsa$sn_dilution))) %>%
    dplyr::mutate(sc_digestion = F) %>%
    # debug_pipe %>%
    dplyr::mutate(stim_group =
      ifelse(grepl('^0 ng/ml', stim_group),
        'Unstimulated in vitro', stim_group)) %>%
    dplyr::mutate(frozen = grepl('frozen', sample_group))
    dtf[18, 'ifn_conc'] <- '100'
    return(add_sc_condition_name(dtf))
  }),

  tar_target(sc_6489_sample_annotation, {
    ## Bulk sample annotation, used for referencing some variables
    bsa <- sample_annotation_exp5029
    out <- tibble::tribble(~HT1, ~HT2, ~stim_group,
        10, 1, 'Unstimulated in vitro',
        10, 2, 'SN',
        10, 3, '1/5 SN',
        10, 4, '1/25 SN',
        10, 5, '1/125 SN',
        10, 6, '1/625 SN',
        10, 7, 'Unstimulated in vitro',
        10, 8, 'SN',
        10, 9, '1/5 SN',
        8, 1, '1/25 SN',
        8, 2, '1/125 SN',
        8, 3, '1/625 SN',
        8, 4, 'Unstimulated in vitro',
        8, 5, 'SN',
        8, 6, '1/5 SN',
        8, 7, '1/25 SN',
        8, 9, '1/125 SN',
        9, 1, '1/625 SN',
        9, 2, 'Unstimulated in vitro',
        9, 3, 'SN',
        9, 4, '1/5 SN',
        9, 5, '1/25 SN',
        9, 6, '1/125 SN',
        9, 7, '1/625 SN') %>%
    dplyr::mutate(across(c(HT1, HT2), as.integer)) %>%
    dplyr::mutate(duration = rep(c(2, 6, 12, 24), each = 6)) %>%
    dplyr::mutate(sample_origin = factor(c('in_vitro'),
        levels = levels(bsa$sample_origin))) %>%
    dplyr::mutate(sample_type = factor(c('sc'),
        levels = levels(bsa$sample_type))) %>%
    dplyr::mutate(
      sn_dilution = factor(
        rep(c(0, 1, 1/5, 1/25, 1/125, 1/625), length.out = 24))) %>%
    dplyr::mutate(sc_digestion = T) %>%
    add_sc_condition_name()
    return(out)
  }),

  tar_target(sc_6493_sample_annotation, {
    bsa <- sample_annotation_exp5029
    data.frame(
      hash_tag = as.integer(1:10),
      sample_type = factor('sc',
        levels = levels(bsa$sample_type)),
      sample_origin = factor('in_vivo',
        levels = levels(bsa$sample_origin)),
      stim_group = factor(c(
          rep('Exposed to T-cells in vivo', 3),
          rep('Unexposed in vivo', 2),
          rep('Exposed to T-cells in vivo', 3),
          rep('Unexposed in vivo', 2))),
      duration = factor(c(rep(16, 5), rep(44, 5))),
      sc_digestion = T) %>%
    dplyr::mutate(mouse = paste('6493', 1:10, sep = '_')) %>%
    add_sc_condition_name() %>%
    order_condition_name()
  }),

  # exp 6600:
  #       stim   tumor  # combo
  # 1    medium 1      2/5
  # 2    TNFa   1      2/6
  # 3    combi  1      2/7
  # 4    medium 2      2/8
  # 5    TNFa   2      2/10
  # 6    combi  2      4/6
  # 7    medium 3      4/7
  # 8    TNFa   3      4/8
  # 9    combi  3      4/9
  # 10   medium 4      4/10
  # 11   TNFa   4      3/5
  # 12   combi  4      3/6
  # 13   medium 5      3/7
  # 14   TNFa   5      3/8
  # 15   combi  5      3/9
  tar_target(sc_6600_sample_annotation, {
    bsa <- sample_annotation_exp5029
    tribble(~HT1, ~HT2, ~stim_group,
        2, 5, 'Unstimulated in vitro',
        2, 6, '10 ng/ml TNFa',
        2, 7, '100 ng/ml IFNy 10 ng/ml TNFa',
        2, 8, 'Unstimulated in vitro',
        2, 10, '10 ng/ml TNFa',
        4, 6, '100 ng/ml IFNy 10 ng/ml TNFa',
        4, 7, 'Unstimulated in vitro',
        4, 8, '10 ng/ml TNFa',
        4, 9, '100 ng/ml IFNy 10 ng/ml TNFa',
        4, 10, 'Unstimulated in vitro',
        3, 5, '10 ng/ml TNFa',
        3, 6, '100 ng/ml IFNy 10 ng/ml TNFa',
        3, 7, 'Unstimulated in vitro',
        3, 8, '10 ng/ml TNFa',
        3, 9, '100 ng/ml IFNy 10 ng/ml TNFa') %>%
    dplyr::mutate(across(c(HT1, HT2), as.integer)) %>%
    dplyr::mutate(
      duration = factor(rep(6, n()),
        levels = levels(bsa$duration))) %>%
    dplyr::mutate(
      sample_origin = factor(rep('ex_vivo_in_vitro', n()),
        levels = unique(c(
          levels(bsa$sample_origin),
          'ex_vivo_in_vitro')))) %>%
    dplyr::mutate(sample_type = factor(rep('sc', n()),
        levels = levels(bsa$sample_type))) %>%
    dplyr::mutate(mouse = paste('6600', rep(1:5, each = 3),
      sep = '_')) %>%
    dplyr::mutate(ifn_conc = factor(rep(c(0, 0, 100), 5),
      levels = levels(bsa$ifn_conc))) %>%
    dplyr::mutate(tnf_conc = factor(rep(c(0, 10, 10), 5),
      levels = levels(bsa$ifn_conc))) %>%
    dplyr::mutate(sc_digestion = T) %>%
    add_sc_condition_name() %>%
    { . }
  }),

  # Exp 6601 # combo
  # 1    2/5   TNFa 6h
  # 2    2/6   TNFa 6h
  # 3    2/7   TNFa 6h
  # 4    2/8   TNFa + IFNg 6h
  # 5    2/9   TNFa + IFNg 6h
  # 6    4/6   HBSS 6h
  # 7    4/7   HBSS 6h
  # 8    4/8   IFNg 24h
  # 9    4/9   IFNg 24h
  # 10   3/9   HBSS 24h
  # 11   3/5   HBSS 24h
  # 12   3/6   IFNg 48h
  # 13   3/7   HBSS 48h
  tar_target(sc_6601_sample_annotation, {
    bsa <- sample_annotation_exp5029
    tribble(~HT1, ~HT2, ~stim_group,
        2, 5, '10 ng/ml TNFa',
        2, 6, '10 ng/ml TNFa',
        2, 7, '10 ng/ml TNFa',
        2, 8, '100 ng/ml IFNy 10 ng/ml TNFa',
        2, 9, '100 ng/ml IFNy 10 ng/ml TNFa',
        4, 6, 'Unexposed in vivo',
        4, 7, 'Unexposed in vivo',
        4, 8, '100 ng/ml IFNy',
        4, 9, '100 ng/ml IFNy',
        3, 9, 'Unexposed in vivo',
        3, 5, 'Unexposed in vivo',
        3, 6, '100 ng/ml IFNy',
        3, 7, 'Unexposed in vivo') %>%
    dplyr::mutate(across(c(HT1, HT2), as.integer)) %>%
    dplyr::mutate(mouse = paste('6601', 1:13, sep = '_')) %>%
    dplyr::mutate(
      duration = factor(c(rep(6, 7), rep(24, 4), rep(48, 2)),
        levels = unique(c(levels(bsa$duration), 48)))) %>%
    dplyr::mutate(
      sample_origin = factor(rep('intratumoral_injection', n()),
        levels = unique(c(
          levels(bsa$sample_origin),
          'intratumoral_injection')))) %>%
    dplyr::mutate(sample_type = factor(rep('sc', n()),
        levels = levels(bsa$sample_type))) %>%
    dplyr::mutate(tnf_conc = factor(c(10, 10, 10, 10, 10, 
          0, 0, 0, 0, 0, 0, 0, 0), levels = levels(bsa$tnf_conc))) %>%
    dplyr::mutate(ifn_conc = factor(c(0, 0, 0, 100, 100, 0, 0, 
          100, 100, 0, 0, 100, 0), 
        levels = levels(bsa$ifn_conc))) %>%
    dplyr::mutate(sc_digestion = T) %>%
    add_sc_condition_name()
  }),

  tar_target(
    name = sc_6743_sample_annotation, 
    command = {
      bsa <- sample_annotation_exp5029
      # source('~/MirjamHoekstra/R/init.R')
      # bsa <- tar_read(sample_annotation_exp5029)
      tribble(~HT1, ~HT2, ~stim_group,
          1, 2, 'Ag-GAS + T',
          1, 3, 'Ag-GAS + T',
          2, 4, 'Ag-GAS + T',
          2, 5, 'Ag-GAS + T',
          3, 6, 'Ag-GAS + T',
          3, 7, 'Ag-GAS + T',
          4, 8, 'Mix PBS',
          4, 9, 'Mix PBS',
          5, 3, 'Mix PBS',
          5, 7, 'Mix PBS',
          6, 9, 'Mix PBS',
          6, 10, 'Mix PBS',
          7, 6, 'Mix + T',
          7, 10, 'Mix + T',
          8, 2, 'Mix + T',
          8, 10, 'Mix + T',
          9, 1, 'Mix + T',
          9, 3, 'Mix + T') %>%
      dplyr::mutate(across(c(HT1, HT2), as.integer)) %>%
      dplyr::mutate(
        duration = factor(rep(44, n()),
          levels = levels(bsa$duration))) %>%
      dplyr::mutate(
        sample_origin = factor(rep('in_vivo', n()),
          levels = levels(bsa$sample_origin))) %>%
      dplyr::mutate(sample_type = factor(rep('sc', n()),
          levels = levels(bsa$sample_type))) %>%
      dplyr::mutate(mouse = paste('6743', 1:n(), sep = '_')) %>%
      dplyr::mutate(sc_digestion = T) %>%
      add_sc_condition_name() %>%
      { . }
    }
  ),

  NULL
)


gen_sc_grid <- function(experiments = sc_e) {
  tibble::tibble(
    experiment = experiments,
    sc_sample_annotation = rlang::syms(purrr::map_chr(experiment,
        ~glue::glue('sc_{.x}_sample_annotation'))),
    fc_so = rlang::syms(purrr::map_chr(experiment,
        ~glue::glue('filtered_cleaned_so_{.x}')))
  )
}

### ------------     Single-cell      ------------ ###
seurat_object_targets <- list(

  tarchetypes::tar_map(
    values = gen_sc_grid(c(sc_e, '6743')),
    names = experiment,

    tar_target(vanilla_so, {
      raw_dir <- file.path(data_root, 
        # glue::glue('raw_exp_{experiment}'))
        paste0('raw_exp_', experiment))
      so <- create_seurat_object(raw_dir)
      if (experiment == '6493') so <- add_6493_strep_counts(so)
      return(so)
    }),

    tar_target(vanilla_HTO_ann_so, {
      # if (vanilla_so@meta.data$exp[1] == '6600') browser()
      seurat_process_HTO(
        so = vanilla_so,
        sa = sc_sample_annotation, 
        ncores = 1
      )
    ## This target kept being retriggered by a dependency on a
    ## static function
    }, cue = tar_cue(mode = 'never')),
    # }),

    tar_target(HTO_QC_so, {
      filter_HQ_cells(
        so = vanilla_HTO_ann_so,
        hashtag = default_filtering_opts[['hashtag']],
        max_percent_mt = 100,
        min_UMI = 0,
        min_fd_cc = default_filtering_opts[['min_fd_cc']],
        max_hashtag_evenness =
          default_filtering_opts[['max_hashtag_evenness']]
      )
    }),

    tar_target(filtered_so, {
      filter_HQ_cells(
        so = HTO_QC_so,
        hashtag = default_filtering_opts[['hashtag']],
        max_percent_mt = default_filtering_opts[['max_percent_mt']],
        min_UMI = default_filtering_opts[['min_UMI']]
      ) %>%
      process_so(sct = 'vanilla') %>%
      identify_outlying_clusters()
    }),

    tar_target(filtered_cleaned_so, {
      so_f_c <- filtered_so[,
        filtered_so@meta.data$outlying_dbscan_cluster == F]
      process_so(so_f_c, redo_SCT = TRUE)
    }),

    tar_target(sc_pseudobulk, {
      pseudo_bulk(
        filtered_cleaned_so, 
        u_var = c('condition_name', 'mouse')
      )
    }),

    tar_target(sc_pseudobulk_high_fidelity, {
      so <- filtered_cleaned_so
      so <- so[, so@meta.data$fd_cc > 4]
      pseudo_bulk(
        so,
        u_var = c('condition_name', 'mouse')
      )
    })
  ),

  tarchetypes::tar_map(
    values = gen_sc_grid(c(sc_e)),
    names = experiment,
    
    tar_target(filtered_cleaned_informative_so, {
      # sel_genes <- tar_read(gs_data_step) %>%
      #   dplyr::filter(!geneset %in%
      #     c('blacklist', 'none', 'lowly expressed'))
      so <- subset_feats(fc_so, read_geneset('informativeV15'))
      VariableFeatures(so) <- read_geneset('informativeV15')
      so <- process_so(so, redo_variable_genes = FALSE)
      return(so)
    }),

    tar_target(filtered_cleaned_mono_informative_so, {
      selection_crit_table_u <- tar_read(gs_data_step)
      # table(selection_crit_table_u$geneset)
      gs_data <- tar_read(gs_data_step)
      sel_genes <- gs_data %>%
        dplyr::filter(geneset %in% c('mono synergy', 'mono reporter'))
      so <- subset_feats(fc_so, sel_genes$gene)
      VariableFeatures(so) <- sel_genes$gene
      so <- process_so(so, redo_variable_genes = FALSE)
      return(so)
    })

    # tarchetypes::tar_map(
    #   values = tidyr::expand_grid(min_avg_exp = c(5, 25, 50)),

    #   tar_target(filtered_cleaned_magic_so, {
    #     magic_so <- add_MAGIC(
    #       so = filtered_cleaned_so,
    #       genelist = 'informativeV15',
    #       min_avg_exp = min_avg_exp
    #     )
    #     # magic_so@assays$SCT@counts <-
    #     #   magic_so@assays$MAGIC@counts
    #     # DefaultAssay(magic_so) <- 'SCT'
    #     return(magic_so)
    #   }),

    #   NULL
    # )
  ),

  NULL
)
