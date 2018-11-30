library(targets)
gene_selection_targets <- list(
  ### ------------    Gene selection    ------------ ###

  tar_target(gs_data_step_pre, load_gs_data(),
    iteration = 'vector', format = 'rds'),

  tar_target(gs_data_step, {
    update_selection_crit_table_gs(
      gs_data_step_pre,
      version_idx = 10L
    )
  }, iteration = 'vector', format = 'rds'),

  tar_target(default_gene_sets, {
    mono_reps <- gs_data_step %>%
      dplyr::filter(geneset == 'mono reporter')
     list(
      TNFa = mono_reps %>%
        dplyr::filter(TNFa_bias >= .5) %>%
        dplyr::pull(gene),
      IFNy = mono_reps %>%
        dplyr::filter(TNFa_bias < .5) %>%
        dplyr::pull(gene),
      synergy =
        gs_data_step %>%
        dplyr::filter(geneset %in% c('synergy', 'mono synergy')) %>%
        dplyr::filter(synergy_score >= 1.2) %>%
        dplyr::pull(gene),
      CXCL10 = 'CXCL10',
      HK_genes =
        gs_data_step %>%
        dplyr::filter(
          max_t < 1 &
          Amean >= 3 &
          timecor_median < .8
        ) %>%
        dplyr::pull(gene)
    )
  }),

  tar_target(mirjam_genes,
    file.path(data_dir,
      'proposed genes for time scores to Maarten.xlsx') %>%
      maartenutils::read_excel() %>%
      tidyr::pivot_longer(cols = everything(),
        values_to = 'gene', names_to = 'geneset') %>%
      dplyr::filter(!is.na(gene)) %>%
      dplyr::mutate(
        geneset = stringr::str_replace(geneset, 'IFNg', 'IFNy')) %>%
      dplyr::mutate(
        geneset = stringr::str_replace(
          geneset, 'tafelberg', 'plateau')) %>%
      dplyr::nest_by(geneset) %>%
      tibble::deframe() %>%
      purrr::map(~unname(unlist(.x)))),

  tar_target(time_informative_genesets, {
    tidyr::expand_grid(
      'duration' = c(2, 6, 12, 24),
      'gs' = names(default_gene_sets)[1:3]
    ) %>%
    dplyr::mutate(genes =
      purrr::pmap(., function(duration, gs) {
        gs_data_step %>%
          dplyr::filter(gene %in% default_gene_sets[[gs]]) %>%
          dplyr::filter(!is.na(max_timepoint)) %>%
          dplyr::filter(max_timepoint == .env[['duration']]) %>%
          dplyr::filter(time_class == 'dominant_tp') %>%
          dplyr::filter(gene %nin% c('IRAK2', 'PARP9', 'WDR25',
              'SEMA3C', 'ETS1', 'LCP1', 'ZC3H12C', 'BTG3', 'LCP1',
              'KIF21B', 'RND3')) %>%
          dplyr::pull(gene)
      })
    ) %>%
    {
       .[purrr::map_int(.$genes, length) > 0, ]
    } %>%
    { rlang::set_names(.$genes, paste0(.$gs, '-', .$duration)) } %>%
    { . }
  }),

  tar_target(clustered_mono_genes, {
    mono_genes <- unlist(tar_read(default_gene_sets)[1:3])
    bulk_so <- bulk_5029_so
    out <- 
      cluster_genes(
        genes = mono_genes,
        k = 4,
        so = bulk_so[, bulk_so@meta.data$sn_dilution == '0']
      ) %>%
      { .$genesets } %>%
      { setNames(., paste0('mono_', names(.))) }
    return(out)
  }),

  tar_target(
    name = clustered_IFNy_genes,
    command = {
      N_gene_clusters = 5
      MV = .1
      bulk_so <- bulk_5029_so
      all_genesets <- c(default_gene_sets, mirjam_genes, 
        time_informative_genesets)

      cg <- cluster_genes(
        so = bulk_so[, bulk_so@meta.data$sn_dilution == '0'], 
        min_var = MV,
        k = N_gene_clusters,
        genes = names(all_genesets) %>%
          stringr::str_subset('IFNy') %>%
          # setdiff('IFNy') %>%
          { all_genesets[.] } %>%
          unlist
      )

      cg$genesets %>%
        { setNames(., paste0('IFNy_', names(.))) }
    }
  ),

  tar_target(
    name = clustered_TNFa_genes,
    command = {
      N_gene_clusters = 4
      MV = .1
      bulk_so <- bulk_5029_so
      all_genesets <- c(default_gene_sets, mirjam_genes, 
        time_informative_genesets)

      cg <- cluster_genes(
        so = bulk_so[, bulk_so@meta.data$sn_dilution == '0'], 
        min_var = MV,
        k = N_gene_clusters,
        genes = names(all_genesets) %>%
          stringr::str_subset('TNFa') %>%
          { all_genesets[.] } %>%
          unlist
      )

      cg$genesets %>%
        { setNames(., paste0('TNFa_', names(.))) }
    }
  ),

  tar_target(
    name = clustered_limma_MR_genes,
    command = {
      human_MR_gene_group <- readRDS(file.path(rds_dir, 
          'human_MR_gene_group.rds'))
      clust2list(human_MR_gene_group) %>%
        { setNames(., paste0('limma_MR', names(.))) }
    }
  ),

  tar_target(
    name = human_HM_genesets,
    command = {
      library(genesets)
      hallmark_gs <- filter_gmt(
        pattern = 'HALLMARK', gmt_pattern = 'msigdb', 
        type = 'symbols')
      caplist <- caplist_def %>% c('via')
      names(hallmark_gs) <-
        names(hallmark_gs) %>%
        stringr::str_replace('HALLMARK_', '') %>%
        tolower() %>%
        # stringr::str_replace('_', ' ') %>%
        simple_cap(cap_first_word_only = T, caplist = caplist) %>%
        { paste0('HM ', .) } %>%
        { . }
      return(hallmark_gs)
    }
  ),

  tar_target(all_genesets, 
    command = {
      c(default_gene_sets, mirjam_genes, time_informative_genesets, 
        clustered_mono_genes, clustered_IFNy_genes,
        clustered_TNFa_genes, clustered_limma_MR_genes,
        human_HM_genesets)
    }
  ),

  tar_target(non_reproducible_genes, identify_non_reducible_genes(),
    iteration = 'vector'),

  tar_target(gene_reproducibility_score, {
    tar_read(non_reproducible_genes) %>%
      dplyr::filter(tnf_conc == 10 | ifn_conc == 100) %>%
      dplyr::group_by(gene) %>%
      # dplyr::summarize(score = sum(weight_vec[q9q1om_b], na.rm = T))
      dplyr::summarize(score = mean(q9q1om)) %>%
      dplyr::mutate(score_group = as.integer(cut(score,
          include.lowest = TRUE, ordered_result = TRUE,
          breaks = quantile(score, probs = c(0, .25, .5, .75, 1))))
      )
  })
)
