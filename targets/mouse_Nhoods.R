mouse_Nhoods_targets <- tarchetypes::tar_map(
  names = type,
  values = tibble(type = 'mouse'),
  unlist = T,

  tar_target(
    name = agg_neighbourhoods,
    command = {
      so <- tar_read(filtered_so_6743)
      limma_MRs <- tar_read(limma_MRs_mouse)

      agg_neighbourhoods <-
        perform_agg_neighbourhoods(
          so, genes = limma_MRs, genelist = NULL)

      if (TRUE) {
        # agg_neighbourhoods <- tar_read(agg_neighbourhoods_mouse)
        # agg_neighbourhoods <- runPCA(agg_neighbourhoods, ncomponents=30)
        # agg_neighbourhoods <- runUMAP(agg_neighbourhoods)
        print_plot_eval(
          print(plotUMAP(agg_neighbourhoods, 
              colour_by = 'stim_group')),
          width = 17.4, height = 10,
          filename = file.path(Sys.getenv('img_dir'),
            glue::glue('milo_UMAP_mouse.pdf')))
      }

      return(agg_neighbourhoods)
    },
    iteration = 'list'
  ),


  # tar_target(NH_test_DA,
  #   command = {
  #     agg_neighbourhoods <- tar_read(agg_neighbourhoods_mouse)

  #     if (interactive())
  #       source('~/MirjamHoekstra/R/init.R')
  #     NH_test_DA <- test_milo_DA_Nhoods(
  #       test_mode = 'stim_group_only',
  #       sce = agg_neighbourhoods
  #     )

  #     return(NH_test_DA)
  #   },
  #   iteration = 'list'
  # ),

  tar_target(
    name = milo_plots,
    command = {

      ## Plot fraction of each condition in each neighbourhood
      library(miloR)
      library(SummarizedExperiment)
      sums <- apply(nhoodCounts(agg_neighbourhoods), 2, sum) %>%
        { . / sum(.) }
      plots <-
        1:ncol(nhoodCounts(agg_neighbourhoods)) %>%
        purrr::map(function(ci) {
          idxs <- unlist(nhoodIndex(agg_neighbourhoods))
          colData(agg_neighbourhoods)$col_val <- NA_real_
          colData(agg_neighbourhoods)$col_val[idxs] <-
            apply(nhoodCounts(agg_neighbourhoods), 1,
              function(x) x[ci] / sum(x)) / sums[ci]
          p <- plotNhoodGraph(agg_neighbourhoods,
            colour_by = 'col_val') +
            ggtitle(colnames(nhoodCounts(agg_neighbourhoods))[ci])
          # print_plot_eval(
          #   print(),
          #   width = 17.4, height = 15,
          #   filename = of
          # )
          # return(of)
          return(p)
        })
      test_mode <- 'stim_group_only'
      experiment <- colData(agg_neighbourhoods)$exp[1]
      o_fn <- file.path(Sys.getenv('img_dir'),
        glue::glue('milo_embedding-CN_frac-{experiment}\\
          {make_flag(test_mode)}.pdf'))
      maartenutils::plot_panel_layout(
        plots = plots,
        labels = NULL,
        plot_direct = test_rendering(),
        nrow = 2, ncol = 1,
        filename = o_fn)

      ## Plot DA results
      plot_DA_Nhoods(
        sce = agg_neighbourhoods,
        DA_Nhoods = NH_test_DA,
        test_mode = 'stim_group_only',
        # fn_app = paste(genelist, k, d, sep='-')
        fn_app = ''
      )
    },
    format = 'file'
  ),

  # tar_target(
  #   name = NH_so,
  #   command = {
  #     so <- tar_read(filtered_so_6743)

  #     if (interactive())
  #       source('~/MirjamHoekstra/R/init.R')
  #     NH_so <- gen_milo_so(
  #       so = so,
  #       sce = agg_neighbourhoods,
  #       min_neighbourhood_size = 0
  #     )

  #     return(NH_so)
  #   },
  #   iteration = 'list'
  # ),

  tar_target(
    name = so_score,
    command = {

      ## Investigate are driven by the neighbourhood level
      ## pseudobulking or by the now included genes
      MR_clusters_mouse <- tar_read(MR_clusters_mouse)
      mouse_genesets <- clust2list(MR_clusters_mouse)

      library(genesets)
      HM_mouse <- filter_gmt(pattern='.*', 
        gmt_pattern='mh.all.v0.3', 
        type='symbols')
      names(HM_mouse) <- stringr::str_replace(
        names(HM_mouse), 'HALLMARK_', '')
      names(HM_mouse) <- stringr::str_replace_all(
        names(HM_mouse), '_', ' ')

      so <- tar_read(filtered_so_6743)

      if (interactive())
        source('~/MirjamHoekstra/R/init.R')
      so_score <- ComputeGeneSetScores(
        so = so,
        assay = 'SCT',
        datatype = 'counts',
        weighting_funcs = all_weighting_funcs['unweighted'],
        genesets = c(mouse_genesets, HM_mouse, 
          mirjam_genes_220516, mirjam_genes_220518, 
          mirjam_genes_220519),
        # tV = apply(tM_human, 1, max),
        simplify_names = F
      )
      so_score <-
        min_max_scale_assay(so_score, assay = 'GS') %>%
        { . }

      return(so_score)
    },
    iteration = 'list'
  ),

  # tar_target(
  #   name = NH_gamma_titration_cosine,

  #   command = {
  #     gamma_tit <- perform_TRANSACT_gamma_titration(
  #       reference_M = NH_harmonized_Ms_obj[['reference']],
  #       query_M = NH_harmonized_Ms_obj[['query']],
  #       TRANSACT_params = list(
  #         N_source_components = N_source_components,
  #         N_target_components = N_target_components,
  #         N_PV = min(N_PV, N_source_components, N_target_components),
  #         kernel = 'mallow'
  #       ),
  #       gamma_resolution = 0.5,
  #       # gamma_resolution = 4,
  #       gamma_bounds = c(-7, -3),
  #       ncores = 1L
  #     )
  #     return(gamma_tit)
  #   }
  # ),

  # tar_target(
  #   name = so_score_hallmark,
  #   command = {

  #     ## Investigate are driven by the neighbourhood level
  #     ## pseudobulking or by the now included genes
  #     pathways <- filter_gmt(pattern='.*', 
  #       gmt_pattern='mh.all.v0.3', 
  #       type='symbols')
  #     names(pathways) <- stringr::str_replace(names(pathways), 'HALLMARK_', '')
  #     names(pathways) <- stringr::str_replace_all(names(pathways), '_', ' ')
  #     so <- tar_read(filtered_so_6743)

  #     if (interactive())
  #       source('~/MirjamHoekstra/R/init.R')
  #     so_score_hallmark <- ComputeGeneSetScores(
  #       so = so,
  #       assay = 'SCT',
  #       datatype = 'counts',
  #       weighting_funcs = all_weighting_funcs['unweighted'],
  #       genesets = pathways,
  #       simplify_names = F
  #     )
  #     so_score_hallmark <-
  #       min_max_scale_assay(so_score_hallmark, assay = 'GS') %>%
  #       { . }

  #     return(so_score_hallmark)
  #   },
  #   iteration = 'list'
  # ),

  tar_target(
    name = so_score_mirjam_set,
    command = {
      so <- tar_read(filtered_so_6743)
      so_s <- ComputeGeneSetScores(
        so = so,
        assay = 'SCT',
        datatype = 'counts',
        weighting_funcs = all_weighting_funcs['unweighted'],
        genesets = c(mirjam_genes_220516, mirjam_genes_220518, 
          mirjam_genes_220519),
        simplify_names = T
      )
      so_s <-
        min_max_scale_assay(so_s, assay = 'GS') %>%
        { . }
      return(so_s)
    },
    iteration = 'list'
  ),

  # tar_target(
  #   name = NH_so_score,
  #   command = {

  #     mouse_genesets <- clust2list(MR_clusters_mouse)

  #     if (interactive() && !test_rendering())
  #       source('~/MirjamHoekstra/R/init.R')
  #     NH_so_score <- ComputeGeneSetScores(
  #       so = NH_so,
  #       assay = 'RNA',
  #       datatype = 'counts',
  #       weighting_funcs = all_weighting_funcs['unweighted'],
  #       genesets = mouse_genesets,
  #       # tV = apply(tM_human, 1, max),
  #       simplify_names = F
  #     )
  #     NH_so_score <-
  #       min_max_scale_assay(NH_so_score, assay = 'GS') %>%
  #       { . }

  #     return(NH_so_score)
  #   },
  #   iteration = 'list'
  # ),

  # tar_target(
  #   name = NH_so_umaps,
  #   command = {
  #     library(miloR)
  #     library(SummarizedExperiment)
  #     agg_neighbourhoods <- tar_read(agg_neighbourhoods_mouse)

  #     ## First plot geneset scores on top of neighbourhood graph
  #     fns <-
  #       rownames(NH_so_score[['GS']]) %>%
  #       purrr::map(function(fn) {
  #         # fn <- rownames(NH_so_score[['GS']])[1]
  #         of <- file.path(Sys.getenv('img_dir'),
  #           glue::glue('Nhoods_{fn}.pdf'))

  #         idxs <- unlist(nhoodIndex(agg_neighbourhoods))
  #         colData(agg_neighbourhoods)$col_val <- NA_real_
  #         colData(agg_neighbourhoods)$col_val[idxs] <-
  #           NH_so_score[['GS']][fn, as.character(idxs)]

  #         print_plot_eval(
  #           print(
  #             plotNhoodGraph(agg_neighbourhoods,
  #               colour_by = 'col_val')),
  #           width = 17.4, height = 15,
  #           filename = of
  #         )

  #         return(of)
  #       })

  #     ## Next plot individual genes on top of neighbourhood graph
  #     limma_MRs <- tar_read(limma_MRs_mouse)
  #     fns <- map(limma_MRs, function(fn) {
  #       colData(agg_neighbourhoods)$col_val <-
  #         NH_so_score[['RNA']][fn, ]
  #       colData(agg_neighbourhoods)$col_val <-
  #         as.numeric(colData(agg_neighbourhoods)$col_val)
  #       of <- file.path(
  #         Sys.getenv('img_dir'), glue::glue('Nhoods_{fn}.pdf')
  #       )
  #       p <- plotNhoodGraph(agg_neighbourhoods, colour_by = 'col_val')
  #       print_plot_eval(print(p),
  #         width = 17.4, height = 15, filename = of)
  #       return(of)
  #     })

  #     return(fns[[length(fns)]])
  #   }
  # ),

  tar_target(
    name = SC_score_vln,
    command = {
      # so <- tar_read(filtered_cleaned_so_6743) %>%
      #   format_6743()
      so_score <- tar_read(so_score_mirjam_set_mouse) %>%
        format_6743()
      # so_score@meta.data$stim_group

      p_dat <- assay2dtf(so_score, assay = 'GS') %>%
        format_6743()
      # p_dat$stim_group
      fns <- rownames(so_score@assays$GS[,])
      p_dat$grp_i <- as.integer(p_dat$stim_group) %>%
        factor(., levels = rev(sort(unique(.))))

      if (interactive() && !test_rendering())
        source('~/MirjamHoekstra/R/init.R')
      plots <- fns %>%
        c(.[1], .) %>%
        purrr::imap(function(fn, i) {
          p <- 
            mk_ann_GS_vln(
              GS_p_dat = p_dat,
              fn = fn,
              cn_mode = 'all',
              experiment = '6743',
              p_values = 'all',
              # return_mode = 'plot',
              plot_lgd = i == 1
            ) +
            scale_x_discrete(breaks = c()) +
            theme()
          return(p)
        })

      experiment <- so_score@meta.data$exp[1]
      o_fn <- file.path(
        Sys.getenv('img_dir'),
        glue::glue('SC_GS_scores-{experiment}.pdf')
      )
      maartenutils::plot_panel_layout(
        plots = plots,
        labels = NULL,
        plot_direct = test_rendering(),
        nrow = 4, ncol = 1,
        filename = o_fn
      )

      return(o_fn)
    }
  ),

  # tar_target(
  #   name = SC_score_vln_hallmark,
  #   command = {

  #     p_dat <- assay2dtf(so_score_hallmark, assay = 'GS')
  #     fns <- rownames(so_score_hallmark@assays$GS[,])
  #     plots <- map(fns, function(fn) {
  #       p <-
  #         plot_vln(
  #           p_dat = p_dat,
  #           fn = fn,
  #           group_var = 'stim_group',
  #           fill_var = 'stim_group',
  #           assay = 'GS',
  #           p_values = 'reference'
  #         ) +
  #         scale_x_discrete(breaks = c()) +
  #         theme()
  #       return(p)
  #     })

  #     sa <- so_score_hallmark@meta.data[, c('stim_group'), drop = F]
  #     M <- so_score_hallmark@assays$GS[,]
  #     M <- t(scale(t(M)))
  #     M[M > 5] <- 5
  #     M[M < -5] <- -5
  #     des <- stringr::str_subset(rownames(M), 
  #       'INTERFERON|TNF|APOP|TGF|WNT|IL2')
  #     M <- M[des, ]
  #     rn <- rownames(M)
  #     rn <- stringr::str_replace(rn, 'HALLMARK_', '')
  #     rn <- stringr::str_replace_all(rn, '_', ' ')
  #     rn <- stringr::str_replace_all(rn, '.vanilla.unweighted.sum', ' ')
  #     rownames(M) <- rn
  #     HM <- gen_HM(M, ca = sa, name = 'GS score')
  #     print_plot_eval(
  #       draw(HM, merge_legend = T, heatmap_legend_side = 'bottom'),
  #       width = 17.4, height = 15,
  #       filename = file.path(Sys.getenv('img_dir'),
  #         glue::glue('mouse_hallmark_HM.pdf')))

  #     experiment <- so_score_hallmark@meta.data$exp[1]
  #     o_fn <- file.path(Sys.getenv('img_dir'),
  #       glue::glue('SC_Hallmark_GS_scores-{experiment}.pdf'))
  #     maartenutils::plot_panel_layout(
  #       plots = plots,
  #       labels = NULL,
  #       plot_direct = test_rendering(),
  #       nrow = 2, ncol = 1,
  #       filename = o_fn
  #     )

  #     return(o_fn)
  #   }
  # ),

  # tar_target(
  #   name = ,
  #   command = {
  #     so <- tar_read(filtered_so_6743)
  #   }
  # ),

  # tar_target(
  #   name = ,
  #   command = {
  #     so_s[['sga']] <- so_s[['stim_group']]

  #     fns <- rownames(so_s[['GS_MM']])
  #     # plots <-
  #     purrr::walk(fns, function(fn) {
  #       p <- plot_vln(so_s, fn,
  #         # assay = 'GS_MM',
  #         assay = 'GS',
  #         group_var = 'sga',
  #         pt.size = .02) +
  #         scale_x_discrete(expand = c(0, 0)) +
  #         ggtitle(with(long_gene_cluster_names, gcn[which(gc == fn)])) +
  #         ggpubr::stat_compare_means(
  #           label = "p.signif",
  #           method = "wilcox.test",
  #           ref.group = 'Mix PBS',
  #           label.x.npc = .5,
  #           label.y.npc = 1,
  #           colour = 'indianred3'
  #         )
  #       pe <- rlang::expr({ print(p) })
  #       print_plot_eval(!!pe, width = 12, height = 10,
  #         filename = file.path(Sys.getenv('img_dir'),
  #           glue::glue('violin_{fn}.pdf')))
  #     })
  #   }
  # ),

  NULL
)
