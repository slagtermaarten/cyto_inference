murine_targets <- list(
  tar_map(
    names = type,
    values = tibble(type = 'mouse'),

    tar_target(
      name = expM,
      command = {
        expM <-
          readr::read_tsv(file.path(data_dir, 'raw_exp_5892', 'readcounts.txt')) %>%
          rename_with(~basename(.x), matches('5892')) %>%
          rename_with(~gsub('5892_\\d{1,2}_(.*)_[ATCG]{7}_.*', '\\1', .x),
                      matches('5892')) %>%
          dplyr::select(-ensembl_gene_id, -start_position, -end_position,
            -gene_biotype, -chromosome_name, -description) %>%
          dplyr::rename(gene = external_gene_id)
        return(expM)
      }
    ),

    tar_target(
      name = sample_annotation,
      command = {
        num_cols <- tolower(setdiff(colnames(expM), 'gene'))
        duration <- factor(as.numeric(gsub('(\\d{1,2})h_.*$', '\\1', num_cols)))
        tnf_conc <- gsub('\\d{1,2}h.*_(\\d{1,2}\\.*\\d*)_ng_ml_tnfa$', '\\1', num_cols)
        tnf_conc[tnf_conc == num_cols] <- 0
        ifn_conc <- gsub('\\d{1,2}h.*_(\\d{1,2}\\.*0*1*)_ng_ml_ifn.*', '\\1', num_cols)
        ifn_conc[ifn_conc == num_cols] <- 0
        lt_conc <- gsub('\\d{1,2}h.*_(\\d{1,2})ng_ml_lymphotoxin_a1_b2', '\\1', num_cols)
        lt_conc[lt_conc == num_cols] <- 0
        tgf_conc <- gsub('\\d{1,2}h.*_(\\d{1,2})_ng_ml_tgfb', '\\1', num_cols)
        tgf_conc[tgf_conc == num_cols] <- 0
        ifna_conc <- gsub('\\d{1,2}h.*_(\\d{1,9})_u_ml_ifna1', '\\1', num_cols)
        ifna_conc[ifna_conc == num_cols] <- 0
        il2_conc <- gsub('\\d{1,2}h.*_(\\d{1,9})_u_ml_il-2', '\\1', num_cols)
        il2_conc[il2_conc == num_cols] <- 0
        stim_groups <- gsub('_', ' ', num_cols) %>%
          { gsub('ng ml', 'ng/ml', .) } %>%
          { gsub('(\\d) (\\d+) sn', '\\1/\\2 SN', .) } %>%
          { gsub('tnfa', 'TNFa', .) } %>%
          { gsub('il-2', 'IL-2', .) } %>%
          { gsub('a1 b2', 'A1 B2', .) } %>%
          { gsub('ifna1', 'IFNA1', .) } %>%
          { gsub('ifn[y|g]', 'IFNy', .) } %>%
          { gsub('unstim', 'Unstimulated in vitro', .) } %>%
          { gsub('^\\d{1,2}h ', '', .) } %>%
          { gsub('u ml', 'U/ml', .) } %>%
          { gsub('2ng/ml', '2 ng/ml', .) } %>%
          { gsub('0 01', '0.01', .) } %>%
          { gsub('0 1', '0.1', .) } %>%
          { gsub('tgfb', 'TGFb', .) } %>%
          { . }
        sample_annotation <-
          data.frame(
            'sample_name' = setdiff(colnames(expM), 'gene'),
            'duration' = factor(factor_to_numeric(duration)),
            'stim_group' = stim_groups,
            'tnf_conc' = factor(factor_to_numeric(tnf_conc)),
            'ifn_conc' = factor(factor_to_numeric(ifn_conc)),
            'lt_conc' = factor(factor_to_numeric(lt_conc)),
            'tgf_conc' = factor(factor_to_numeric(tgf_conc)),
            'ifna_conc' = factor(factor_to_numeric(ifna_conc)),
            'il2_conc' = factor(factor_to_numeric(il2_conc)),
            sample_origin = include_na_lev(v = 'in_vitro', el = 'in_vivo'),
            sample_type = include_na_lev('bulk', el = 'sc')
          ) %>%
          dplyr::arrange(duration, ifn_conc, tnf_conc, lt_conc, tgf_conc,
            ifna_conc, il2_conc) %>%
          dplyr::mutate(stim_group = factor(stim_group,
              levels = unique(stim_group))) %>%
          add_condition_name()
        return(sample_annotation)
      }
    ),

    tar_target(
      name = M,
      command = {
        M <- as.matrix(dplyr::select(expM, -gene)) %>%
          { .[, match(
            sample_annotation$sample_name, 
            colnames(.)
           )] } %>%
          { set_colnames(., tolower(colnames(.))) } %>%
          # { . }
          set_rownames(expM$gene) %>%
          edgeR::DGEList() %>%
          edgeR::calcNormFactors(method = 'TMM')
        return(M)
      }
    ),

    tar_target(bulk_5892_so, {
      sample_annotation <- tar_read(sample_annotation_mouse)
      M <- tar_read(M_mouse)
      stopifnot(tolower(sample_annotation$sample_name) == colnames(M))
      colnames(M) <- sample_annotation$sample_name
      colnames(M$counts)
      exp5892 <- CreateSeuratObject(
        counts = M$counts,
        meta.data = sample_annotation %>%
          as.data.frame() %>%
          set_rownames(NULL) %>%
          dplyr::mutate(sample_name = tolower(sample_name)) %>%
          tibble::column_to_rownames('sample_name')
      )
      # exp5892@meta.data
      exp5892 <- SCTransform(exp5892)
      exp5892 <- FindVariableFeatures(exp5892)
      exp5892 <- RunPCA(exp5892, assay = 'SCT',
        npcs = 20, verbose = F)
      return(exp5892)
    }, iteration = 'list'),

    tar_target(
      name = M_cpm,
      command = {
        M_cpm <- edgeR::cpm(M)
        colnames(M_cpm) <- sample_annotation$sample_name
        return(M_cpm)
      }
    ),

    tar_target(
      name = gene_var,
      command = {
        library(scran)
        gene_var <- modelGeneVar(log2(M_cpm + 1),
          block = t(sample_annotation[, 'duration', drop = F])
        )
        gene_var <- as.data.frame(gene_var)
        return(gene_var)
      }
    ),

    tar_target(
      name = f_gene_var,
      command = {
        desirable <- gene_var$bio > 0 & gene_var$FDR <= 1e-7
        gene_var[desirable, ]
      }
    ),

    tar_target(
      name = example_fano_lineplots,
      command = {
        # f_gene_var <- tar_read(f_gene_var_mouse)
        # sample_annotation <- tar_read(sample_annotation_mouse)
        # M_cpm <- tar_read(M_cpm_mouse)

        fns <- exp_dynamics_panel(
          features = rownames(f_gene_var) %>%
            { sample(., 4*8L-1) },
          meta_data = sample_annotation,
          merge_cn = 'sample_name',
          version = 'sample_fano',
          lookup_data = M_cpm,
          redo = F,
          leave_out_sn = T
        )

        return(fns[[length(fns)]])
      },
      format = 'file'
    ),

    tar_target(
      name = efit,
      command = {
        library(limma)
        library(edgeR)
        d0 <- M
        gene_max_e <- apply(edgeR::cpm(d0), 1, max)
        gene_sd <- apply(edgeR::cpm(d0), 1, sd)
        gene_max_e_rc <- apply(M, 1, max)
        drop <- which(gene_max_e < 25 | sqrt(gene_sd) < 0.05)
        mm <- with(sample_annotation,
          model.matrix(~ duration + stim_group))
        d0 <- d0[-drop, ]
        print_plot_eval(
          { v <<- limma::voom(d0, mm, plot=T, span = .3) },
          filename = file.path(Sys.getenv('img_dir'),
            glue::glue('mouse_limma.png'))
        )
        vfit <- lmFit(v, mm)
        efit <- eBayes(vfit)
        return(efit)
      }
    ),

    tar_target(
      name = tM,
      command = {
        efit <- tar_read(efit_mouse)
        tM <- efit$t
        tM <- tM[, paste0('stim_group', c('6000 U/ml IL-2',
              '5000 U/ml IFNA1', '10 ng/ml TGFb',
              '2 ng/ml lymphotoxin A1 B2', '10 ng/ml TNFa',
              '100 ng/ml IFNy', '100 ng/ml IFNy 10 ng/ml TNFa'))]
      }
    ),

    tar_target(
      name = tM_HM,
      command = {
        tM <- tar_read(tM_mouse)
        colnames(tM) <- stringr::str_replace(colnames(tM),
          'stim_group', '')
        of <- file.path(Sys.getenv('img_dir'),
            glue::glue('mouse_tM.pdf'))
        max_t <- apply(tM, 1, max)
        idxs <- which(max_t >= 3)
        tM_s <- tM[idxs, , drop=F]

        HM <- GE_vs_stim_dur(
          M = tM_s,
          cell_width_mm = 10,
          cell_height_mm = 10*(20-4)/nrow(tM_s),
          scale_data = F,
          N_hl_genes = 30,
          row_dist_f = 'euclidean',
          row_cluster_k = 6,
          split_cols = FALSE,
          # cluster_method = 'complete',
          show_column_names = TRUE,
          # show_column_dend = FALSE,
          # N_hl_genes = 10,
          name = 't-value'
        )

        if (FALSE) {
          M_cpm <- tar_read(M_cpm_mouse)
          # cl_ <- 
          #   sort(gen_tree_split(gene_clust_obj, cluster_k = 13)) %>%
          #   { .[intersect(names(.), attr(HM, 'highlight_genes'))] }
          fns <- exp_dynamics_panel(
            features = attr(HM, 'highlight_genes'),
            meta_data = tar_read(sample_annotation_mouse),
            merge_cn = 'sample_name',
            version = 'murine_limma_t_hl_genes',
            lookup_data = M_cpm,
            redo = F,
            leave_out_sn = F
          )
        }

        print_plot_eval(
          draw(HM, heatmap_legend_side = 'top'),
          width = 17.4, height = 20,
          filename = of
        )

        return(of)
      },
      format = 'file'
    ),
  
    tar_target(
      name = limma_MR_params,
      command = {
        limma_MR_params <- list(LL = 4.2, UL = 1, SD = .5, fd = 4)
        return(limma_MR_params)
      }
    ),

    tar_target(
      name = stringent_limma_MRs,
      command = {
        tM <- tar_read(tM_mouse)

        limma_MR_params <- tar_read(limma_MR_params_mouse)
        RC  <- tM >= limma_MR_params$LL
        NRC <- tM < limma_MR_params$UL
        MRs <- names(which(
          apply(RC[, 1:6], 1, sum) == 1 &
          apply(NRC[, 1:6], 1, sum) >= 5 &
          abs(apply(tM[, 5:6], 1, sum) - tM[, 7]) <= limma_MR_params$SD
        ))
        IFNy_MRs <- 
          names(which(
            RC[, 6] &
            apply(NRC[, c(1, 3:6)], 1, sum) >= 4 &
            abs(apply(tM[, 5:6], 1, sum) - tM[, 7]) <= limma_MR_params$SD
          ))
        # easy_MRs <- 
        #   names(which(
        #     apply(RC[, c(1, 3:6)], 1, sum) == 1 &
        #     apply(NRC[, c(1, 3:6)], 1, sum) >= 4 &
        #     abs(apply(tM[, 5:6], 1, sum) - tM[, 7]) <= limma_MR_params$SD
        #   ))

        synergy_genes <-
          names(which(tM[, 7] >
              limma_MR_params$fd * apply(abs(tM[, 5:6]), 1, sum)))

        f_gene_var <- tar_read(f_gene_var_mouse)
        stringent_limma_MRs <- c(IFNy_MRs, MRs, synergy_genes) %>%
          intersect(rownames(f_gene_var))
        # list('IFNy' = IFNy_MRs, 'MRs' = MRs, 'synergy' = synergy_genes)
        return(stringent_limma_MRs)
      }
    ),

    tar_target(
      name = stringent_MR_clusters,
      command = {
        M_cpm <- tar_read(M_cpm_mouse)
        lMs <- get_lMs(M = M_cpm, genes = stringent_limma_MRs)
        clust_obj <-
          gen_clust_object(t(lMs),
            dist_f = 'pearson',
            clust_method = 'complete'
          )
        stringent_MR_clusters <- 
          gen_tree_split(clust_obj, cluster_k = 13) %>%
          sort()
        # length(tar_read(MR_clusters_mouse))
        return(stringent_MR_clusters)
      }
    ),

    tar_target(
      name = HM,
      command = {
        out_dir <- file.path(Sys.getenv('img_dir'), 'mouse')
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

        for (ref_mode in c('all', 'non-homeopathic', 'tnfa_ifny')) {
          for (gene_sel in c('mirjam_mouse', 'stringent_limma_MRs',
              'most_var')) {
          # gene_sel = 'mirjam_mouse'
          # gene_sel = 'most_var'
          # ref_mode = 'tnfa_ifny'

          sa <- tar_read(sample_annotation_mouse) %>%
            dplyr::select(sample_name, duration, stim_group) %>%
            order_stim_group() %>%
            dplyr::rename(m_sg = stim_group)
          if (ref_mode == 'non-homeopathic') {
            idxs <- which(!stringr::str_detect(sa$m_sg, 'IL-2|A1 B2'))
          } else if (ref_mode == 'tnfa_ifny') {
            idxs <- which(!stringr::str_detect(sa$m_sg, 
                'IL-2|A1 B2|TGFb|IFNA1'))
          } else if (ref_mode == 'all') {
            idxs <- 1:nrow(sa)
          }
          sa <- sa[idxs, ]
          sa$m_sg <- droplevels(sa$m_sg)
          M_cpm <- tar_read(M_cpm_mouse) %>% { .[, idxs] }

          if (gene_sel == 'stringent_limma_MRs') {
            M_cpm <- M_cpm %>%
              subset_feats(stringent_limma_MRs)
          } else if (gene_sel == 'mirjam_mouse') {
            M_cpm <- M_cpm %>%
              subset_feats(unlist(mouse_genesets))
          } else if (gene_sel == 'most_var') {
            f_gene_var <- tar_read(f_gene_var_mouse)
            M_cpm <- M_cpm %>%
              subset_feats(rownames(f_gene_var))
          }

          gene_clust_obj <-
            #  %>%
            M_cpm %>%
            t() %>%
            gen_clust_object(
              dist_f = 'euclidean',
              clust_method = 'complete'
            )

          if (interactive() && !test_rendering())
            source('~/MirjamHoekstra/R/init.R')
          of <- file.path(out_dir,
            glue::glue('mouse_MR{make_flag(ref_mode)}{make_flag(gene_sel)}.pdf'))
          stopifnot(colnames(M_cpm) == sa$sample_name)
          HM <- GE_vs_stim_dur(
            M = as.matrix(M_cpm),
            cell_width_mm = 1.5,
            sa = sa %>%
              as.data.frame() %>%
              {
                set_rownames(dplyr::select(., -sample_name), .$sample_name)
              },
            cell_height_mm = 10*(20-6)/nrow(M_cpm),
            scale_data = T,
            # row_clust_object = gene_clust_obj,
            row_dist_f = 'euclidean',
            row_cluster_k = 4,
            col_cluster_k = 4,
            # split_cols = TRUE,
            # split_rows = TRUE,
            split_cols = FALSE,
            split_rows = FALSE,
            row_clust_method = 'complete',
            show_column_names = FALSE
            # show_column_dend = FALSE,
            # N_hl_genes = 10,
            # name = 'Z-scaled\ngene expression'
          )
          print_plot_eval(
            draw(HM, heatmap_legend_side = 'top'),
            width = 17.4, height = 20,
            filename = of
          )

          if (TRUE) {
            if (interactive())
              source('~/MirjamHoekstra/R/init.R')
            M_cpm <- tar_read(M_cpm_mouse)
            cl_ <- 
              sort(gen_tree_split(gene_clust_obj, cluster_k = 13)) %>%
              { .[intersect(names(.), attr(HM, 'highlight_genes'))] }
            fns <- exp_dynamics_panel(
              features = names(cl_),
              plot_titles = paste0(cl_, ' - ', names(cl_)),
              meta_data = tar_read(sample_annotation_mouse),
              merge_cn = 'sample_name',
              # y_scale = 'linear',
              version = glue::glue('murine_limma_t_hl_genes\\
                {make_flag(ref_mode)}{make_flag(gene_sel)}'),
              lookup_data = M_cpm,
              out_dir = out_dir,
              redo = T,
              leave_out_sn = F
            )
          }

        }}
        return(of)
      },
      format = 'file'
    ),

    # tar_target(
    #   name = MR_cluster_gsea,
    #   command = {
    #     library(genesets)
    #     library(fgsea)

    #     fns <- pmap(colnames(tM), function(sg) {
    #       # sg <- colnames(tM)[1]
    #       ranks <- tM[, sg]
    #       cn_pf <- stringr::str_replace_all(sg, '/| ', '_')
        
    #       pathways <- filter_gmt(pattern='.*', 
    #         gmt_pattern='mh.all.v0.3', 
    #         type='symbols')
    #       names(pathways) <- stringr::str_replace(names(pathways), 'HALLMARK_', '')
    #       names(pathways) <- stringr::str_replace_all(names(pathways), '_', ' ')

    #       fgseaRes <-
    #         fgsea::fgsea(pathways = pathways, stats = ranks) %>%
    #         dplyr::arrange(desc(NES))
        
    #       print_plot_eval(
    #         {
    #           tgrob <- fgsea::plotGseaTable(
    #             pathways = pathways[fgseaRes$pathway],
    #             fgseaRes = fgseaRes, 
    #             stats = ranks, 
    #             render = F
    #           )
    #           plot(tgrob)
    #         },
    #         width = 17.4, height = 25,
    #         filename = file.path(Sys.getenv('img_dir'),
    #           glue::glue('mouse_gsea_tab_{cn_pf}.pdf'))
    #       )

    #     })
    #     return(fns[[length(fns)]])
    #   },
    #   format = 'file'
    # ),

    # tar_target(
    #   name = MR_cluster_fisher,
    #   command = {
    #     mouse_genesets <- clust2list(stringent_MR_clusters_mouse)
    #     fns <- map(mouse_genesets, function(gs) {
    #       ranks <- tM[, cn]
    #       if (F) {
    #         ranks <- ranks[ranks >= min(tV[names(co)])]
    #       }
    #       cn_pf <- stringr::str_replace_all(cn, '/| ', '_')
        
    #       fgseaRes <-
    #         fgsea::fgsea(pathways = pathways, stats = ranks) %>%
    #         dplyr::arrange(desc(NES))
        
    #       print_plot_eval(
    #         {
    #           tgrob <- fgsea::plotGseaTable(
    #             pathways = pathways[fgseaRes$pathway],
    #             fgseaRes = fgseaRes, 
    #             stats = ranks, 
    #             render = F
    #           )
    #           plot(tgrob)
    #         },
    #         width = 17.4, height = 25,
    #         filename = file.path(Sys.getenv('img_dir'),
    #           glue::glue('gsea_tab_{cn_pf}.pdf'))
    #       )
    #     })
    #     return(fns[[length(fns)]])
    #   },
    #   format = 'file'
    # ),

    # tar_target(
    #   name = limma_synergy,
    #   command = {
    #     tM <- tar_read(tM_mouse)
    #     names(which(tM[, 7] > limma_MR_params$fd * apply(abs(tM[, 5:6]), 1, sum)))
    #   }
    # ),

    tar_target(
      name = gene_stats,
      command = {
        source(file.path(Sys.getenv('r_dir'), 'mouse.R'))
        extract_gene_stats(efit = efit, M_cpm = M_cpm)
      }
    ),

    tar_target(
      name = limma_genes,
      command = {
        gene_stats <- tar_read(gene_stats_mouse)
        gene_stats %>%
          dplyr::filter(min_p <= .1) %>%
          dplyr::filter(max_estimate >= .5) %>%
          # dplyr::filter(gene %in% rownames(f_gene_var)) %>%
          { . }
      }
    ),

    tar_target(
      name = mono_reporters,
      command = {

        # summary(tar_read(limma_genes_mouse)$synergy)
        limma_genes <- tar_read(gene_stats_mouse)
        # limma_genes %>%
        #   tally(tgfb_mr_l, tnfa_neg_l, ifny_neg_l)
        mono_reporters <- limma_genes %>%
          dplyr::filter(min_p <= .1) %>%
          dplyr::filter(max_estimate >= .5) %>%
          # dplyr::filter(tnfa_mr | ifny_mr | ifny_mr_l | tnfa_mr_l) %>%
          dplyr::filter(ifny_mr_l | tnfa_mr_l | 
            tgfb_mr_l | synergy >= 2) %>%
          # dplyr::arrange(tnfa_mr, ifny_mr, ifny_mr_l, tnfa_mr_l)
          dplyr::arrange(ifny_mr_l, tnfa_mr_l, tgfb_mr_l)

        return(mono_reporters)
      }
    ),

    tar_target(
      name = mono_reporters_HM,
      command = {

        for (include_homeopathic in c(T, F)) {
          # include_homeopathic = T

          sa <- 
            tar_read(sample_annotation_mouse) %>%
            dplyr::select(sample_name, duration, stim_group) %>%
            order_stim_group() %>%
            dplyr::rename(m_sg = stim_group)
          if (!exists('mono_reporters')) {
            mono_reporters <- tar_read(mono_reporters_mouse)
          }
          # mono_reporters$tgfb_mr_l
          if (!include_homeopathic) {
            idxs <- which(!stringr::str_detect(sa$m_sg, 'IL-2|A1 B2'))
          } else {
            idxs <- 1:nrow(sa)
          }
          sa <- sa[idxs, ]
          M_cpm <- tar_read(M_cpm_mouse) %>%
            subset_feats(mono_reporters$gene) %>%
            { .[, idxs] }
          sa$m_sg <- droplevels(sa$m_sg)

          gene_clust_obj <-
            #  %>%
            M_cpm %>%
            t() %>%
            gen_clust_object(
              dist_f = 'pearson',
              clust_method = 'complete'
            )
          # str(gene_clust_obj)
          # gene_clust_obj <- as.hclust(gene_clust_obj)
          # gene_clust_obj$height
          # gene_clust_obj

          if (interactive() && !test_rendering())
            source('~/MirjamHoekstra/R/init.R')
          of <- file.path(Sys.getenv('img_dir'),
            glue::glue('mouse_MR{make_flag(include_homeopathic)}.pdf'))
          stopifnot(colnames(M_cpm) == sa$sample_name)
          HM <- GE_vs_stim_dur(
            M = M_cpm,
            cell_width_mm = 1.5,
            sa = dplyr::select(sa, -sample_name),
            cell_height_mm = 10*(20-4)/nrow(M_cpm),
            scale_data = T,
            row_clust_object = gene_clust_obj,
            row_dist_f = 'euclidean',
            row_cluster_k = 13,
            split_cols = FALSE,
            row_clust_method = 'complete',
            show_column_names = FALSE
            # show_column_dend = FALSE,
            # N_hl_genes = 10,
            # name = 'Z-scaled\ngene expression'
          )
          print_plot_eval(
            draw(HM, heatmap_legend_side = 'right', 
              merge_legends = TRUE),
            width = 17.4, height = 20,
            filename = of
          )

          if (TRUE) {
            if (interactive())
              source('~/MirjamHoekstra/R/init.R')
            M_cpm <- tar_read(M_cpm_mouse)
            cl_ <- 
              sort(gen_tree_split(gene_clust_obj, cluster_k = 13)) %>%
              { .[intersect(names(.), attr(HM, 'highlight_genes'))] }
            fns <- exp_dynamics_panel(
              features = names(cl_),
              plot_titles = paste0(cl_, ' - ', names(cl_)),
              meta_data = tar_read(sample_annotation_mouse),
              merge_cn = 'sample_name',
              # y_scale = 'linear',
              version = glue::glue('murine_limma_t_hl_genes\\
                {make_flag(include_homeopathic)}'),
              lookup_data = M_cpm,
              redo = T,
              leave_out_sn = F
            )
          }

        }

        return(of)
      },
      format = 'file'
    ),

    tar_target(
      name = MR_clusters,
      command = {
        M_cpm <- tar_read(M_cpm_mouse)
        mono_reporters_mouse <- tar_read(mono_reporters_mouse)
        lMs <- get_lMs(M = M_cpm, genes = mono_reporters_mouse$gene)
        clust_obj <-
          gen_clust_object(t(lMs),
            dist_f = 'pearson',
            clust_method = 'complete'
          )
        MR_clusters <- 
          gen_tree_split(clust_obj, cluster_k = 13) %>%
          sort()
        # length(tar_read(MR_clusters_mouse))
        return(MR_clusters)
      }
    ),


    tar_target(
      name = monoreporter_lineplots,
      command = {
        fns <- exp_dynamics_panel(
          features = mono_reporters$gene %>%
            intersect(rownames(f_gene_var)),
          meta_data = sample_annotation,
          merge_cn = 'sample_name',
          version = 'murine_mono_reporters',
          lookup_data = M_cpm,
          redo = F,
          leave_out_sn = T
        )
        return(fns[[length(fns)]])
      },
      format = 'file'
    ),

    tar_target(
      name = monoreporter_lineplots_non_fano,
      command = {
        fns <- exp_dynamics_panel(
          features = mono_reporters$gene %>%
            setdiff(rownames(f_gene_var)),
          meta_data = sample_annotation,
          merge_cn = 'sample_name',
          version = 'murine_mono_reporters_non_fano',
          lookup_data = M_cpm,
          redo = F,
          leave_out_sn = T
        )
        return(fns[[length(fns)]])
      },
      format = 'file'
    ),

    # tar_target(
    #   name = auto_genes,
    #   command = {
    #     intersect(limma_genes$gene, rownames(f_gene_var))
    #   }
    # ),

    # tar_target(
    #   name = synergy_genes,
    #   command = {
    #     limma_genes %>%
    #       dplyr::filter(synergy > 0) %>%
    #       dplyr::filter(ifny > 0 & tnfa > 0) %>%
    #       { . }
    #   }
    # ),

    # tar_target(
    #   name = synergy_lineplots,
    #   command = {
    #     fns <- exp_dynamics_panel(
    #       features = synergy_genes$gene,
    #       meta_data = sample_annotation,
    #       merge_cn = 'sample_name',
    #       version = 'murine_synergy',
    #       lookup_data = M_cpm,
    #       redo = F,
    #       leave_out_sn = T
    #     )
    #     return(fns[[length(fns)]])
    #   },
    #   format = 'file'
    # ),

    # tar_target(
    #   name = sample_group,
    #   command = {
    #     k = 18
    #     sample_cor <- cor(get_lMs(M = M_cpm, genes = auto_genes),
    #       method = 'pearson')
    #     sample_clust <- hclust(as.dist(1 - sample_cor))
    #     sample_group <- cutree(sample_clust, k = k)
    #     return(sample_group)
    #   }
    # ),

    # tar_target(
    #   name = gene_group,
    #   command = {
    #     k = 40
    #     gene_cor <- cor(t(get_lMs(M = M_cpm, genes = auto_genes)),
    #       method = 'pearson')
    #     gene_clust <- hclust(as.dist(1 - gene_cor))
    #     gene_group <- cutree(gene_clust, k = k)
    #     return(gene_group)
    #   }
    # ),

    # tar_target(
    #   name = all_lineplots,
    #   command = {

    #     fns <- purrr::map(unique(gene_group), function(g) {
    #       exp_dynamics_panel(
    #         features = names(which(gene_group == g)),
    #         meta_data = sample_annotation,
    #         merge_cn = 'sample_name',
    #         version = glue::glue('murine_GM{g}'),
    #         lookup_data = M_cpm,
    #         redo = F,
    #         leave_out_sn = T
    #       )
    #     })

    #     return(fns[[length(fns)]])
    #   },
    #   format = 'file'
    # ),

    tar_target(
      name = mirjam_lineplots,
      command = {

        fns <- 
          c(mirjam_genes_220516, mirjam_genes_220518) %>%
          purrr::imap(function(genes, name) {
          exp_dynamics_panel(
            features = genes,
            meta_data = sample_annotation,
            merge_cn = 'sample_name',
            version = name,
            lookup_data = M_cpm,
            redo = F,
            leave_out_sn = T
          )
        })

        return(fns[[length(fns)]])
      },
      format = 'file'
    ),

    # tar_target(
    #   name = fano_only_lineplots,
    #   command = {
    #     fano_genes <- rownames(f_gene_var)
    #     fns <- exp_dynamics_panel(
    #       features = setdiff(fano_genes, auto_genes),
    #       meta_data = sample_annotation,
    #       merge_cn = 'sample_name',
    #       version = glue::glue('murine_fano_only_genes'),
    #       lookup_data = M_cpm,
    #       redo = F,
    #       leave_out_sn = T
    #     )
    #     return(fns[[length(fns)]])
    #   },
    #   format = 'file'
    # ),

    # tar_target(
    #   name = limma_only_lineplots,
    #   command = {
    #     genes <- limma_genes$gene
    #     fns <- exp_dynamics_panel(
    #       features = setdiff(genes, auto_genes),
    #       meta_data = sample_annotation,
    #       merge_cn = 'sample_name',
    #       version = glue::glue('murine_limma_only_genes'),
    #       lookup_data = M_cpm,
    #       redo = F,
    #       leave_out_sn = T
    #     )
    #     return(fns[[length(fns)]])
    #   },
    #   format = 'file'
    # ),

    NULL
  )
)
