human_targets <- list(
  tar_map(
    names = type,
    values = tibble(type = 'human'),

    tar_target(
      name = efit,
      command = {
        fit_limma()
      },
      iteration = 'list'
    ),

    tar_target(
      name = tM,
      command = {
        tM <- efit$t[read_geneset('informativeV15'), ]
        tM[, c(paste0('stim_group', c('10 ng/ml TNFa',
              '100 ng/ml IFNy', '100 ng/ml IFNy 10 ng/ml TNFa')))]
      }
    ),


    tar_target(
      name = gene_var,
      command = {
        M_cpm <- tar_read(kallisto_5029)
        sa <- tar_read(sample_annotation_exp5029)
        stopifnot(all(sa$sample_name == colnames(M_cpm)))
        library(scran)
        gene_var <- modelGeneVar(log2(M_cpm + 1),
          block = t(sa[, 'duration', drop = F])
        )
        gene_var <- as.data.frame(gene_var)
        return(gene_var)
      }
    ),

    tar_target(
      name = f_gene_var,
      command = {
        desirable <- gene_var$bio > 0.25 & gene_var$FDR <= 1e-9
        message(sum(desirable))
        gene_var[desirable, ]
      }
    ),

    tar_target(
      name = example_fano_lineplots,
      command = {
        M_cpm <- as.matrix(tar_read(kallisto_5029))
        sa <- tar_read(sample_annotation_exp5029)
        fns <- exp_dynamics_panel(
          features = rownames(f_gene_var) %>%
            { sample(., 4*8L-1) },
          meta_data = sa,
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
      name = limma_MR_params,
      command = list(LL = 4.2, UL = 1, SD = .5, fd = 2)
    ),

    tar_target(
      name = limma_MRs,
      command = {
        tM <- tar_read(tM_human)
        params <- tar_read(limma_MR_params_human)
        RC <- tM >= params$LL
        NRC <- tM < params$UL

        MRs <- names(which(
          apply(RC[, 1:2], 1, sum) == 1 &
          apply(NRC, 1, sum) >= 1 &
          abs(apply(tM[, 1:2], 1, sum) - tM[, 3]) <= SD
        ))
        tM[MRs, ]

        synergy_genes <- 
          names(which(tM[, 3] > params$fd * apply(abs(tM[, 1:2]), 1, sum)))
        tM[synergy_genes, ]
        return(MRs, synergy_genes)
      }
    ),

    tar_target(
      name = MR_clusters,
      command = {
        # tM <- tar_read(tM_human)
        # params <- tar_read(limma_MR_params_human)
        M_cpm <- as.matrix(tar_read(kallisto_5029))
        sa <- tar_read(sample_annotation_exp5029)
        lM <- subset_feats(as.matrix(M_cpm), limma_MRs)
        k = 40
        gene_cor <- cor(t(lMs), method = 'pearson')
        gene_clust <- hclust(as.dist(1 - gene_cor))
        gene_group <- cutree(gene_clust, k = k)
      }
    ),

    # tar_target(
    #   name = gene_group,
    #   command = {
    #     k = 40
    #     gene_cor <- cor(t(lMs), method = 'pearson')
    #     gene_clust <- hclust(as.dist(1 - gene_cor))
    #     gene_group <- cutree(gene_clust, k = k)
    #     return(gene_group)
    #   }
    # ),
     
    tar_target(
      name = monoreporter_lineplots,
      command = {
        M_cpm <- as.matrix(tar_read(kallisto_5029))
        sa <- tar_read(sample_annotation_exp5029)
        MRs <- 
          ## rmd/select_genes_6.Rmd
          readRDS(file.path(rds_dir, 'human_MR_gene_group.rds')) %>%
          names()
        fns <- exp_dynamics_panel(
          features = MRs,
          meta_data = sa,
          merge_cn = 'sample_name',
          version = 'human_mono_reporters',
          lookup_data = M_cpm,
          redo = F,
          leave_out_sn = T
        )
        return(fns[[length(fns)]])
      },
      format = 'file'
    ),

    NULL
  )
)
