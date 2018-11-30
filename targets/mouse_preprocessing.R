mouse_so_preprocessing_targets <- list(
  tar_target(
    name = bulk_5892_so,
    command = {
      M <- tar_read(M_mouse)
      sample_annotation <- tar_read(sample_annotation_mouse)
      stopifnot(tolower(sample_annotation$sample_name) == colnames(M))
      colnames(M) <- sample_annotation$sample_name
      # colnames(M$counts)
      exp5892 <- CreateSeuratObject(
        # counts = log2(M$counts + 1) %>%
        counts = M$counts %>%
          { set_colnames(., tolower(colnames(.))) },
        meta.data = sample_annotation %>%
          as.data.frame() %>%
          set_rownames(NULL) %>%
          dplyr::mutate(sample_name = tolower(sample_name)) %>%
          tibble::column_to_rownames('sample_name')
      )
      exp5892 <- SCTransform(exp5892)
      exp5892 <- FindVariableFeatures(exp5892)
      exp5892 <- RunPCA(exp5892, assay = 'SCT',
        npcs = 20, verbose = F)
      return(exp5892)
    }
  ),

  tar_target(
    name = outlier_markers_mouse,
    command = {
      # Compare the right and left islands. Turns out the right island is
      # lower in almost ALL genes
      so <- tar_read(filtered_so_6743)
      
      dbscan_cl <- 
        so@reductions$umap@cell.embeddings %>%
        fpc::dbscan(eps = 2, MinPts = 20)
      
      so[['dbscan']] <- dbscan_cl$cluster

      # ## 8-fold difference in UMIs
      so@meta.data[, c('nCount_RNA', 'dbscan_cluster')] %>%
        dplyr::group_by(dbscan_cluster) %>%
        dplyr::summarize(
          'median' = median(nCount_RNA), 
          'mean' = mean(nCount_RNA))

      ## Cluster 2 (dying cells) more prevalent in the T-cell
      ## condition
      so@meta.data[, c('stim_group', 'dbscan_cluster')] %>%
        tally(stim_group, dbscan_cluster) %>%
        dplyr::group_by(stim_group) %>%
        dplyr::summarize(across(), f_2 = freq[2] / (freq[2] + freq[1]))
      
      ## Pseudobulk based on mouse and DBSCAN cluster
      so_pb <- 
        Seurat::AggregateExpression(so, 
          group.by = c('dbscan', 'condition_i'), 
        return.seurat = TRUE)
      # colnames(so_pb@assays$SCT[,])
      
      markers <- FindMarkers(so_pb, 
        ident.1 = rownames(so_pb@meta.data)[1:18], 
        ident.2 = rownames(so_pb@meta.data)[19:36], 
        test.use = 'roc'
      )
      
      markers <- markers %>%
        dplyr::arrange(desc(myAUC), desc(avg_log2FC))
      
      markers <- markers %>%
        tibble::rownames_to_column('gene')

      return(markers)
    }
  ),

  tar_target(
    name = plot_outlier_markers_mouse,
    command = {
      markers <- tar_read(outlier_markers_mouse)
      out_dir <- file.path(Sys.getenv('img_dir'), 'mouse_TRANSACT')
      of <- file.path(out_dir,
        glue::glue('characterize_island.pdf'))

      so[['logUMI']] <- log(so[['nCount_RNA']])
      p1 <- FeaturePlot(so, features = 'logUMI') +
        theme_cyto_inf() +
        gg_tabula_rasa

      p2 <- markers %>%
        ggplot(aes(x = avg_log2FC, y = myAUC, label = gene)) + 
        geom_point(alpha = .4) +
        ggrepel::geom_label_repel(data = markers[1:2, ], 
          min.segment.length = 0) +
        ylab('ROC-value of gene')
      
      p3 <- FeaturePlot(so, combine = T, features = 'Taco1')[[1]] +
        theme_cyto_inf() +
        gg_tabula_rasa
      p4 <- FeaturePlot(so, combine = T, features = 'Rbp1')[[1]] +
        theme_cyto_inf() +
        gg_tabula_rasa
      
      comp <- (p1 + p2) / (p3 + p4)
      print_plot_eval(print(comp),
        width = 17.4, height = 15,
        filename = of)

      return(of)
    },
    format = 'file'
  ),

  tar_target(
    name = outlier_markers_mouse_informative_only,
    command = {
      ## When looking at hand-picked 'informative' genes only, the
      ## islands no longer appear. Apparently, the divisive genes are
      ## among the most informative ones with respect to the observed
      ## cytokines (at least the ones we're interested in)
      so <- tar_read(filtered_so_6743)

      so <- so[unname(unlist(all_mouse_genesets)), ]
      VariableFeatures(so) <- unname(unlist(all_mouse_genesets))
      so <- process_so(so, redo_variable_genes = F)

      so[['logUMI']] <- log(so[['nCount_RNA']])
      p1 <- FeaturePlot(so, features = 'logUMI') +
        theme_cyto_inf() +
        gg_tabula_rasa

      p1 <- FeaturePlot(so, features = 'logUMI') +
        theme_cyto_inf() +
        gg_tabula_rasa

      p2 <- DimPlot(
        so, 
        group.by = 'stim_group', 
        label.size = 2,
        label = FALSE,
        raster = TRUE) +
        theme_cyto_inf() +
        ggtitle('') +
        guides(color = guide_legend(ncol = 2)) +
        gg_tabula_rasa

      of <- file.path(out_dir,
        glue::glue('characterize_island_mirjam_genes.pdf'))
      print_plot_eval(print(p1 + p2),
        width = 10, height = 12,
        filename = of)

      return(of)
    },
    format = 'file'
  ),

  ## As can be seen in the result of `plot_outlier_markers_mouse`,
  ## a large clump of cells has considerably lower UMI cound (~50x).
  ## They are filtered out here and the data is reprocessed
  tar_target(
    name = filtered_cleaned_so_6743_manual,
    command = {
      so <- tar_read(filtered_so_6743)
      
      dbscan_cl <- 
        so@reductions$umap@cell.embeddings %>%
        fpc::dbscan(eps = 2, MinPts = 20)
      # table(dbscan_cl$cluster)
      
      so[['dbscan']] <- dbscan_cl$cluster
      # table(so[['dbscan']]) %>% { . / sum(.) }
      so_f <- so[, so@meta.data$dbscan == 1]
      so_f <- process_so(so_f)
    }
  ),

  tar_target(
    name = sce_mouse,
    command = {
      limma_MRs <- tar_read(stringent_limma_MRs_mouse)
      so <- tar_read(filtered_cleaned_so_6743_manual)
      sce <- perform_agg_neighbourhoods(
        so = so,
        k = 10,
        d = 20,
        genes = limma_MRs,
        genelist = NULL
      )
      sce <- runUMAP(sce, dimred = 'PCA', name = 'umap')
      colData(sce)$mouse <- colData(sce)$condition_i <- 
        colData(sce)$dominant_hashtag
      return(sce)
    }
  ),

  tar_target(
    name = sce_mouse_var_genes,
    command = {
      so <- tar_read(filtered_cleaned_so_6743_manual)
      
      sce <- perform_agg_neighbourhoods(
        so = so,
        k = 10,
        d = 20,
        genes = VariableFeatures(so),
        genelist = NULL
      )
      sce <- runUMAP(sce, dimred = 'PCA', name = 'umap')
      colData(sce)$mouse <- colData(sce)$condition_i <- 
        colData(sce)$dominant_hashtag
      return(sce)
    }
  )
)
