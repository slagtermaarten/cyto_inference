plotNhoodGraphDA <- function (x, milo_res, alpha = 0.05,
  res_column = 'logFC', p_column = 'SpatialFDR', ...) {
  library(SummarizedExperiment)
  if (!miloR:::.valid_graph(nhoodGraph(x))) {
    stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
  }
  if (is.character(layout)) {
    if (!layout %in% names(reducedDims(x))) {
      stop(layout, "is not in readucedDim(x) - choose a different layout")
    }
  }
  signif_res <- milo_res
  signif_res[signif_res[[p_column]] > alpha, res_column] <- 0

  colData(x)[res_column] <- NA
  ## Original line
  # colData(x)[unlist(nhoodIndex(x)[signif_res$Nhood]), res_column] <- signif
  ## Replacement
  colData(x)[unlist(nhoodIndex(x)), res_column] <-
    signif_res[match(unlist(nhoodIndex(x)), signif_res$Nhood), res_column]

  plotNhoodGraph(x, colour_by = res_column, ...)
}


# plotNhoodGraph <- function (
#   x, layout = "UMAP", colour_by = NA, subset.nhoods = NULL,
#     size_range = c(0.5, 3), node_stroke = 0.3, ...) {

#     if (!.valid_graph(nhoodGraph(x))) {
#         stop("Not a valid Milo object - neighbourhood graph is missing. Please r
# un buildNhoodGraph() first.")
#     }
#     if (is.character(layout)) {
#         if (!layout %in% names(reducedDims(x))) {
#             stop(layout, "isn't in readucedDim(x) - choose a different layout")
#         }
#     }
#     nh_graph <- nhoodGraph(x)
#     if (!is.null(subset.nhoods)) {
#         nh_graph <- igraph::induced_subgraph(nh_graph, vids = which(as.numeric(V
# (nh_graph)$name) %in%
#             unlist(nhoodIndex(x)[subset.nhoods])))
#     }
#     nh_graph <- permute(nh_graph, order(vertex_attr(nh_graph)$size,
#         decreasing = TRUE))
#     if (is.character(layout)) {
#         redDim <- layout
#         layout <- reducedDim(x, redDim)[as.numeric(vertex_attr(nh_graph)$name),
#             ]
#         if (!any(class(layout) %in% c("matrix"))) {
#             warning("Coercing layout to matrix format")
#             layout <- as(layout, "matrix")
#         }
#     }
#     if (!is.na(colour_by)) {
#         if (colour_by %in% colnames(SummarizedExperiment::colData(x))) {
#             col_vals <- SummarizedExperiment::colData(x)[as.numeric(vertex_attr(nh_graph)$name),
#                 colour_by]
#             if (!is.numeric(col_vals)) {
#                 col_vals <- as.character(col_vals)
#             }
#             V(nh_graph)$colour_by <- col_vals
#         }
#         else {
#             stop(colour_by, "is not a column in
#             SummarizedExperiment::colData(x)")
#         }
#     }
#     else {
#         V(nh_graph)$colour_by <- V(nh_graph)$size
#         colour_by <- "Nhood size"
#     }
#     if (colour_by %in% c("logFC")) {
#         plot.g <- simplify(nh_graph)
#         pl <- ggraph(simplify(nh_graph), layout = layout) + geom_edge_link0(aes(
# width = weight),
#             edge_colour = "grey66", edge_alpha = 0.2) + geom_node_point(aes(fill
#  = colour_by,
#             size = size), shape = 21, stroke = node_stroke) +
#             scale_size(range = size_range, name = "Nhood size") +
#             scale_edge_width(range = c(0.2, 3), name = "overlap size") +
#             theme_classic(base_size = 14) + theme(axis.line = element_blank(),
#             axis.text = element_blank(), axis.ticks = element_blank(),
#             axis.title = element_blank())
#     }
#     else {
#         pl <- ggraph(simplify(nh_graph), layout = layout) + geom_edge_link0(aes(
# width = weight),
#             edge_colour = "grey66", edge_alpha = 0.2) + geom_node_point(aes(fill
#  = colour_by,
#             size = size), shape = 21, stroke = node_stroke) +
#             scale_size(range = size_range, name = "Nhood size") +
#             scale_edge_width(range = c(0.2, 3), name = "overlap size") +
#             theme_classic(base_size = 14) + theme(axis.line = element_blank(),
#             axis.text = element_blank(), axis.ticks = element_blank(),
#             axis.title = element_blank())
#     }
#     if (is.numeric(V(nh_graph)$colour_by)) {
#         pl <- pl + scale_fill_gradient2(name = colour_by)
#     }
#     else {
#         mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(V
# (nh_graph)$colour_by)))
#         pl <- pl + scale_fill_manual(values = mycolors, name = colour_by,
#             na.value = "white")
#     }
#     pl
# }


run_milo <- function(...) UseMethod('run_milo')


run_milo.Seurat <- function(
  so = NULL, stim_group_only = F) {

  run_milo(so_sce, bi = bi, stim_group_only = stim_group_only)
}


#' 
#'
#' @param sce Milo-SCE object as generated using
#' perform_agg_neighbourhoods
test_milo_DA_Nhoods <- function(
  sce = NULL, verbose = T, debug_f = F) {

  library(miloR)
  library(SingleCellExperiment)
  library(SummarizedExperiment)

  experiment <- SummarizedExperiment::colData(sce)$exp[1] %||%
    SummarizedExperiment::colData(sce)$experiment[1]
  SummarizedExperiment::colData(sce) <- 
    SummarizedExperiment::colData(sce) %>%
    as.data.frame() %>%
    order_duration() %>%
    DataFrame()

  if (!'condition_i' %in% colnames(SummarizedExperiment::colData(sce)) &&
      'dominant_hashtag' %in% colnames(SummarizedExperiment::colData(sce))) {
    SummarizedExperiment::colData(sce)$condition_i <-
      SummarizedExperiment::colData(sce)$dominant_hashtag
  }

  if (unique(SummarizedExperiment::colData(sce)$exp) == '6493') {
    SummarizedExperiment::colData(sce)$mouse <-
      colData(sce)$dominant_hashtag %>%
      as.integer() %>%
      factor()
  }

  if ('mouse' %in% colnames(colData(sce)) &&
      is.character(colData(sce)$mouse)) {
    colData(sce)$mouse <- colData(sce)$mouse %>%
      stringr::str_replace('\\d{4}_', '') %>%
      as.integer() %>%
      factor()
  }

  if (T && !'mouse' %in% colnames(colData(sce))) {
    # table(colData(sce)$condition_name)
    colData(sce)$mouse <- colData(sce)$condition_i
    # colData(sce)$mouse <- NULL
  }

  design_M <-
    data.frame(colData(sce)) %>%
    # tibble::rownames_to_column('sample_name') %>%
    dplyr::select(
      condition_i, any_of(c('duration', 'stim_group'))
    ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(condition_i) %>%
    set_rownames(NULL) %>%
    { . }

  if ('duration' %in% colnames(design_M)) {
    design_M <- dplyr::mutate(
      design_M, 
      duration = factor(duration)
    )
  }

  form <- 
    paste(colnames(design_M)[2:ncol(design_M)], collapse = ' + ') %>%
    { paste0('~', .) } %>%
    { as.formula(.) }

  design_M <- 
    model.matrix(
      form
      , data = design_M
    ) %>%
    as.data.frame()
  colnames(design_M) <- colnames(design_M) %>%
    stringr::str_replace_all(' |\\+|-|\\/', '_')

  form <- 
    paste(colnames(design_M)[2:ncol(design_M)], collapse = ' + ') %>%
    { paste0('~', .) } %>%
    { as.formula(.) }

  all_NH_test_DA <-
    colnames(design_M)[2:ncol(design_M)] %>%
    auto_name() %>%
    purrr::map(function(cn) {
      miloR::testNhoods(
        sce,
        model.contrasts = cn,
        design = form,
        design.df = as.data.frame(design_M)
      ) %>%
      dplyr::mutate(Nhood = unlist(nhoodIndex(sce))) %>%
      { set_rownames(., .$Nhood) } %>%
      # dplyr::arrange(desc(logFC)) %>%
      { . }
    })

  return(all_NH_test_DA)
}


plot_DA_Nhoods <- function(
  sce, DA_Nhoods, 
  fn_app='', 
  out_dir = Sys.getenv('img_dir'),
  test_mode = NA) {

  library(SummarizedExperiment)
  library(miloR)
  library(SingleCellExperiment)

  exp_string <- gen_exp_string(colData(sce)$exp)
  o_fn <- file.path(out_dir,
    glue::glue('milo_embedding-exp{exp_string}{fn_app}\\
      {make_flag(test_mode)}.pdf'))

  # source('~/MirjamHoekstra/R/init.R')

  # p <-
  #   sce %>%
  #   # SingleCellExperiment() %>%
  #   tidy() %>%
  #   as_tibble() %>%
  #   dplyr::mutate(mouse = factor(mouse)) %>%
  #   plot_dim_reduc(
  #     colour_var = 'mouse',
  #     print_to_file = F,
  #     annotation_code =
  #       rlang::expr(p <- p + facet_wrap(~condition_name, ncol = 4))) +
  #     guides(colour = guide_legend(title = 'Mouse'))

  # plotNhoodGraph
  # rm(plotNhoodGraphDA)

  p_dat <-
    as.data.frame(cbind(colData(sce), reducedDim(sce, 'UMAP'))) %>%
    order_stim_group() %>%
    order_duration() %>% 
    { . }

  UMAPs <- map(c('duration', 'stim_group'), function(cv) {
    plot_dim_reduc(
      p_dat,
      coord_regex = 'UMAP_',
      colour_var = cv,
      use_stim_group_cols = FALSE,
      shape_var = NULL
    ) +
    theme(legend.position = 'bottom') +
    guides(
      color = guide_legend(
        ncol = 1,
        override.aes = list(alpha = 1)
      )
    ) +
    theme()
  })

  plots <-
    c(
      UMAPs,
      list(
        plotNhoodSizeHist(sce) + 
          theme_cyto_inf(
            legend.position = 'right', 
            legend.direction = 'vertical'
          )
      ),
      purrr::imap(DA_Nhoods, function(.x, .y) {
        plotNhoodGraphDA(sce, .x, alpha = 0.05) +
          ggtitle(.y) +
          theme_cyto_inf() +
          gg_tabula_rasa +
          xlab('') + ylab('') +
          guides(size = 'none', edge_width = 'none')
      })
    ) %>%
    maartenutils::plot_panel_layout(
      labels = NULL,
      plot_direct = test_rendering(),
      nrow = 3, ncol = 2,
      filename = o_fn)

  return(o_fn)
}


perform_agg_neighbourhoods <- function(
  so,
  assay = 'SCT',
  datatype = 'data',
  k = 20,
  d = 10,
  prop = .1,
  genes = NULL,
  genelist = 'informativeV15',
  avg_assay = 'logcounts') {

  library(SingleCellExperiment)
  library(scater)
  # devtools::install_github("MarioniLab/miloR", ref="devel")
  library(miloR)
  # BiocManager::install("tidySingleCellExperiment")
  # library(tidySingleCellExperiment)

  so <- order_duration(so)
  so <- order_condition_name(so)
  so <- separate_duration(so)
  # so <- numerify_regressors(so)

  if (k >= floor(dim(so)[2] / 1)) 
    return(NULL)

  ## Compute PCA on the features of interest
  if (!is.null(genelist)) {
    genes <- union(genes, read_geneset(genelist))
  } else {
    stopifnot(!is.null(genes))
  }
  VariableFeatures(so) <- genes
  so <- RunPCA(so, npcs = d)
  sce <- as.SingleCellExperiment(so)

  # sce <- subset_feats(so,
  #   genes = genes, genelist = genelist) %>%
  #   as.SingleCellExperiment()

  ## logcounts contains the 'data' field from the Seurat object
  ## counts contains the 'counts' field from the Seurat object
  # M1 <- logcounts(sce)
  # M1a <- subset_feats(GetAssayData(so, assay = 'SCT', 'data'),
  #   genes = read_geneset(genelist))
  # M1 - M1a
  # M1 <- counts(sce)
  # M1a <- subset_feats(GetAssayData(so, assay = 'SCT', 'counts'),
  #   genes = read_geneset(genelist))
  # M1 - M1a
  # print_plot_eval(
  #   {
  #     plot(as.vector(counts(sce)) ~ as.vector(logcounts(sce)))
  #   },
  #   width = 17.4, height = 15,
  #   filename = file.path(Sys.getenv('img_dir'),
  #     glue::glue('test.pdf')))
  ## Try to figure out how counts and logcounts relate to each other
  if (F) {
    tibble(
      y = as.vector(logcounts(sce)),
      x = as.vector(counts(sce))
    ) %>%
    dplyr::mutate(x = x + runif(length(x), max = 1e-12)) %>%
    { nls(y ~ log(x+b, a), data = ., start = list(a = 2, b = 1),
      trace = T) }
    logcounts(sce) - log(as.vector(counts(sce)) + 1.091, 3.910)
  }

  experiment <- colData(sce)$exp[1]

  # colData(sce)$condition_name

  sce <- Milo(sce)
  sce <- buildGraph(sce, k = k, d = d)
  # sce <- buildFromAdjacency(sce, k = k, d = d)
  # assay(sce, 'logcounts')
  sce <- makeNhoods(
    sce,
    prop = prop, k = k, d = d, refined = TRUE
  )
  sce <- countCells(
    sce,
    meta.data = data.frame(colData(sce)),
    samples = c('condition_name', 'duration')[1]
  )
  sce <- calcNhoodExpression(sce, assay = avg_assay)
  sce <- buildNhoodGraph(sce)
  metadata(sce) <- list(
    experiment = experiment,
    k = k, d = d, prop = prop,
    genelist = genelist
  )
  sce <- calcNhoodDistance(
    sce, d = d, 
    reduced.dim = 'pca.corrected'
  )

  sce <- runUMAP(sce, dimred = 'PCA', name = 'umap')

  return(sce)
}


#'  Create a Seurat object of Neighbhorhoods, rather than cells
#'
#'
gen_milo_so <- function(
  so, 
  sce, 
  min_neighbourhood_size = 0, 
  ## Useless argument, remove before next rerun of results
  condition_levels = levels(so@meta.data$condition_name),
  query = so@meta.data$exp[1]) {


  library(miloR)
  library(SingleCellExperiment)
  ## Ensure the NHC matrix will be sorted by condition_i correctly
  colData(sce)$condition_i <- 
    factor(colData(sce)$condition_i,
      levels = sort(unique(colData(sce)$condition_i))
    )

  # colData(sce)$mouse <- colData(sce)$condition_i <- 
  #   colData(sce)$dominant_hashtag

  if (is.null(so@meta.data$condition_i)) {
    idxs <- match(rownames(so@meta.data), rownames(colData(sce)))
    so@meta.data$condition_i <- colData(sce)[idxs, 'condition_i']
    # tibble(colData(sce)[, 'condition_i']
    # so@meta.data$condition_i <- as.integer(so@meta.data$condition_name)
  }

  sce <- 
    countCells(
      sce,
      meta.data = data.frame(colData(sce)),
      samples = c('condition_i')
      # samples = c('stim_group')
    )
  NHC <- miloR::nhoodCounts(sce)
  NHC_s <- apply(NHC, 1, sum)
  # NHC_n <- NHC / NHC_s

  allowed_Nhoods <- which(NHC_s >= min_neighbourhood_size)
  NHC <- as.data.frame(NHC[allowed_Nhoods, ])
  NHC_s <- NHC_s[allowed_Nhoods]

  # sce <- calcNhoodExpression(sce)
  NE <- miloR::nhoodExpression(sce)[, allowed_Nhoods]

  ## SC sample annotation
  so <- order_condition_name(so)
  sc_sa <- extract_sa(so,
    meta_fields = c('condition_i', 'stim_group', 'duration')) %>%
    dplyr::distinct() %>%
    # order_condition_name() %>%
    order_duration() %>%
    dplyr::arrange(condition_i) %>%
    { . }
  # sc_sa <- 
  #   dplyr::arrange_at(sc_sa, intersect(c('condition_i'), 
  #       colnames(sc_sa)))
  
  sa <- as.data.frame(NHC)
  colnames(sa) <- paste0('CN', 1:ncol(sa))
  # cidx <- match(condition_levels, sc_sa$condition_name)
  # sa <- sa[, cidx]
  # sa <- sa[cidx, ]
  # NHC <- NHC[, cidx]

  # colnames(NHC) <- paste0('CN', 
  #   match(colnames(NHC), condition_levels))
  sa$N <- NHC_s
  rownames(sa) <- colnames(NE)
  sa$experiment <- sa$exp <- query

  q_so <- SeuratObject::CreateSeuratObject(NE, meta.data = sa)

  ## Copy over RNA slot to SCT slot for compatibility with the
  ## bulk/reference object
  M <- as.matrix(GetAssayData(q_so, 'counts', assay = 'RNA')) 
  q_so[['SCT']] <- CreateAssayObject(counts = M)
  q_so <- FindVariableFeatures(q_so)

  return(q_so)
}


gen_NH_CN_HM <- function(dtf = NULL, sce = NULL, side = 'top', 
  meta_fields = c('duration', 'stim_group'),
  allowed_Nhoods = NULL, sa_ann = NULL, 
  cluster_columns = F,
  forego_sa = FALSE, ...) {

  library(miloR)
  library(SummarizedExperiment)
  library(SingleCellExperiment)

  if (!forego_sa && is.null(sa_ann)) {
    sa <-
      SummarizedExperiment::colData(sce) %>%
      { .[, c(meta_fields, 'condition_i'), drop = F] } %>%
      as.data.frame() %>%
      dplyr::distinct() %>%
      # dplyr::arrange(condition_i) %>%
      # dplyr::select(-condition_i) %>%
      # dplyr::arrange(stim_group, )
      { . }

    o_vars <- c('stim_group', 'duration', 
          'tnfa_conc', 'ifny_conc', 'sn_dilution') %>%
      intersect(colnames(sa))

    sa <- sa %>%
      order_duration() %>%
      dplyr::arrange_at(o_vars, .keep_all = T)

    if ('condition_name' %in% colnames(sa)) {
      sa <- sa %>%
        order_condition_name() %>%
        dplyr::arrange(condition_name) %>%
        dplyr::select(any_of(c(meta_fields, 'condition_name'))) %>%
        { . }
    }

    # sa_ann <- dplyr::select(sa, any_of(o_vars))
  }

  sa_ann <- sa %>%
    dplyr::select(-any_of('condition_i'))
  # M <- dtf %>%
  #   dplyr::select(matches('CN\\d+')) %>%
  #   as.matrix() %>%
  #   t() %>%
  #   # set_rownames(CN_levs) %>%
  #   # { .[levels(sa$condition_name), ] } %>%
  #   { . }
  library(SummarizedExperiment)
  M <- as.matrix(nhoodCounts(sce))
  stopifnot(!is.null(sa$condition_i))
  M <- M[, naturalsort::naturalsort(colnames(M))]
  M <- M[, sa$condition_i]
  rownames(M) <- unlist(nhoodIndex(sce))

  if (!is.null(allowed_Nhoods)) {
    # dtf <- dtf[match(allowed_Nhoods, dtf$rn), ]
    # dtf <- dtf[match(allowed_Nhoods, dtf$rn), ]
    M <- M[match(allowed_Nhoods, unlist(nhoodIndex(sce))), ]
  }
  M <- t(M)
  dn <- dimnames(M)
  M <- M %*% diag(1/apply(M, 2, sum))
  dimnames(M) <- dn
  # rownames(M) <- levels(sa$condition_name)

  if (side == 'left') {
    M <- t(M)
    ra <- NULL
    ca <- sa_ann
  } else if (side == 'top') {
    ra <- sa_ann
    ca <- NULL
  }

  if (forego_sa)
    ra <- ca <- NULL

  # HS <- c(17.4, 20)
  HM <- gen_HM(
    M,
    # width = unit(HS[1], 'cm'),
    # height = unit(nrow(M)*.5, 'cm'),
    # ra = sa_ann,
    ra = ra,
    ca = ca,
    # cluster_columns = get_obj('NH_expression_clustering')[['reference']],
    cluster_columns = cluster_columns,
    # cluster_rows = get_obj('NH_expression_clustering')[['query']],
    cluster_rows = T,
    # show_column_names = F,
    name = 'Neighbourhood abundance\nof condition',
    ...
  )

  return(HM)
}


#' Visualize gene expression in neighbourhoods
#'
#' Plots the average gene expression in neighbourhoods, sorted by DA fold-change
#'
#' @param x A \code{\linkS4class{Milo}} object
#' @param da.res a data.frame of DA testing results
#' @param features a character vector of features to plot (they must be in rownames(x))
#' @param alpha significance level for Spatial FDR (default: 0.1)
#' @param subset.nhoods A logical, integer or character vector indicating a subset of nhoods to show in plot
#' (default: NULL, no subsetting)
#' @param cluster_features logical indicating whether features should be clustered with hierarchical clustering.
#' If FALSE then the order in \code{features} is maintained (default: FALSE)
#' @param assay A character scalar that describes the assay slot to use for calculating neighbourhood expression.
#' (default: logcounts)
#' Of note: neighbourhood expression will be computed only if the requested features are not in the \code{nhoodExpression} slot
#' of the milo object. If you wish to plot average neighbourhood expression from a different assay, you should run
#' \code{calcNhoodExpression(x)} with the desired assay.
#' @param scale_to_1 A logical scalar to re-scale gene expression values between 0 and 1 for visualisation.
#' @param show_rownames A logical scalar whether to plot rownames or not. Generally useful to set this to
#' \code{show_rownames=FALSE} when plotting many genes.
#' @param highlight_features A character vector of feature names that should be highlighted on the right side of
#' the heatmap. Generally useful in conjunction to \code{show_rownames=FALSE}, if you are interested in only a few
#' features
#' @return a \code{ggplot} object
#'
#' @author Emma Dann
#'
#' @examples
#' NULL
#'
#' @export
#' @rdname plotNhoodExpressionDA
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr mutate left_join filter percent_rank first group_by summarise
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @importFrom stats hclust
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_replace
plotNhoodExpressionDA <- function(x, 
  da.res, 
  features, alpha=0.1,
  subset.nhoods=NULL, cluster_features=FALSE, assay="logcounts",
  scale_to_1 = FALSE,
  show_rownames=TRUE,
  highlight_neighbourhoods = NULL,
  highlight_features = NULL) {

  if (length(features) <= 0 | is.null(features)) {
    stop("features is empty")
  }
  ## Check if features are in rownames(x)
  if (!all(features %in% rownames(x))) {
    stop("Some features are not in rownames(x)")
  }
  ## Check if nhood expression exists
  if (dim(nhoodExpression(x))[2] == 1){
    warning("Nothing in nhoodExpression(x): computing for requested features...")
    x <- calcNhoodExpression(x, assay = assay, subset.row = features)
  }
  ## Check if all features are in nhoodExpression
  if (!all(features %in% rownames(nhoodExpression(x)))) {
    warning("Not all features in nhoodExpression(x): recomputing for requested features...")
    x <- calcNhoodExpression(x, assay = assay, subset.row = features)
  }

  expr_mat <- nhoodExpression(x)[features, ]
  ## Bug fix?!
  # colnames(expr_mat) <- seq_len(ncol(nhoods(x)))
  colnames(expr_mat) <- unlist(nhoodIndex(x))

  ## Get nhood expression matrix
  if (!is.null(subset.nhoods)) {
      expr_mat <- expr_mat[,subset.nhoods, drop=FALSE]
  }

  if (!isFALSE(scale_to_1)) {
      expr_mat <- t(apply(expr_mat, 1, function(X) (X - min(X))/(max(X)- min(X))))
      # force NAs to 0?
      if(sum(is.na(expr_mat)) > 0){
          warning("NA values found - resetting to 0")
          expr_mat[is.na(expr_mat)] <- 0
      }
  }

  rownames(expr_mat) <- 
    sub(pattern = "-", replacement = ".", rownames(expr_mat)) ## To avoid problems when converting to data.frame

  my_percent_rank <- function(x) {
    (dplyr::min_rank(x))/(sum(!is.na(x)) + 1)
  }
  # my_percent_rank(1:4)
  # my_percent_rank(4:1)
  # browser()
  if (!is.data.frame(da.res) && is.list(da.res)) {
    ## Complete me, allow plotting of multiple DAs
    da.res.l <- da.res
    da.res <- imap_dfr(da.res.l, ~mutate(.x, name = .y))
    # da.res <- da.res.l[[1]]
  } else {
    da.res$name = ''
  }

  pl_df <- data.frame(t(expr_mat)) %>%
    tibble::rownames_to_column('Nhood') %>%
    mutate(Nhood = as.double(Nhood)) %>%
    left_join(da.res, by = 'Nhood') %>%
    mutate(logFC_rank = my_percent_rank(logFC))
  eps <- min(pl_df$logFC_rank)/2

  ## Top plot: nhoods ranked by DA log FC
  pl_top <- pl_df %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, paste0('SpatialFDR < ', alpha), NA)) %>%
    ggplot(aes(logFC_rank, logFC)) +
    geom_hline(yintercept = 0, linetype=2) +
    geom_point(size=0.2, color="grey") +
    geom_point(data=.%>% filter(!is.na(is_signif)), aes(color=is_signif), size=1) +
    theme_bw(base_size=8) +
    ylab("DA logFC") +
    scale_color_manual(values="red", name="") +
    scale_x_continuous(expand = c(0.00, 0)) +
    coord_cartesian(xlim = c(0+eps, 1-eps)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())


  if (!is.null(highlight_neighbourhoods)) {
    pl_df_hl <- pl_df %>%
      dplyr::arrange(logFC_rank) %>%
      dplyr::mutate(d = c(0, diff(logFC_rank))) %>%
      dplyr::filter(Nhood %in% highlight_neighbourhoods)
    pl_top <- pl_top + geom_rect(data = pl_df_hl, 
      mapping = aes(xmin = logFC_rank - d * .5, xmax = logFC_rank + d * .5,
        ymin = -Inf, ymax = Inf), 
      alpha = .3, colour = 'grey50')
  }

  ## Bottom plot: gene expression heatmap
  if (isTRUE(cluster_features)) {
      row.order <- hclust(dist(expr_mat))$order # clustering
      # row.order <- as.hclust(gen_clust_object(t(expr_mat)))$order
      ordered_features <- rownames(expr_mat)[row.order]
  } else {
    ordered_features <- rownames(expr_mat)
  }

  # this code assumes that colnames do not begin with numeric values
  # add 'X' to feature names with numeric first characters
  rownames(expr_mat) <- stringr::str_replace(rownames(expr_mat), pattern="(^[0-9]+)", replacement="X\\1")

  pl_df <- pl_df %>%
    tidyr::pivot_longer(cols=rownames(expr_mat), names_to='feature', values_to="avg_expr") %>%
    mutate(feature=factor(feature, levels=ordered_features))

  if (!is.null(highlight_features)) {
    if (!all(highlight_features %in% pl_df$feature)){
      missing <- highlight_features[which(!highlight_features %in% pl_df$feature)]
      warning('Some elements of highlight_features are not in features and will not be highlighted. \nMissing features: ', paste(missing, collapse = ', ') )
    }
    pl_df <- pl_df %>%
      mutate(label=ifelse(feature %in% highlight_features, as.character(feature), NA))
  }

  pl_bottom <- pl_df %>%
    ggplot(aes(logFC_rank, feature, fill=avg_expr)) +
    geom_tile() +
    scale_fill_viridis_c(option="magma", name="Avg.Expr.") +
    xlab("Neighbourhoods") + ylab("Features") +
    scale_x_continuous(expand = c(0.00, 0)) +
    scale_y_discrete(expand = c(0.00, 0)) +
    theme_classic(base_size = 8) +
    coord_cartesian(clip="off") +
    theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          panel.spacing = margin(2, 2, 2, 2, "cm"),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(10,10,10,10)
     )

  if (!is.null(highlight_features)) {
    hl_df <- pl_df %>%
      dplyr::filter(!is.na(label)) %>%
      dplyr::group_by(label) %>%
      summarise(logFC_rank=max(logFC_rank), avg_expr=mean(avg_expr), 
        feature=as.character(feature)[1])
    pl_bottom <- pl_bottom +
      ggrepel::geom_text_repel(data=hl_df,
        aes(label=label, x=logFC_rank),
        size=1,
        xlim = c(max(pl_df$logFC_rank) + 0.01, max(pl_df$logFC_rank) + 0.02),
        force = 10,
        min.segment.length = 0,
        max.overlaps=Inf,
        seed=42)

  }

  if (isFALSE(show_rownames)) {
    pl_bottom <- pl_bottom + theme(axis.text.y=element_blank())
  }

  ## Assemble plot
  (pl_top / pl_bottom) +
    plot_layout(heights = c(1,4), guides = "collect") &
    theme(legend.justification=c(0, 1),
          legend.margin = margin(0,0,0,50))
}
# rm('plotNhoodExpressionDA')

findNhoodGroupMarkers <- function (x, da.res, assay = "logcounts", 
  aggregate.samples = FALSE,
    sample_col = NULL, subset.row = NULL, gene.offset = TRUE,
    subset.nhoods = NULL, subset.groups = NULL, 
    na.function = "na.pass") {
  if (!is(x, "Milo")) {
      stop("Unrecognised input type - must be of class Milo")
  }
  else if (any(!assay %in% assayNames(x))) {
      stop("Unrecognised assay slot: ", assay)
  }
  if (is.null(na.function)) {
      warning("NULL passed to na.function, using na.pass")
      na.func <- get("na.pass")
  }
  else {
      tryCatch({
          na.func <- get(na.function)
      }, warning = function(warn) {
          warning(warn)
      }, error = function(err) {
          stop("NA function ", na.function, " not recognised")
      }, finally = {
      })
  }
  if (isTRUE(aggregate.samples) & is.null(sample_col)) {
      stop("if aggregate.samples is TRUE, the column storing sample information must be spec
ified by setting 'sample_col'")
  }
  if (!"NhoodGroup" %in% colnames(da.res)) {
      stop("'NhoodGroup' columns is missing from da.res. Please run groupNhoods() or define
neighbourhood groupings otherwise.")
  }
  nhs.da.gr <- da.res$NhoodGroup
  names(nhs.da.gr) <- da.res$Nhood
  if (!is.null(subset.nhoods)) {
      nhs.da.gr <- nhs.da.gr[subset.nhoods]
  }
  nhood.gr <- unique(nhs.da.gr)
  nhs <- nhoods(x)
  if (!is.null(subset.nhoods)) {
      nhs <- nhs[, subset.nhoods]
  }
  fake.meta <- data.frame(CellID = colnames(x)[rowSums(nhs) >
      0], Nhood.Group = rep(NA, sum(rowSums(nhs) > 0)))
  rownames(fake.meta) <- fake.meta$CellID
  for (i in seq_along(nhood.gr)) {
      nhood.x <- which(nhs.da.gr == nhood.gr[i])
      nhs <- nhs[rowSums(nhs) > 0, ]
      nhood.gr.cells <- rowSums(nhs[, nhood.x, drop = FALSE]) > 0
      fake.meta[nhood.gr.cells, "Nhood.Group"] <- 
        ifelse(is.na(fake.meta[nhood.gr.cells, "Nhood.Group"]), 
          nhood.gr[i], NA)
  }
  fake.meta <- fake.meta[!is.na(fake.meta$Nhood.Group), ]
  if (!is.null(subset.row)) {
      x <- x[subset.row, , drop = FALSE]
  }
  exprs <- assay(x, assay)[, fake.meta$CellID]
  marker.list <- list()
  i.contrast <- c("TestTest - TestRef")
  if (length(nhood.gr) == 1) {
      if (sum(fake.meta$Nhood.Group == nhood.gr[1]) == nrow(fake.meta)) {
          warning("All graph neighbourhoods are in the same group - cannot perform DGE testi
ng. Returning NULL")
          return(NULL)
      }
  }
  if (isTRUE(aggregate.samples)) {
      fake.meta[, "sample_id"] <- 
        SummarizedExperiment::colData(x)[fake.meta$CellID,
          sample_col]
      fake.meta[, "sample_group"] <- paste(fake.meta[, "sample_id"],
          fake.meta[, "Nhood.Group"], sep = "_")
      sample_gr_mat <- matrix(0, nrow = nrow(fake.meta), ncol = length(unique(fake.meta$sample_group)))
      colnames(sample_gr_mat) <- unique(fake.meta$sample_group)
      rownames(sample_gr_mat) <- rownames(fake.meta)
      for (s in colnames(sample_gr_mat)) {
          sample_gr_mat[which(fake.meta$sample_group == s),
              s] <- 1
      }
      exprs_smp <- matrix(0, nrow = nrow(exprs), ncol = ncol(sample_gr_mat))
      if (assay == "counts") {
          summFunc <- rowSums
      }
      else {
          summFunc <- rowMeans
      }
      for (i in seq_len(ncol(sample_gr_mat))) {
          if (sum(sample_gr_mat[, i]) > 1) {
              exprs_smp[, i] <- summFunc(exprs[, which(sample_gr_mat[,
                i] > 0)])
          }
          else {
              exprs_smp[, i] <- exprs[, which(sample_gr_mat[,
                i] > 0)]
          }
      }
      rownames(exprs_smp) <- rownames(exprs)
      colnames(exprs_smp) <- colnames(sample_gr_mat)
      smp_meta <- unique(fake.meta[, c("sample_group", "Nhood.Group")])
      rownames(smp_meta) <- smp_meta[, "sample_group"]
      fake.meta <- smp_meta
      exprs <- exprs_smp
  }
  if (!is.null(subset.groups)) {
      nhood.gr <- subset.groups
  }
  for (i in seq_along(nhood.gr)) {
      i.meta <- fake.meta
      i.meta$Test <- "Ref"
      i.meta$Test[fake.meta$Nhood.Group == nhood.gr[i]] <- "Test"
      if (ncol(exprs) > 1 & nrow(i.meta) > 1) {
          i.design <- as.formula(" ~ 0 + Test")
          i.model <- model.matrix(i.design, data = i.meta)
          # colSums(i.model)
          rownames(i.model) <- rownames(i.meta)
      }
      if (assay == "logcounts") {
          i.res <- miloR:::.perform_lognormal_dge(exprs, i.model, 
            model.contrasts = i.contrast,
              gene.offset = gene.offset)
      } else if (assay == "counts") {
          i.res <- miloR:::.perform_counts_dge(exprs, i.model, 
            model.contrasts = i.contrast, gene.offset = gene.offset)
          # i.res[which.min(i.res$FDR), ]
          colnames(i.res)[ncol(i.res)] <- "adj.P.Val"
      }
      else {
          warning("Assay type is not counts or logcounts - assuming (log)-normal distribution. Use these results at your peril")
          i.res <- miloR:::.perform_lognormal_dge(exprs, i.model, 
            model.contrasts = i.contrast, gene.offset = gene.offset)
      }
      i.res$adj.P.Val[is.na(i.res$adj.P.Val)] <- 1
      i.res$logFC[is.infinite(i.res$logFC)] <- 0
      i.res <- i.res[, c('logFC', 'adj.P.Val', 'logCPM')]
      colnames(i.res) <- paste(colnames(i.res), nhood.gr[i],
          sep = "_")
      marker.list[[paste0(nhood.gr[i])]] <- i.res
  }
  concord.rownames <- Reduce(x = marker.list, f = function(x,
      y) all(rownames(x) == rownames(y)))
  if (isTRUE(all(concord.rownames))) {
      marker.df <- do.call(cbind.data.frame, marker.list)
      colnames(marker.df) <- gsub(colnames(marker.df), pattern = "^[0-9]+\\.",
          replacement = "")
      marker.df$GeneID <- rownames(marker.df)
  }
  else {
      warning("Rownames of DGE results are reordered")
      n.rows <- lapply(marker.list, nrow)
      if (length(unique(n.rows)) > 1) {
          warning("Not all DGE results contain the same features - results may be truncated"
)
      }
      marker.df <- Reduce(x = marker.list, f = function(x,
          y) merge(x, y, by = 0))
      colnames(marker.df) <- gsub(colnames(marker.df), pattern = "^[0-9]+\\.",
          replacement = "")
      marker.df$GeneID <- rownames(marker.df)
  }
  return(marker.df)
}


#' What fraction of cells is included if we subselect Nhoods?
#'
#'
frac_included_cells <- function(sce, allowed_Nhoods = 1:10,
  allowed_cells = 1:ncol(sce)) {
  if (length(allowed_Nhoods) == 0) return(0)
  tryCatch(
    mean(rowSums(sce@nhoods[allowed_cells, allowed_Nhoods, drop = F]) >= 1), 
    error = function(e) { NA_real_ }) 
}
