`%|||%` <- function(a, b) {
  if (is.null(a) || is.na(a) || (is.vector(a) &&
      length(a) == 0) ||
      (is.data.frame(a) && nrow(a) == 0))
    return(b)
  else
    return(a)
}
# c() %|||% 4 %|||% c()
# 4 %|||% c()

# bm_bulk_experiments <- bulk_experiments <- c('5029', '6434')
# bulk_experiments <- c('4910', '5029', '6434', '6623')
# in_vivo_sc_experiments <- c('5310_in_vivo', '6493', '6601')


gen_exp_string <- function(experiments) {
  sort(unique(as.character(experiments))) %>%
    paste(collapse = '-')
}


exp_plot_dir <- function(experiment) {
  dir_loc <- file.path(Sys.getenv('img_dir'),
    glue::glue('exp{experiment}'))
  dir.create(dir_loc, showWarnings = F)
  return(dir_loc)
}


report_combinations <- function(N, instance = 'KL_div_matrix') {
  N_work <- format(N, big.mark=',')
  mymessage(msg = glue::glue('{N_work} combinations to test'),
            instance = instance)
}


get_cur_date <- function()
  gsub('-|:', '_', format(Sys.time(), '%F_%R'))


if (F && !exists('exp5310')) {
  delayedAssign('exp5310',
                readRDS(file.path(rds_dir, 'exp5310_sc-proc.rds')))
  # exp5310
  if (F) {
    exp5310[['MAGIC']] <-
      retrieve_feature_mat('transact-linear-informative_filtered-bu', 'sc') %>%
      set_colnames(colnames(so)) %>%
      CreateAssayObject()
    DefaultAssay(exp5310) <- 'MAGIC'
  }
}


MAD <- function(x) median(abs(x - median(x)))


#' Set panel size of ggplot grob in gtable
#'
#'
set_panel_size <- function(p=NULL, g=ggplotGrob(p),
  width=unit(5, 'cm'), height=unit(4, 'cm')) {
  panel_index_w<- g$layout$l[g$layout$name=='panel']
  panel_index_h<- g$layout$t[g$layout$name=='panel']
  g$widths[[panel_index_w]] <- width
  g$heights[[panel_index_h]] <- height
  class(g) <- c('fixed', class(g), 'ggplot')
  return(g)
}


query_cytokines <- c('ifn', 'tnf')


axis_labels <- c(
  'duration'='Exposure time [h]',
  'sn_dilution'='Super-natant [dilution]',
  'tnf_conc'='TNFa [ng/ml]',
  'ifn_conc'='IFNy [ng/ml]',
  'tnf_conc_norm'='normalized TNFa score',
  'ifn_conc_norm'='normalized IFNy score',
  'tnf_duration_norm'='normalized TNFa exposure duration',
  'ifn_duration_norm'='normalized IFNy exposure duration',
  'tnf_rank'='TNFa [ng/ml]',
  'ifn_rank'='IFNy [ng/ml]',
  'pc_score'='Factor score',
  'rho' = 'Concentration [ng/ml]',
  'tau' = 'Time [h]',
  'tau_rank' = 'Time [h]',
  'tau_ifn_rank' = 'Time IFNy [h]',
  'tau_tnf_rank' = 'Time TNFa [h]',
  'lik' = 'Likelihood',
  'eta_mu' = 'Gene expression',
  'gene_expression_cpm' = 'CPM',
  'log_gene_expression_cpm' = 'log2(CPM)',
  'gene_expression' = 'TMM-scaled read counts',
  'log_gene_expression' = 'log2(TMM + 1)'
)


cyto_breaks <- list(
  'ifn' = c(0, .01, 1, 100),
  'tnf' = c(0, .1, 1, 10),
  'ifn_rank' = c(0, .01, 1, 100),
  'tnf_rank' = c(0, .1, 1, 10),
  'tau' = c(2, 6, 12, 24),
  'tau_rank' = c(2, 6, 12, 24),
  'tau_tnf_rank' = c(2, 6, 12, 24),
  'tau_ifn_rank' = c(2, 6, 12, 24)
)


cyto_names <- list(
  'ifn' = c(0, .01, 1, 100),
  'tnf' = c(0, .1, 1, 10),
  'ifn_rank' = c(0, .01, 1, 100),
  'tnf_rank' = c(0, .1, 1, 10),
  'tau' = c(2, 6, 12, 24),
  'tau_rank' = c(2, 6, 12, 24)
)


add_flag <- function(fn, flag = '-inf_method=point') {
  base <- gsub('(.*)(\\..+)', '\\1', fn)
  ext <- gsub('(.*)(\\..+)', '\\2', fn)
  mod_fn <- paste0(base, flag, ext)
}


#' Modify an existing filename by adding a flag
#'
#'
add_flag_to_fn <- function(fn, flag = '-inf_method=point') {
  mod_fn <- add_flag(fn = fn, flag = flag)
  if (file.exists(fn) && !file.exists(mod_fn)) {
    file.rename(fn, mod_fn)
  }
  return(mod_fn)
}


gen_summary_string <- function(v, sep = ', ') {
  sum <- summary(v)
  purrr::map_chr(1:6,
    ~glue::glue('{names(sum)[.x]}: {round(sum[.x], 1)}')) %>%
    paste(collapse = sep)
}


load_exp <- function(exp) {
  fn <- file.path(rds_dir, glue::glue('exp{exp}_sc.rds'))
  if (exp != '5092') {
    fn <- fn %>% add_flag('-proc')
  }
  so <- readRDS(fn)
}


object_size_overview <- function() {
  objects <- ls(parent.frame())
  sort(sapply(objects, function(x) { pryr::object_size(get(x)) / 1e6 }))
}
# object_size_overview()


reset_plotting_device <- function() {
  while (length(dev.list()) > 0) {
    dev.off()
  }
}


test_grob <- function(grob, doit = F) {
  if (!doit) return(invisible())
  reset_plotting_device()
  dev.new()
  pushViewport(viewport(width = 0.5, height = 0.5))
  grid.rect()
  draw(grob)
  popViewport()
}


stim_group_levels <- c(
  'Unstimulated in vitro',
  '1/20000 SN',
  '1/200 SN',
  '1/2 SN',
  '0.01 ng/ml IFNy',
  '1 ng/ml IFNy',
  '10 ng/ml IFNy',
  '100 ng/ml IFNy',
  '0.1 ng/ml TNFa',
  '1 ng/ml TNFa',
  '10 ng/ml TNFa',
  '1 ng/ml IFNy 0.1 ng/ml TNFa',
  '1 ng/ml IFNy 1 ng/ml TNFa',
  '1 ng/ml IFNy 10 ng/ml TNFa',
  '10 ng/ml IFNy 0.1 ng/ml TNFa',
  '10 ng/ml IFNy 1 ng/ml TNFa',
  '10 ng/ml IFNy 10 ng/ml TNFa',
  '100 ng/ml IFNy 0.1 ng/ml TNFa',
  '100 ng/ml IFNy 1 ng/ml TNFa',
  '100 ng/ml IFNy 10 ng/ml TNFa',
  'Unexposed in vivo',
  'Exposed to T-cells in vivo',
  'In vivo 10 ng/ml TNFa',
  'In vivo 100 ng/ml IFNy',
  'In vivo 100 ng/ml IFNy 10 ng/ml TNFa',
  '5000 U/ml IFNA1'
)


string_compare <- function(source, target) {
  source_split <- strsplit(source, '')[[1]]
  target_split <- strsplit(target, '')[[1]]
  print(length(source_split) == length(target_split))
  return(all(source_split == target_split))
}


load_genesets <- function() {
  if (!exists('reactome_genesets')) {
    library(genesets)
    ## Make genesets easier on the eye
    reactome_genesets <-
      filter_gmt(pattern = '.*REACTOME.*', gmt_pattern = 'msigdb')
    names(reactome_genesets) <- gsub('REACTOME_', '',
      names(reactome_genesets))
    if (F) {
      names(reactome_genesets) <-
        cyto_inf_cap(names(reactome_genesets),
          cap_first_word_only = T) %>%
        gsub(' Sig', ' sig', .) %>%
        gsub('(T|t)cr', 'TCR', .) %>%
        gsub('Microrna', 'MicroRNA', .) %>%
        gsub('Mrna', 'mRNA', .) %>%
        gsub('^(E|e)r ', 'ER ', .) %>%
        gsub('(M|m)hc', 'MHC', .) %>%
        gsub('b cell', 'B cell', .) %>%
        gsub('class ii', 'class II', .) %>%
        gsub('class i', 'class I', .) %>%
        gsub('Interferon gamma', 'Interferon Gamma', .) %>%
        { . }
    }
    reactome_genesets <<- reactome_genesets
  }

  if (F && !exists('cellmarker_genesets')) {
    cellmarker_genesets <<- suppressWarnings(read_CellMarker('Human'))
  }
}


#' Compile filenames for Seurat objects
#'
#'
compile_so_fns <- function(
  experiment = '5310',
  lookup_mode = F,
  filtering_opts = list(
    min_UMI = 1000, max_percent_mt = 100, min_fd_cc = 2,
    max_hashtag_evenness = 1),
  ...) {
  dots <- list(...)
  filtering_opts <- filtering_opts %>% replace(names(dots), dots)
  fns <- with(filtering_opts,
    tibble(
      vanilla = file.path(rds_dir,
        glue::glue('exp{experiment}_sc.rds')),
      HTO_QC = add_flag(vanilla,
        glue::glue('-HTO_QC{make_flag(min_fd_cc)}\\
          {make_flag(max_hashtag_evenness)}')),
      ## --- SC QC filtered variations
      filtered = add_flag(HTO_QC,
        glue::glue('-filtered{make_flag(max_percent_mt)}\\
          {make_flag(min_UMI)}')),
      filtered_ig = add_flag(filtered, '-ig'),
      filtered_ig_mono = add_flag(filtered, '-ig_mono'),
      filtered_cleaned = add_flag(filtered, '-cleaned'),
      filtered_ig_cleaned = add_flag(filtered_ig, '-cleaned'),
      filtered_ig_mono_cleaned = add_flag(filtered_ig_mono,
        '-cleaned'),
      ## --- Mito reg variations
      filtered_mito_reg =
        add_flag(filtered, '-mito_reg'),
      filtered_ig_mito_reg =
        add_flag(filtered_ig, '-mito_reg'),
      filtered_ig_mono_mito_reg =
        add_flag(filtered_ig_mono, '-mito_reg'),
      filtered_cleaned_mito_reg =
        add_flag(filtered_cleaned, '-mito_reg'),
      filtered_ig_cleaned_mito_reg =
        add_flag(filtered_ig_cleaned, '-mito_reg'),
      filtered_ig_mono_cleaned_mito_reg =
        add_flag(filtered_ig_mono_cleaned, '-mito_reg'),
      ## --- Markers
      marker = add_flag(filtered_cleaned, '-outlying_cluster_markers')
    )
  )

  if (lookup_mode) {
    ## When looking up files, cleaned files that do not exist (because
    ## nothing was to be cleaned from their source files) should be
    ## replaced with their uncleaned counterparts. Symlinks
    ## would have probably been more pretty but alas
    cleaned_fns <- fns %>%
      dplyr::select(matches('cleaned'))

    non_existent_cleaned_files <- cleaned_fns %>%
      purrr::map_lgl(~file.exists(.x)) %>%
      { names(.)[. == F] }

    for (fn in non_existent_cleaned_files) {
      fns[[fn]] <- str_replace(fns[[fn]], '-cleaned', '')
    }

    second_check <- fns[non_existent_cleaned_files] %>%
      purrr::map_lgl(~file.exists(.x)) %>%
      all

    if (!second_check) browser()
    # stopifnot(second_check)
  }

  fns <- fns %>%
    # dplyr::select(-filtered_cleaned, -filtered_ig_cleaned,
    #   -filtered_ig_mono_cleaned) %>%
    as.list() %>%
    { . }

  fns
}
# compile_so_fns(experiment = '5310')
# compile_so_fns(experiment = '5310', min_UMI = 10)


call_with <- function(f,
  args = purrr::map(ls(), ~get(.x, envir = parent.frame()))) {
  if (missing(args) || is.null(args) || !is.list(args))
    stop('No valid args')
  overlapping_args <- intersect(names(formals(f)), names(args))
  # ignored_args <- setdiff(names(formals(f)), names(args))
  ignored_args <- setdiff(names(args), names(formals(f)))
  if (length(ignored_args) > 0) {
    message('Following list items from call to', ' are ignored:',
      paste0(ignored_args, collapse = ', '))
  }
  purrr::exec(f, !!!args[overlapping_args])
}


read_kallisto_counts <- function(fn) {
  bl_columns <- c('chromosome_name', 'gene_biotype', 'start_position',
                  'end_position', 'external_gene_id', 'description')
  dtf <- readr::read_tsv(fn, show_col_types = F) %>%
    as.data.table %>%
    maartenutils::normalize_colnames() %>%
    { set_colnames(., gsub('5029_\\d{1,2}_(.*)_\\w{7}', '\\1',
                           colnames(.))) } %>%
    { set_colnames(., gsub('0_1', '0.1', colnames(.))) } %>%
    { set_colnames(., gsub('0_01', '0.01', colnames(.))) } %>%
    { set_colnames(., gsub('10ng', '10_ng', colnames(.))) } %>%
    cond_setnames('ensg', 'ensembl_gene_id') %>%
    .[, .SD, .SDcols = setdiff(colnames(.), bl_columns)] %>%
    .[!is.na(ensembl_gene_id)]
  dtf
}


tally <- function(x, ...) {
  col_q <- rlang::enquos(...)
  total_n <- nrow(x)
  x %>%
    dplyr::group_by(!!!col_q) %>%
    dplyr::summarise(n = n(), freq = n() / total_n)
}
# rm(tally)
# dplyr:::tally.data.frame
# dplyr:::group_vars
# dplyr:::group_vars.data.frame
# dplyr:::group_data.data.frame
# if (F) {
#   tibble(
#       'letter' = base::sample(LETTERS[1:2], 1000, T),
#       'number' = base::sample(1:2, 1000, T)
#     ) %>%
#     tally(letter, number)
# }

pseudo_bulk <- function(
  so,
  u_var = 'condition_name',
  assay = 'SCT',
  datatype = 'counts',
  TPM = FALSE,
  annotate_clusters = FALSE) {

  u_var <- intersect(u_var, colnames(so@meta.data))
  if (length(u_var) == 0) {
    rlang::warn('No variables to pseudobulk by, returning NULL')
    return(NULL)
  }
  message(
    'Detected the following variables to pseudobulk by: ', 
    paste(u_var, collapse = ', ')
  )

  # if (u_var %nin% colnames(so@meta.data))
  #   stop('u_var not recognized')

  sc_M <- so2M(so, assay = assay, datatype = datatype)
  if (is.null(sc_M)) return(NULL)
  stopifnot(length(so@meta.data[[u_var[1]]]) == ncol(sc_M))

  ## Collapse all 'unique variables/uvars' to a single variable:
  ## 'pb_var'. Order its levels according to the orderings of the
  ## uvars
  so@meta.data$pb_var <- 
    so@meta.data %>%
    dplyr::select(any_of(u_var)) %>%
    apply(1, paste0, collapse = ' - ')

  pb_var_levs <- 
    so@meta.data %>%
    dplyr::distinct(across(all_of(u_var)), .keep_all = T) %>%
    dplyr::arrange_at(u_var) %>%
    dplyr::select(all_of(u_var), pb_var) %>%
    # dplyr::pull(pb_var)
    { . }

  so@meta.data$pb_var <- factor(
    so@meta.data$pb_var, 
    levels = pb_var_levs$pb_var
  )

  ## The rows of sc_M will be ordered according to the levels of
  ## so@meta.data$pb_var
  sc_M <- apply(sc_M, 1, 
    function(x) tapply(x, so@meta.data$pb_var, mean))
  ## Rows are samples here, remove all NA ones
  all_NA_samples <- which(apply(sc_M, 1, function(x) all(is.na(x))))
  if (length(all_NA_samples) > 0) {
    sc_M <- sc_M[-all_NA_samples, ]
    pb_var_levs <- pb_var_levs[-all_NA_samples, ]
  }

  if (TPM) {
    ## Make lib size 1e4, not strictly necessary when TPM is done
    ## afterwards. Rows will features again (and columns samples),
    ## like they ought to be
    sc_M <- apply(sc_M, 1, function(x) 1e4 / sum(x) * x)
    stopifnot(maartenutils::eps(colSums(sc_M),  1e4))
  } else {
    sc_M <- t(sc_M)
  }

  # sum_var <- setdiff(c('condition_name', 'seurat_clusters'), u_var)
  # if (length(sum_var) > 1) browser()
  # u_var_annotation <-
  #   so@meta.data %>%
  #   dplyr::mutate(u_var = .data[[u_var]]) %>%
  #   dplyr::mutate(sum_var = .data[[sum_var]]) %>%
  #   tally(u_var, sum_var) %>%
  #   dplyr::group_by(u_var) %>%
  #   dplyr::mutate(freq = freq / sum(freq)) %>%
  #   { . }

  if (annotate_clusters &&
      'seurat_clusters' %in% colnames(so@meta.data)) {
    ## How do clusters and condition names interrelate?
    cluster_annotation <-
      so@meta.data %>%
      tally(condition_name, seurat_clusters) %>%
      dplyr::group_by(condition_name) %>%
      dplyr::mutate(freq = freq / sum(freq)) %>%
      dplyr::select(-n) %>%
      { . }

    if (u_var == 'condition_name') {
      cluster_annotation_w <-
        tidyr::pivot_wider(cluster_annotation,
          names_from = seurat_clusters,
          values_from = freq, values_fill = 0)
      clean_names <- c()
    } else if (u_var == 'seurat_clusters') {
      # stim_group_ranks <-
      #   dplyr::select(so@meta.data, condition_name, matches('rank')) %>%
      #   unique()
      stim_group_ranks <-
        dplyr::select(so@meta.data, condition_name,
          matches('_conc|duration')) %>%
      unique()

    # cluster_rank_annotation_certainty <- cluster_annotation %>%
    #   group_by(seurat_clusters) %>%
    #   dplyr::mutate(certainty = freq / sum(freq, na.rm = T))

    cluster_annotation %>%
      dplyr::filter(seurat_clusters == 1)

    cluster_rank_annotation <- cluster_annotation %>%
      dplyr::right_join(stim_group_ranks, by = 'condition_name') %>%
      numerify_ranks() %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::summarize(
        # across(matches('rank'),
        across(matches('duration|_conc$'),
          list('med' = function(x) {
            # if (length(ind) > 1) browser()
            if (T) {
              ind <- which(!is.na(x))
              x[ind][which.max(freq[ind])]
            } else{
              ## Weighted sum
              out <- sum(x[ind] * freq[ind], na.rm = T) /
                sum(freq[ind], na.rm = T)
              return(out)
            }
        },
        'certainty' = function(x) {
          ind <- which(!is.na(x))
          return(freq[ind][which.max(freq[ind])] /
            sum(freq, na.rm = T))
        }))
      )

      cluster_annotation_w <-
        tidyr::pivot_wider(cluster_annotation,
          names_from = condition_name,
          values_from = c(freq),
          values_fill = 0) %>%
      dplyr::left_join(cluster_rank_annotation,
        by = 'seurat_clusters')

      clean_names <- c('stim_group', 'condition_name')
    }

    cluster_annotation_w <-
      dplyr::rename_with(cluster_annotation_w,
        ~gsub(' |-', '_', .x))

    # apply(so2M(sc_so), 2, function(x) !all(is.na(x)))
    # so@meta.data %>% dplyr::select(any_of(u_var))
    # so@meta.data[[u_var]]
  }

  new_meta <-
    so@meta.data %>%
    # dplyr::arrange_at(u_var) %>%
    dplyr::distinct(across(any_of(u_var)), .keep_all = T) %>%
    dplyr::select(
      # pb_var,
      -matches('HTO'),
      -matches('dbscan|hashtag|cluster|snn|SCT'),
      -matches('ribo'),
      -matches('HT'),
      -matches('bm_exp'),
      -matches('nCount|nFeature|fd_cc|hash'),
      -matches('Score|Phase|N_UMI|percent.mt')
    ) %>%
    { dplyr::inner_join(pb_var_levs, ., by = colnames(pb_var_levs)) } %>%
    dplyr::select(-pb_var) %>%
    { . }

  if (length(all_NA_samples) > 0) {
    new_meta <- new_meta[-all_NA_samples, ]
  }

  if ('seurat_clusters' %in% colnames(so@meta.data) &&
      exists('cluster_annotation_w')) {
    new_meta <-
      new_meta %>%
      dplyr::left_join(cluster_annotation_w) %>%
      { . }
  }

  if (!(all(u_var %in% c('condition_name', 'mouse', 'tnf_conc',
          'ifn_conc', 'duration')))) {
    new_meta <- dplyr::select(new_meta,
      -matches('rank|conc|dilution'),
      -matches('duration'))
  }

  new_meta <- new_meta[, purrr::map_lgl(new_meta, ~!all(is.na(.x)))]
  # new_meta <- dplyr::select(new_meta, -any_of(clean_names))
  # new_meta <- bind_cols(new_meta, cluster_annotation_w)
  # new_meta <- dplyr::right_join(new_meta, cluster_annotation_w,
  #   by = u_var)

  u_count <-
    so@meta.data %>%
    dplyr::group_by_at(u_var) %>%
    dplyr::summarize(n = n())
  new_meta <- dplyr::right_join(new_meta, u_count, by = u_var)

  new_meta <- set_rownames(new_meta, pb_var_levs$pb_var)
  new_meta <- as.data.frame(new_meta)

  sc <- CreateSeuratObject(
    counts = sc_M, assay = assay,
    meta.data = new_meta
  )

  sc <- sc %>%
    order_duration() %>%
    order_stim_group() %>%
    order_condition_name() %>%
    { . }

  return(sc)
}


find_shared_genes <- function(object_list, ...) {
  object_list <- c(object_list, list(...))
  object_list <- purrr::discard(object_list, 
    # function(x) is.null(x) || all(is.na(x)))
    function(x) is.null(x))

  if (length(object_list) == 0) return(NULL)

  feature_list <- purrr::map(object_list, function(o) {
    if (is.character(o)) {
      return(o)
    } else if (is.matrix(o) || inherits(o, 'dgCMatrix') ||
               inherits(o, 'Seurat')) {
      return(detected_genes(o))
    } else {
      stop('Unexpected data type')
    }
  })

  purrr::reduce(feature_list, intersect)
}
if (F) {
  if (!exists('ol')) {
    ol <- tar_read(filtered_cleaned_so, branch = 1:2)
  }
  feats_1 <- find_shared_genes(object_list = ol)
  feats_2 <- find_shared_genes(object_list = NULL, ol[[1]], ol[[2]])
  all(feats_1 == feats_2)
}


detach_package <- function(pkg, character.only = FALSE) {
  if(!character.only) {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search()) {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


unload_cytokineinference <- function() {
  detach_package(cytokineinference)
  DLL_ov <- getLoadedDLLs()
  for (fp in purrr::map_chr(
      DLL_ov[names(DLL_ov) == 'cytokineinference'],
      ~.x[['path']])) {
    dyn.unload(fp)
  }
}
# unload_cytokineinference()


gen_random_M <- function(r = 10000L, c = 10L, s = 2) {
  out <- matrix(rnorm(r*c) * s, ncol = c, nrow = r, byrow = T)
  rownames(out) <- purrr::map_chr(1:r,
    ~paste0(base::sample(LETTERS, 10L, TRUE), collapse = ''))
  colnames(out) <- purrr::map_chr(1:c,
    ~paste0(base::sample(LETTERS, 10L, TRUE), collapse = ''))
  abs(out)
}


factor_to_rank <- function(v, exp_levels = unique(v)) {
  if (is.character(v)) {
    v <- as.numeric(v)
  }
  if (is.factor(v)) {
    v_n <- levels(v)
    numeric_levels <- suppressWarnings(sort(as.numeric(levels(v))))
    v <- factor(v, levels = c(numeric_levels, setdiff(v_n, numeric_levels)))
    NA_idx <- which(is.na(levels(v)))
    v_i <- as.integer(v)
    v_i[v == NA_idx] <- NA
  } else if (is.numeric(v)) {
    v_i <- frank(v, ties.method = 'dense')
  }
  return(v_i)
}
# v <- structure(c(1L, 3L, 2L, 1L, 4L, 2L, 1L, 1L, 1L),
#                .Label = c("0", "10", "0.1", "1", "Unknown", NA),
#                class = "factor")
# factor_to_rank(v)
# v <- c(100, 0, 10, 0, 10)
# factor_to_rank(v)


#' For moments of mental fog, where is the first difference between
#' two character vectors that are supposedly identical?
#'
#'
find_character_diff <- function(a, b, error_overlap = 4) {
  a_s <- strsplit(a, '')[[1]]
  b_s <- strsplit(b, '')[[1]]
  la <- length(a_s)
  lb <- length(b_s)
  min_l <- min(la, lb)
  idx <- which(a_s[1:min_l] != b_s[1:min_l])
  if (length(idx) == 0) {
    message('No diff')
    return(NULL)
  } else {
    out <- idx[1] %>% {
      start_i <- max(.-error_overlap, 1)
      end_i <- min(.+error_overlap, min_l)
      list(
        paste0(a_s[start_i:end_i], collapse = ''),
        paste0(b_s[start_i:end_i], collapse = '')
      )
    }
    return(out)
  }
}


print_traceback <- function(depth = 3) {
  if (length(.Traceback) == 0) return(NULL)
  idxs <- rev(1:(min(depth, length(.Traceback))))
  # print(.Traceback[idxs])
  # print(rlang::last_trace())
  # print(rlang::last_error())
}
options(error = function(x) traceback(3, max.lines = 2))


#' Identify the dimension for which the names correspons with the
#' expected feats the best
#'
#'
match_dimnames <- function(M, feats) {
  match_count <- purrr::map_int(dimnames(M),
    ~length(intersect(.x, feats)))
  which.max(match_count)
}


min_max_scaling <- function(x) {
  (x - min(x, na.rm = T)) / diff(range(x, na.rm = T))
}


Z_scaling <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
}


cleanup_vctrs_bug <- function(dtf) {
  test_accessibility <- function(x) tryCatch(any(!is.na(x)),
    error = function(e) { FALSE })
  dtf %>%
    dplyr::select(where(test_accessibility))
}


M2tibble <- function(M) {
  M %>%
    {
      rn = rownames(.);
      dplyr::mutate(as_tibble(.), gene = rn)
    } %>%
    { . }
}


#' Replace outliers with maximum within range values
#'
#'
outlier_norm <- function(M, axis = 1) {
  out <- apply(M, axis, function(x) {
    Q <- quantile(x, probs = c(.25, .75))
    thresh <- Q[2] + 1.5 * diff(Q)
    idxs <- which(x > thresh)
    message(sum(length(idxs)))
    x[idxs] <- thresh
    x
  })
  if (axis == 1) {
    out <- t(out)
  }
  return(out)
}


add_MAGIC <- function(
  so,
  assay = DefaultAssay(so),
  datatype = 'counts',
  genelist = 'informativeV15',
  min_avg_exp = 3) {
  library(Seurat)

  M <- so2M(so, assay = assay, datatype = datatype)
  M <- subset_feats(M, read_geneset(genelist))

  if (!is.null(min_avg_exp) &&
      is.numeric(min_avg_exp) &&
      min_avg_exp > 0) {
    gene_stats <- compute_gene_stats(so, genes = rownames(M))

    thresh <- quantile(gene_stats$mean,
      probs = min_avg_exp/100)
    detected_genes <-
      gene_stats %>%
      dplyr::filter(mean >= thresh) %>%
      dplyr::pull(gene)
    M <- subset_feats(M, detected_genes)
  }

  if (any(dim(M) == 0)) return(NULL)

  M_magic <- tryCatch(
    t(as.matrix(Rmagic::magic(t(M), n.jobs = 1L)$result)),
    error = function(e) { print(e); NULL })

  # dimnames(M)
  # dimnames(M_magic)

  if (!is.null(M_magic)) {
    so[['MAGIC']] <- CreateAssayObject(M_magic)
    DefaultAssay(so) <- 'MAGIC'
  }

  return(so)
}


my_render_report <- function(source_fn, params = NULL,
  object_list = list) {
  if (!is.null(params)) {
    params <-
      params[!sapply(params, function(x) is.null(x) || is.na(x))]
    flags <- purrr::imap_chr(params, ~glue::glue('-{.y}={.x}'))
    comb_flags <- paste(flags, collapse = '')
  } else {
    comb_flags <- ''
  }
  doc_id <- str_replace(basename(source_fn), '\\..*$', '')
  o_fn <- file.path(reports_dir,
    glue::glue('{doc_id}{comb_flags}.pdf'))

  expect_params <- rmarkdown::yaml_front_matter(source_fn)$params %>%
    names %>% setdiff(NA)
  rmarkdown::render(
    input = source_fn,
    knit_root_dir = p_root,
    # params = params[expect_params],
    output_file = o_fn,
    envir = base::list2env(object_list)
  )
  return(o_fn)
}
if (F) {
  my_render_report(
    source_fn = file.path(rmd_dir, 'test.rmd'),
    params = list(a = 2, b = 3))
}


#' Seurat object to matrix
#'
#'
## Seurat SCT
## counts -> (corrected) counts,
## data -> log1p(counts),
## scale.data -> pearson residuals;
## sctransform::vst intermediate results are saved in misc
## slot of the new assay.
so2M <- function(so, assay = DefaultAssay(so), datatype = 'counts') {
  if (is.null(so)) return(NULL)
  tryCatch(as.matrix(slot(so@assays[[assay]], datatype)),
    error = function(e) { NULL })
}


se2tibble <- function(obj, append_meta = T) {
  if (is.null(obj)) return(NULL)
  if (!inherits(obj, 'SummarizedExperiment') &&
      !inherits(obj, 'SingleCellExperiment')) return(obj)

  library(SummarizedExperiment)
  library(SingleCellExperiment)
  if (length(assays(obj)) > 1) {
    out <- dplyr::bind_cols(
      tibble::as_tibble(SummarizedExperiment::colData(obj)),
      purrr::exec(bind_cols, purrr::map(
          assays(obj),
          ~tibble::as_tibble(t(.x))))
    )
  } else {
    out <- dplyr::bind_cols(
      tibble::as_tibble(SummarizedExperiment::colData(obj)),
      as_tibble(t(assays(obj)[[1]]))
    )
  }
  if (append_meta) {
    for (f in names(metadata(obj))) {
      out[[f]] <- paste(metadata(obj)[[f]], collapse = '_')
    }
  }
  return(out)
}


drop_levels_everywhere <- function(dtf) {
  factor_vars <-
    names(which(purrr::map_lgl(dtf, ~(class(.x) == 'factor')[1])))
  for (fac in factor_vars) {
    dtf[[fac]] <- droplevels(dtf[[fac]])
  }
  return(dtf)
}


subselect_ann <- function(
  dtf, merge_id = attr(dtf, 'merge_id')) {
  # stopifnot(merge_id %in% colnames(dtf))
  # stopifnot(!is.null(merge_id) && !is.na(merge_id))
  if (any(class(dtf) == 'data.frame')) {
    out <- dtf %>%
      dplyr::select(
        any_of(c(merge_id, 'sample_name', 'stim_group',
            'condition_name')),
        matches('conc|dilution|duration')
      )
    out <- out[, purrr::map_lgl(out, ~!all(is.na(.x)))]
    out <- drop_levels_everywhere(out)
    return(out)
  } else if (any(class(dtf) == 'SummarizedExperiment')) {
    return(as_tibble(SummarizedExperiment::colData(dtf)))
  }
}


recover_ann <- function(
  query,
  lookup_data,
  merge_id = attr(lookup_data, 'merge_id')) {

  stopifnot(is.character(query))

  ann_dat <- subselect_ann(lookup_data)
  tibble::tibble({{merge_id}} := query) %>%
    dplyr::left_join(ann_dat, by = merge_id) %>%
    set_attr('merge_id', merge_id)
}


apply_data_filtering <- function(...)
  UseMethod('apply_data_filtering')


apply_data_filtering.data.frame <-
apply_data_filtering.tbl_df <-
  function(dtf, filtering_settings = NULL, verbose = T) {
  stopifnot(
    is.null(filtering_settings) || is.list(filtering_settings))

  library(rlang)
  if (verbose && !is.null(filtering_settings)) {
    cat(glue::glue('{nrow(dtf)} rows before filtering'), '\n')
  }

  for (fs in filtering_settings) {
    # fs = filtering_settings[[1]]
    ## Not all experiments have these fields, and this is an easy way
    ## of preventing errors
    dtf <- tryCatch(dplyr::filter(dtf, !!fs),
      error = function(e) { dtf })
    if (verbose) {
      cat(glue::glue('{nrow(dtf)} remaining after \\
          {rlang::as_label(fs)}'), '\n')
    }
    if (nrow(dtf) == 0) break
  }

  dtf <- drop_levels_everywhere(dtf)

  return(dtf)
}


apply_data_filtering.Seurat <- function(so,
  filtering_settings = NULL, verbose = T) {

  so@meta.data$tmp_i <- 1:nrow(so@meta.data)
  tmp <- apply_data_filtering(so@meta.data,
    filtering_settings = filtering_settings, verbose = verbose)

  ## Subset the indexes that remain
  return(so[, tmp$tmp_i])
}


#' Convert ENSG to gene symbol and vice versa in (gene) expression
#' matrix
#'
#' @param M Matrix with ensembl gene IDs as rownames
convert_ensg_to_gene_symbol <- function(arg = rna_data, reverse = F,
                                        verbose = F) {
  df_test <-
    any(c('matrix', 'data.frame', 'data.table') %in% class(arg))
  mapping <- fread(
    file.path(raw_10x_support_dir, 'gene_symbol_to_gene.tsv'),
    header = F)
  setnames(mapping, c('ensembl_gene_id', 'external_gene_id'))
  if (reverse) {
    id_col <- 'external_gene_id'
  } else {
    id_col <- 'ensembl_gene_id'
  }
  setkeyv(mapping, id_col)
  value_col <- setdiff(colnames(mapping), id_col)
  if (df_test) {
    if (id_col %in% colnames(arg)) {
      query <- arg[, get(id_col)]
      setDT(arg)
      arg[, (id_col) := NULL]
    } else {
      query <- rownames(arg)
    }
    M_dt <- as.data.table(arg)
  } else if (class(arg) == 'character') {
    query <- arg
  }
  if (mean(query %in% mapping[, get(value_col)]) > .5) {
    N_overlap <- sum(query %in% mapping[, get(value_col)])
    mymessage(
      sprintf('Found %d query genes in target column, indicating that
      conversion has already occurred', N_overlap))
    return(arg)
  }
  found_targets <- mapping[query, , mult = 'all']
  if (any(duplicated(found_targets$external_gene_id))) {
    if (verbose) {
      print(found_targets[, .N, external_gene_id][N > 1])
    }
    found_targets <- mapping[query, , mult = 'first']
  }
  if (df_test) {
    M_dt$gene_id <- found_targets[, get(value_col)]
    M_dt <- M_dt[, lapply(.SD, sum), by = gene_id]
    # M_res <- M_dt[, .SD, .SDcols = setdiff(query, get(id_col))] %>%
    M_res <- M_dt[, !"gene_id"] %>%
      as.matrix %>%
      set_rownames(M_dt$gene_id)
    # tail(M_res, n=1)
    return(M_res)
  } else {
    return(found_targets[, setNames(get(value_col), get(id_col))])
  }
}
# convert_ensg_to_gene_symbol(M)
# convert_ensg_to_gene_symbol(arg = c('ENSG00000238009', 'ENSG00000237491',
#                                     'ENSG00000225880', 'ENSG00000230368'))


extract_var_from_string <- function(char, var) {
  if (grepl(var, char)) {
    res <- gsub(sprintf('^.*%s=(.*)', var), '\\1', char)
    res <- gsub('-.*$', '', res)
    res <- maartenutils::infer_class(res)
  } else {
    res <- NA
  }
  return(res)
}


#' Make a factor object that has levels for values that are not
#' encountered upon initiation
#'
#' @param el Expected labels
include_na_lev <- function(v, el = NULL) {
  v <- as.character(v)
  levs <- unique(c(v, el, 'Unknown', NA))
  labels <- levs
  levs[length(levs)] <- 'Unknown'
  factor(v, levels = levs, labels = labels)
}
# include_na_lev(v = 'in_vitro', el = 'in_vivo')
# include_na_lev(v = 'in_vitro')


#' Reorder factor according to ordering in exp_levels
#'
#'
order_factor <- function(
  v,
  exp_levels = naturalsort::naturalsort(unique(as.character(v)))) {

  missed_levels <-
    unique(as.character(v)) %>%
    setdiff(exp_levels) %>%
    naturalsort::naturalsort()

  factor(v, levels = c(exp_levels, missed_levels))
}


numeric2factor <- function(v,
  exp_levels = unique(as.character(sort(unique(v))))) {
  if (is.factor(v)) return(v)
  factor(as.character(v), levels = exp_levels)
}


#' Reorder factor according to ordering in exp_levels
#'
#'
length_order_factor <- function(
  v,
  exp_levels = sort_by_length(unique(as.character(v)))) {

  missed_levels <-
    unique(as.character(v)) %>%
    setdiff(exp_levels) %>%
    naturalsort::naturalsort()

  factor(v, levels = c(exp_levels, missed_levels))
}


observe_order_factor <- function(v) {
  factor(v, levels = unique(v))
}


numerify_factor <- function(v) {
  if (is.factor(v) || is.character(v)) {
    out <- suppressWarnings(as.numeric(as.character(v)))
    out_I <- as.integer(out)

    if (all(abs(out_I - out) <= 1e-32, na.rm = T)) {
      return(out_I)
    } else {
      return(out)
    }
  } else {
    return(v)
  }
}


fill_in_tnf <- function(dtf) {
  if ('tnf_conc' %in% colnames(dtf)) {
    # Complement missing in_vitro metadata
    tnf_idxs <- !is.na(dtf$ifn_conc) &
      (is.na(dtf$tnf_conc) | all(is.null(dtf$tnf_conc)))

    dtf$tnf_conc[which(tnf_idxs)] <- 0
  }
  return(dtf)
}


numerify_regressors <- function(dtf,
  regressor_vars = c('duration', 'tnf_conc', 'ifn_conc',
    'sn_dilution', 'tnf_duration', 'ifn_duration', 'sn_duration')) {

  if (maartenutils::null_dat(dtf)) return(dtf)
  stopifnot(is.character(regressor_vars))

  dtf <- fill_in_tnf(dtf)

  dtf <- dplyr::mutate(dtf,
    across(any_of(regressor_vars), numerify_factor))

  if ('sn_dilution' %in% colnames(dtf) &&
      !all(is.na(dtf$sn_dilution)) &&
      any(dtf$sn_dilution > 1, na.rm = T)) {
    dtf$sn_dilution <- ifelse(dtf$sn_dilution == 0, 0,
      1 / dtf$sn_dilution)
  }
  return(dtf)
}


add_binary_regressors <- function(dtf,
  regressor_vars = c('duration', 'tnf_conc', 'ifn_conc',
    'sn_dilution', 'tnf_duration', 'ifn_duration', 'sn_duration')) {

  if (maartenutils::null_dat(dtf)) return(dtf)
  stopifnot(is.character(regressor_vars))

  dtf <- fill_in_tnf(dtf)

  for (cn in intersect(regressor_vars, colnames(dtf))) {
    new_cn <- glue::glue('{cn}_bin')
    dtf[[new_cn]] <- as.numeric(dtf[[cn]] > 0)
  }

  return(dtf)
}


reg_ranges <- list(
  'duration' = c(2, 24),
  'tnf_duration' = c(0, 24),
  'ifn_duration' = c(0, 24),
  'sn_duration' = c(0, 24),
  'tnf_conc' = c(0, 10),
  'ifn_conc' = c(0, 100),
  'sn_dilution' = c(0, 1)
)


reg_breaks <- list(
  'duration' = c(2, 6, 12, 24, 48),
  'ifn_conc' = c(0, 1, 10, 100),
  'tnf_conc' = c(0, .1, 1, 10)
)


# tar_read(sample_annotation_exp6434) %>%
#   dplyr::select(matches('conc|sn|rank'))
norm_regressor_method <- 'regular_steps'


norm_regressor <- function(vn, v) {
  if (!vn %in% names(reg_ranges)) stop('Unexpected regressor')

  if (
    is.null(v) || all(is.na(v)) ||
    all(maartenutils::eps(range(v, na.rm = T), c(0, 1))) ||
    all(maartenutils::eps(range(v, na.rm = T), c(0, (48-2)/22)))
  ) {
    return(v)
  }

  if (norm_regressor_method == 'min_max') {
    v <- (v - reg_ranges[[vn]][1]) / diff(reg_ranges[[vn]])
  } else if (norm_regressor_method == 'regular_steps') {
    ## Map unnormalized values to evenly spaced intervals in domain,
    ## rather than unevenly spaced as happens in the transformation
    ## implemented above
    L <- length(reg_breaks[[vn]])
    # new_breaks <- seq(0, 1, length.out = L)
    new_breaks <- seq(1, L, length.out = L)
    v_t <- approx(reg_breaks[[vn]], new_breaks, v)$y
    v_t[is.na(v_t) & v > reg_breaks[[vn]][L]] <- reg_breaks[[vn]][L]
    v_t[is.na(v) & v < reg_breaks[[vn]][1]] <- reg_breaks[[vn]][1]
    v <- v_t
  }
  return(v)
}
norm_regressor('tnf_conc', c(0, 1, 1.5, 10, 100, 90))
# norm_regressor('ifn_conc', c(0, 1, 1.5, 10, 100, 90))


denorm_regressor <- function(vn, v) {
  if (!vn %in% names(reg_ranges)) stop('Unexpected regressor')

  if (norm_regressor_method == 'min_max') {
    v <- v * diff(reg_ranges[[vn]]) + reg_ranges[[vn]][1]
  } else if (norm_regressor_method == 'regular_steps') {
    L <- length(reg_breaks[[vn]])
    new_breaks <- seq(1, L, length.out = L)
    v_t <- approx(new_breaks, reg_breaks[[vn]], v)$y
    v_t[is.na(v_t) & v > reg_breaks[[vn]][L]] <- reg_breaks[[vn]][L]
    v_t[is.na(v) & v < reg_breaks[[vn]][1]] <- reg_breaks[[vn]][1]
    v <- v_t
  }

  v <- factor(v, levels = sort(unique(v)))

  return(v)
}
# denorm_regressor('ifn_conc', c(0, 0.33, 0.3518519, 0.66, 1.00, 0.96))
norm_regressor('ifn_conc', c(0, 1, 1.5, 10, 100, 90)) %>%
  { denorm_regressor('ifn_conc', .) }


norm_regressors <- function(
  dtf,
  regressor_vars = c('duration', 'tnf_duration', 'ifn_duration',
    'sn_duration', 'tnf_conc', 'ifn_conc', 'sn_dilution')) {

  if (maartenutils::null_dat(dtf)) return(dtf)

  for (vn in unique(regressor_vars)) {
    dtf[[vn]] <- norm_regressor(vn, dtf[[vn]])
  }

  if (F) {
    dtf %>%
      dplyr::filter(experiment == '5310_in_vitro') %>%
      dplyr::pull(tnf_duration) %>%
      table
    N_wrong <- dtf %>%
      dplyr::filter(experiment == '5310_in_vitro' &
        tnf_conc > 0 & tnf_conc < 1) %>%
      nrow()
  }

  return(dtf)
}


denorm_regressors <- denumerify_regressors <-
  function(
    dtf,
    regressor_vars = c('duration', 'tnf_conc', 'ifn_conc',
      'sn_dilution', 'tnf_duration', 'ifn_duration', 'sn_duration')) {
  for (vn in intersect(regressor_vars, colnames(dtf))) {
    dtf[[vn]] <- denorm_regressor(vn, dtf[[vn]])
  }

  return(dtf)
}


recover_experiment <- function(dtf) {
  if (is.null(dtf$experiment) &&
      all(stringr::str_detect(dtf$sample_name, '\\d{4}'))) {
    dtf$experiment <- stringr::str_replace(dtf$sample_name,
      '^(\\d{4}).*', '\\1')
  }
  return(dtf)
}


namepair_to_string <- function(v, sep = ' ') {
  purrr::imap_chr(v, ~glue::glue('{.y}={.x}')) %>%
    paste(collapse = sep)
}
namepair_to_string(c('a' = 3, 'b' = 1))


#' Beware! Because of the use of tar_read_raw here, 'targets' does not
#' track dependencies beyond this function. Only use for off-hand
#' exploration/experiments 
#'
#'
read_preproc_experiments <- function(
  experiments, 
  sc_mode = 'SCT') {

  sc_mode <- match.arg(
    sc_mode,
    choices = c('pseudobulk', 'pseudobulk_high_fidelity', 'SCT',
      'magic'),
    several.ok = FALSE
  )

  # shrink_so <- function(so, norm_method = 'CPM') {
  #   so <- lib_normalize(so, norm_method = norm_method)
  #   if (!is.null(genelist))
  #     so <- subset_feats(so, genes = read_geneset(genelist))
  #   if (!is.null(genes))
  #     so <- subset_feats(so, genes = genes)
  #   return(so)
  # }

  if (F) {
    sample_list <- experiments %>%
      intersect(c(sc_e, '5310_in_vitro', '5310_in_vivo')) %>%
      maartenutils::auto_name() %>%
      purrr::map(function(experiment) {
        ei <- experiment2index(e)
        if (sc_mode == 'pseudobulk') {
          so <- tar_read(sc_pseudobulk, branch = ei)[[1]]
        } else if (sc_mode == 'pseudobulk_high_fidelity') {
          so <- tar_read(sc_pseudobulk_high_fidelity, branch = ei)[[1]]
        } else if (sc_mode == 'SCT') {
          so <- tar_read(filtered_cleaned_so, branch = ei)[[1]]
        } else if (sc_mode == 'magic') {
          so <- tar_read(filtered_cleaned_magic_so, branch = ei)[[1]]
        }
        ## If we don't do TMM norm afterwards, we can already subset
        ## genes to save RAM
        # if (norm_method != 'TMM')
        #   so <- shrink_so(so, norm_method = norm_method)
        return(so)
      })
  } else {
    sample_list <- experiments %>%
      intersect(c(sc_e, '5310_in_vitro', '5310_in_vivo')) %>%
      maartenutils::auto_name() %>%
      purrr::map(function(experiment) {
        if (sc_mode == 'pseudobulk') {
          so <- tryCatch(
            targets::tar_read_raw(glue::glue('sc_pseudobulk_{e}')),
            error = function(e) { NULL }
          )
        } else if (sc_mode == 'pseudobulk_high_fidelity') {
          so <- tryCatch(
            targets::tar_read_raw(
              glue::glue('sc_pseudobulk_high_fidelity_{e}')),
            error = function(e) { NULL }
          )
        } else if (sc_mode == 'SCT') {
          so <- tryCatch(
            targets::tar_read_raw(
              glue::glue('filtered_cleaned_so_{e}')),
            error = function(e) { NULL }
          )
        } else if (sc_mode == 'magic') {
          so <- tryCatch(
            targets::tar_read_raw(
              glue::glue('filtered_cleaned_magic_so_{e}')),
            error = function(e) { NULL }
          )
        }
        return(so)
      })
  }

  bulk_exps <-
    intersect(c('4910', '5029', '6434', '6623'), experiments) %>%
    maartenutils::auto_name() %>%
    purrr::map(function(be) {
      so <- tryCatch(
        targets::tar_read_raw(glue::glue('bulk_{be}_so')),
        error = function(e) { print(e); browser() }
      )
      # if (norm_method != 'TMM')
      #   so <- shrink_so(so, norm_method = norm_method)
      return(so)
    })
  if (length(bulk_exps) > 0) {
    sample_list <- c(sample_list, bulk_exps)
    rm(bulk_exps); gc()
  }

  sample_list <- purrr::discard(sample_list, is.null)
  if (length(sample_list) == 0) return(NULL)

  ## Subset to shared features before library size normalizing
  # shared_genes <- find_shared_genes(object_list = sample_list)
  # if (is.null(shared_genes) || length(shared_genes) == 0) return(NULL)
  # sample_list <- purrr::map(sample_list, 
  #   subset_feats, genes = shared_genes)

  return(sample_list)
}


normalize_sample_list <- function(
  sample_list, 
  norm_method = 'none', 
  genelist = NULL,
  genes = NULL, 
  log_counts = T, 
  gene_merge_method = 'intersection',
  data_mod_code = NULL,
  reduce_combine = FALSE,
  split_5310 = T,
  merge_experiments = FALSE,
  filter_gene_reproducibility = NULL,
  GDR_thresh = NULL,
  Z_scale = FALSE) {

  sample_list <- purrr::discard(sample_list, is.null)
  if (length(sample_list) == 0) return(NULL)

  split_5310 <- split_5310 || 
    any(stringr::str_detect(names(sample_list), '_in_'))

  ## Perform library size normalization
  if (norm_method %in% c('CPM', 'logCPM')) {
    sample_list <- purrr::map(sample_list, 
      lib_normalize, norm_method = norm_method)
  } else if (norm_method == 'TMM') {
    sample_list <- 
      purrr::imap(sample_list, 
        ~SeuratObject::AddMetaData(
          .x, rep(.y, ncol(.x)), 
          col.name = 'split_name'
        )
      )

    if (TRUE) {
      sg <- find_shared_genes(sample_list)
      diagnose_sample_list(sample_list)
      sample_list <- map(sample_list, ~subset_feats(.x, genes = sg))
    }

    ## CombineObject/SeuratObject::merge will silently fail if the
    ## objects harbour different, non-overlapping assays, we gotta
    ## prevent that first
    if ('query' %in% names(sample_list)) {
      so <- sample_list[['query']]
      M <- GetAssayData(object = so, slot = 'counts', assay = 'RNA')
      so[['SCT']] <- CreateAssayObject(counts = M)
      DefaultAssay(so) <- 'SCT'
      sample_list[['query']] <- so
      stopifnot(length(unique(purrr::map_chr(sample_list, 
              ~DefaultAssay(.x)))) == 1)
      stopifnot(all(purrr::map_chr(sample_list, 
            ~DefaultAssay(.x)) == 'SCT'))
      diagnose_sample_list(sample_list, assay = 'SCT')
    }

    temp_so <- CombineObject(sample_list, reduce_combine = T)
    temp_so <- RunTMM(temp_so)
    sample_list <- SplitObject(temp_so, 'split_name')
    rm(temp_so); gc()
  } else if (norm_method == 'extensive_TMM') {
    ## Features to use for TMM normalization
    if (T) {
      gs_data <- targets::tar_read(gs_data_step)
      int_feats <-
        gs_data %>%
        # dplyr::filter(Amean >= 3) %>%
        dplyr::filter(geneset %in% c('none', 'blacklist')) %>%
        # intersect(names(which(GDR >= .5))) %>%
        # intersect(rownames(sc_so@assays$SCT[,])) %>%
        dplyr::pull(gene) %>%
        intersect(scn) %>%
        { . }
    } else {
      int_feats <- scn
    }

    if (F) {
      informative_feats = rownames(so2M(sc_so))
      length(int_feats)
      length(informative_feats)
      list(
        'integration' = int_feats,
        'informative' = informative_feats
        ) %>%
        maartenutils::gen_overlap_matrix() %>%
        maartenutils::overlap_analysis(method = 'corroboration')
    }

    source(file.path(r_dir, 'RunTMM.R'))
    ## These parameters barely have an effect on the median ratio
    ## results that are tested below
    TMM_parm_scan <-
      tidyr::expand_grid(
        logRatioTrim = seq(.40, .45, .05),
        sumTrim = c(.05, .05, .01)
      ) %>%
      dplyr::mutate(
        TMM_res = purrr::pmap(., ~with(list(...), {
            multiple_TMM(
              dfs = list(
                so2M(bulk_so, assay = assay, datatype = datatype),
                so2M(scpb_so, assay = assay, datatype = datatype)
              ),
              integration_feats = int_feats,
              logRatioTrim = logRatioTrim,
              sumTrim = sumTrim
            )
          })
        )
      )

    ## Ratio of bulk and sc median library size
    MR <- purrr::map_dbl(TMM_parm_scan$TMM_res, ~with(list(...),
      median(colSums(.[[1]])) / median(colSums(.[[2]]))))
    min_idx <- which.min(abs(MR - 1))
    TMM_norm_df <- TMM_parm_scan[['TMM_res']][[min_idx]]
    bulk_M <- TMM_norm_df[[1]]
    scpb_M <- TMM_norm_df[[2]]
    # print(summary(colSums(sc_M)))
    # print(summary(colSums(bulk_M)))
  } else if (norm_method == 'TPM_normalization') {
    equal_genes = F
    bulk_M <- so2M(bulk_so, assay = assay, datatype = datatype) %>%
      TPM_M(equal_genes = equal_genes)
    scpb_M <- so2M(scpb_so, assay = assay, datatype = datatype) %>%
      TPM_M(equal_genes = equal_genes)
  } else if (norm_method == 'HK_lib_median') {
    gs_data <- targets::tar_read(gs_data_step)
    gs <- gs_data %>%
      dplyr::filter(geneset %in% c('blacklist', 'none')) %>%
      # dplyr::filter(Amax >= 2) %>%
      dplyr::filter(max_t <= 2) %>%
      dplyr::pull(gene) %>%
      find_shared_genes(bulk_so, scpb_so)

    compute_score <- function(so) {
      apply(so2M(so)[gs, ], 2, quantile, .5) %>%
        tibble::enframe('sample_name', 'libsize')
    }

    lib_sizes <- rbind(
        compute_score(bulk_so),
        compute_score(scpb_so)
      ) %>%
      dplyr::mutate(nf = mean(libsize) / libsize)

    N_bulk <- ncol(so2M(bulk_so))
    bulk_M <- so2M(bulk_so) %>%
      { sweep(., 2, lib_sizes$nf[1:N_bulk], FUN = '*') }
    scpb_M <- so2M(scpb_so) %>%
      { sweep(., 2, lib_sizes$nf %>% { .[(N_bulk+1):length(.)] }, 
        FUN = '*') }
  } else if (norm_method == 'quantile') {
    M <- cbind(so2M(bulk_so), so2M(scpb_so))
    M_norm <- preprocessCore::normalize.quantiles(M)
    dimnames(M_norm) <- dimnames(M)
    # apply(M_norm, 2, summary)
    bulk_M <- M_norm[, colnames(so2M(bulk_so))]
    scpb_M <- M_norm[, colnames(so2M(scpb_so))]
  } else if (norm_method == 'batch_scale') {
    compute_nz_median <- function(so) {
      apply(so2M(so), 2, function(tp) {
        tp <- tp[tp > 0]
        quantile(tp, .5)
      })%>%
      tibble::enframe('sample_name', 'libsize')
    }

    lib_sizes <- rbind(
        compute_nz_median(bulk_so),
        compute_nz_median(scpb_so)
      ) %>%
      dplyr::mutate(nf = mean(libsize) / libsize)

    N_bulk <- ncol(so2M(bulk_so))
    bulk_M <- so2M(bulk_so) %>%
      { sweep(., 2, q_lib_sizes$libsize[1:N_bulk], FUN = '*') }
    scpb_M <- so2M(scpb_so) %>%
      { sweep(., 2, q_lib_sizes$nf %>% 
          { .[(N_bulk+1):length(.)] }, FUN = '*') }
  }

  # purrr::map(sample_list, function(so) {
  #   GetAssayData(so, slot = 'counts')[1:5, 1:5]
  # })
  # purrr::map(sample_list, function(so) {
  #   GetAssayData(so, slot = 'data')[1:5, 1:5]
  # })
  # purrr::map(sample_list, function(so) {
  #   colSums(GetAssayData(so, slot = 'data'))
  # })

  if (log_counts) {
    if (stringr::str_detect(norm_method, 'log')) {
      rlang::warn('Doubly logging the data')
    }
    sample_list <- purrr::map(sample_list, function(so) {
      M <- GetAssayData(so, assay = DefaultAssay(so), slot = 'data')
      SetAssayData(
        object = so,
        slot = 'data',
        new.data = log2(as.matrix(M + 1)),
        assay = DefaultAssay(so)
      )
    })
    diagnose_sample_list(sample_list)
    # purrr::map(sample_list, 
    #   ~colSums(so2M(.x, datatype = 'data')))
    # sample_list <- purrr::map(sample_list, 'exp')
  }

  ## Further subset to genes of interest to save RAM
  if (!is.null(genelist))
    sample_list <- purrr::map(sample_list, subset_feats, 
      genes = read_geneset(genelist))
  if (!is.null(genes))
    sample_list <- purrr::map(sample_list, subset_feats, 
      genes = genes)
  gc()

  det_exp5310_exps <- stringr::str_subset(names(sample_list), '5310')
  ## TODO debug this
  if (split_5310 && length(det_exp5310_exps) > 0) {
    # exp5310_exps <- stringr::str_subset(experiments, '5310')
    if (length(exp5310_exps) > 0) {
      es <- sample_list[exp5310_exps] %>%
        purrr::imap(function(so, experiment) {
          out <- Seurat::SplitObject(so, split.by = 'sample_origin') %>%
            { .[[stringr::str_replace(experiment, '5310_', '')]] } 
          if (is.null(out)) return(NULL)
          out %>%
            # rlang::set_names(paste0('5310_', names(.))) %>%
            { .@meta.data$exp <- experiment; . } %>%
            { . }
      })
      # es <- Seurat::SplitObject(sample_list[['5310_in_vitro']],
      #     split.by = 'sample_origin') %>%
      #   rlang::set_names(paste0('5310_', names(.))) %>%
      #   purrr::imap(function(.x, .y) {
      #     .x@meta.data$exp = .y
      #     .x
      #   })
      sample_list <- 
        c(sample_list[!names(sample_list) %in% det_exp5310_exps], es)
      sample_list <- sample_list[experiments]
    }
  }

  if (!is.null(data_mod_code) && !is.na(data_mod_code)) {
    rlang::eval_tidy(data_mod_code)
    # print(map(sample_list, ncol))
  }

  if (gene_merge_method == 'intersection') {
    shared_genes <- find_shared_genes(object_list = sample_list)
  } else if (gene_merge_method == 'union') {
    stop('Union gene_merge_method not implemented yet')
  }

  if (!is.null(genelist)) {
    if (stringr::str_detect(genelist, '^vst')) {
      temp_so <- CombineObject(sample_list)
      shared_genes <- intersect(
        shared_genes,
        load_genes(genelist = genelist, so = temp_so)
      )
      rm(temp_so)
    } else {
      shared_genes <- intersect(
        shared_genes,
        read_geneset(genelist)
      )
    }
  }

  if (!is.null(genes)) {
    shared_genes <- intersect(shared_genes, genes)
  }

  if (!is.null(filter_gene_reproducibility)) {
    shared_genes <- intersect(
      shared_genes,
      filter_gene_reproducibility_score(filter_gene_reproducibility)
    )
  }

  if (length(shared_genes) == 0) return(NULL)

  ## Apply GDR filtering; require gene to be detected in at least x %
  ## of cells in one or more conditions
  detected_sc_experiments <- 
    purrr::map(sample_list, extract_experiment, 
      expect_single = F) %>%
    purrr::discard(~length(.x) != 1) %>%
    intersect(c(sc_e, '5310_in_vitro', '5310_in_vivo')) %>%
    union(c('query')) %>%
    { . }

  if (!is.null(GDR_thresh) && GDR_thresh > 0 && 
      length(detected_sc_experiments) > 0) {

    source('~/MirjamHoekstra/R/init.R')
    # sample_list[['query']]@meta.data
    detected_genes <-
      sample_list[detected_sc_experiments] %>%
      purrr::discard(is.null) %>%
      purrr::map(get_GDR_table, 
        thresh = GDR_thresh, 
        genes = shared_genes, 
        format = 'passing_genes'
      ) %>%
      purrr::discard(is.null) %>%
      rlang::flatten_chr() %>%
      unique()

    shared_genes <- intersect(shared_genes, detected_genes)
  }

  if (length(shared_genes) == 0) return(NULL)

  sample_list <- purrr::map(sample_list, ~.x[shared_genes, ])
  rm(shared_genes)
  gc()

  if (Z_scale) {
    M <- map(sample_list, ~GetAssayData(.x)) %>%
      { purrr::exec(cbind, !!!.) } %>%
      { t(scale(t(.))) } %>%
      { . }
    Ndims_samples <- map_int(sample_list, ~dim(.x)[[2]])
    
    start_idxs <- c(0, Ndims_samples[-length(Ndims_samples)]) + 1
    names(start_idxs) <- names(sample_list)
    end_idxs <- cumsum(Ndims_samples)
    sample_list_mod <- 
      purrr::map(seq_along(sample_list), function(i) {
        SetAssayData(sample_list[[i]], slot = 'data', 
          new.data = M[, start_idxs[i]:end_idxs[i]])
      }) %>%
      rlang::set_names(names(sample_list))
    sample_list <- sample_list_mod
  }

  if (!merge_experiments) {
    return(sample_list)
  } else {
    out <- CombineObject(
      sample_list,
      reduce_combine = reduce_combine
    )
    rm(sample_list); gc()
    return(out)
  }
}


#' Extract the name from a targets::target
#'
#'
extract_target_name <- function(x) {
  if (is.null(x)) return('')
  else return(tryCatch(x$settings$name, error = function(e) { '' }))
}


CI_l <- function(x) {
  x <- x[!is.na(x)]
  mean(x) - 1.96 * sd(x) / sqrt(length(x))
}
CI_h <- function(x) {
  x <- x[!is.na(x)]
  mean(x) + 1.96 * sd(x) / sqrt(length(x))
}


CV <- function(x) { 
  stopifnot(is.vector(x))
  if (length(x) <= 1) 
    return(NA_real_)

  M <- mean(x, na.rm = T)
  SD <- sd(x, na.rm = T)

  if (maartenutils::eps(M, 0, 1e-32)) {
    if (maartenutils::eps(SD, 0, 1e-32)) {
      return(0)
    } else {
      return(NA_real_)
    }
  }

  return(SD / M)
}


# eps <- function(v1, v2 = 0, epsilon = 0.01) {
#   # (!is.na(v1) | !is.na(v2)) & abs(v1 - v2) < epsilon
#   abs(v1 - v2) < epsilon
# }
# eps(v1 = c(NA, .1, .2), v2 = c(NA, NA, .05))

simplify_condition_name <- function(...) UseMethod('simplify_condition_name')


simplify_condition_name.Seurat <- function(
  so, 
  new_cn = 'condition_name',
  rep_regex = NULL) {
  so@meta.data <- simplify_condition_name(
    so@meta.data, 
    new_cn = new_cn,
    rep_regex = rep_regex
  ) 
  return(so)
}


simplify_condition_name.data.frame <- function(
  dtf,
  new_cn = 'condition_name',
  rep_regex = ' - SC digest| - frozen') {

  orn <- rownames(dtf)
  if (!'condition_name' %in% colnames(dtf)) {
    rlang::warn('condition_name not found in columns')
    return(dtf)
  }

  old_levels <- levels(dtf$condition_name)
  new_levels <- stringr::str_replace(old_levels, rep_regex, '')
  
  out <- dtf %>%
    dplyr::mutate(
      !!new_cn := plyr::revalue(
        condition_name, rlang::set_names(new_levels, old_levels))
    ) %>%
    dplyr::mutate(
      !!new_cn := factor(.data[[new_cn]], 
        levels = unique(new_levels))
    ) %>%
    { . }

  set_rownames(out, orn)
}


add_mouse_condition_name <- function(...) UseMethod('add_mouse_condition_name')


add_mouse_condition_name.Seurat <- function(
  so, 
  new_cn = 'condition_name') {
  so@meta.data <- add_mouse_condition_name(
    so@meta.data, 
    new_cn = new_cn) 
  return(so)
}


add_mouse_condition_name.data.frame <- function(
  dtf, new_cn = 'condition_name') {

  orn <- rownames(dtf)
  if (!all(c('condition_name', 'mouse') %in% colnames(dtf))) {
    rlang::warn('not all required columns found')
    return(dtf)
  }

  dtf$mouse <- stringr::str_replace(dtf$mouse, '^\\d{4}_', '')
  dtf$cn_ag <- paste('mouse', dtf$mouse, '-', dtf$condition_name)

  mapping_tab <- dtf %>%
    dplyr::arrange(condition_name) %>%
    dplyr::distinct(condition_name, cn_ag, .keep_all = FALSE) %>%
    { . }
  
  out <- dtf %>%
    # dplyr::mutate(
    #   !!new_cn := plyr::revalue(
    #     condition_name, rlang::set_names(
    #       mapping_tab$cn_ag, 
    #       as.character(mapping_tab$condition_name)))
    # ) %>%
    dplyr::mutate(
      !!new_cn := factor(cn_ag, levels = mapping_tab$cn_ag)
    ) %>%
    dplyr::select(-cn_ag) %>%
    { . }

  set_rownames(out, orn)
}


add_Ag_condition_name <- function(...) UseMethod('add_Ag_condition_name')


add_Ag_condition_name.Seurat <- function(
  so, 
  new_cn = 'condition_name') {
  so@meta.data <- add_Ag_condition_name(
    so@meta.data, 
    new_cn = new_cn) 
  return(so)
}


add_Ag_condition_name.data.frame <- function(
  dtf, new_cn = 'condition_name') {

  orn <- rownames(dtf)
  if (!any(c('condition_name', 'Ag') %in% colnames(dtf))) {
    rlang::warn('condition_name not found in columns')
    return(dtf)
  }

  dtf$cn_ag <- paste(dtf$condition_name, '-', dtf$Ag)

  mapping_tab <- dtf %>%
    dplyr::arrange(condition_name) %>%
    dplyr::distinct(condition_name, cn_ag, .keep_all = FALSE) %>%
    { . }
  
  out <- dtf %>%
    # dplyr::mutate(
    #   !!new_cn := plyr::revalue(
    #     condition_name, rlang::set_names(
    #       mapping_tab$cn_ag, 
    #       as.character(mapping_tab$condition_name)))
    # ) %>%
    dplyr::mutate(
      !!new_cn := factor(cn_ag, levels = mapping_tab$cn_ag)
    ) %>%
    dplyr::select(-cn_ag) %>%
    { . }

  set_rownames(out, orn)
}


test_target_built <- function(target_names) {
  purrr::map_lgl(target_names, function(.x) {
    meta <- targets::tar_meta(.x)
    out <- tryCatch(
      nrow(meta) >= 1 ||
      targets::tar_progress(.x)$progress == 'built',
      error = function(e) { FALSE })

    if (length(out) == 0)
      out <- F

    return(out)
  })
}


tar_read_regex <- function(
  regex,
  finished_only = TRUE,
  min_date = NULL,
  max_date = NULL,
  sample_N = NULL,
  extract_meta_f = NULL,
  mod_f = NULL,
  verbose = TRUE) {

  stopifnot(is.character(regex))
  stopifnot(is.logical(finished_only))

  target_names <- stringr::str_subset(
    targets::tar_meta()$name, regex)

  if (length(target_names) == 0) return(NULL)

  if (finished_only) {
    are_finished <- test_target_built(target_names)
    target_names <- target_names[are_finished]
  }

  if (length(target_names) == 0) return(NULL)

  if (!is.null(min_date)) {
    target_names <- tar_meta() %>%
      dplyr::right_join(tibble(name = target_names)) %>%
      dplyr::select(name, time) %>%
      dplyr::filter(time >= min_date) %>%
      dplyr::pull(name)
  }

  if (!is.null(max_date)) {
    target_names <- tar_meta() %>%
      dplyr::right_join(tibble(name = target_names)) %>%
      dplyr::select(name, time) %>%
      dplyr::filter(time < max_date) %>%
      dplyr::pull(name)
  }

  if (length(target_names) == 0) return(NULL)

  if (verbose) {
    message('Found', length(target_names), 'targets')
  }

  if (!is.null(sample_N)) {
    target_names <-
      sample(target_names, min(length(target_names), sample_N))
  }

  target_names <- naturalsort::naturalsort(target_names)

  if (verbose)
    message(paste0('Found ', length(target_names), ' targets:\n',
        paste(target_names, collapse = '\n')))

  if (is.null(mod_f)) mod_f <- identity
  out <- maartenutils::auto_name(target_names) %>%
    purrr::map(~tryCatch(
        mod_f(tar_read_raw(.x)), error = function(x) NULL)) %>%
    purrr::discard(is.null)

  if (!is.null(extract_meta_f)) {
    extra_meta <- extract_meta_f(target_names)
    out <- purrr::imap(unname(out), ~cbind(.x, extra_meta[.y, ]))
    names(out) <- target_names
    # purrr::map(out, ~out$feature_weights)
  }

  return(out)
}


sort_by_length <- function(v) {
  if (length(v) == 0) return(v)
  stopifnot(is.character(v) || is.factor(v))

  ordering <-
    tibble::tibble(v = v) %>%
    dplyr::mutate(l = stringr::str_length(v)) %>%
    dplyr::arrange(l, v) %>%
    dplyr::pull(v)

  return(ordering)
}
# sort_by_length(c('aa', 'aaa', 'ab', 'b')) == c('b', 'aa', 'ab', 'aaa')


bring_forward <- function(v, a) {
  vctrs::vec_c(a, setdiff(v, a))
}


extract_stimulus <- function(dtf) {
  dtf
}


#' Maintain attributes of a function's argument in the function's
#' output, without having to adapt the function to do
#'
#'
maintain_attributes <- function(f) {
  wrapper <- function(...) {
    fun_args <- list(...)
    atts <- purrr::map(fun_args, ~attributes(.x)) %>%
      rlang::flatten()
    atts <- atts[!names(atts) %in% c('dim')]
    res <- f(...)
    for (a in names(atts)) {
      attr(res, a) <- atts[[a]]
    }
    return(res)
  }
  return(wrapper)
}


most_frequent <- function(x) {
  out <- tibble(v = x) %>%
    dplyr::group_by(v) %>%
    dplyr::summarize(N = n()) %>%
    dplyr::filter(N == max(N, na.rm = T)) %>%
    dplyr::pull(v)
  if (length(out) > 1)
    out <- out[1]
  return(out)
}
# most_frequent(c(3, 4, 4, 4, 5))
# most_frequent(c(3, 4, 4, 4, 5, 5, 5))
# most_frequent(c(3, 4, 4, 4, 5, 5, 5, 5))


rep_zero <- function(x) {
  if (length(x) > 0) return(x)
  else return(0)
}


#' Subset the rows that will be used for benchmarking purposes
#'
#'
subset_bm_query <- function(bm_dtf) {
  if (maartenutils::null_dat(bm_dtf)) {
    return(NULL)
  }

  query_obj <- extract_sc(bm_dtf) %>%
    dplyr::filter(
      group_id == 'sc-in_vitro' |
      experiment %in% c(6600, 6601)
    )

  if ('subsampled' %in% colnames(query_obj)) {
    query_obj <- query_obj[which(query_obj$subsampled == TRUE), ]
  }

  return(query_obj)
}


factor_to_numeric <- function(v) {
  if (!is.factor(v)) return(as.numeric(v))
  as.numeric(levels(v))[v]
}


#' Take a multiline character object and extract the last non-empty
#' line
#'
extract_last_line <- function(v, verbose = T) {
  out <- v %>%
    strsplit('\n') %>%
    purrr::pluck(1) %>%
    stringr::str_replace_all(' ', '') %>%
    stringr::str_replace_all('\t', '') %>%
    { .[. != ''] } %>%
    dplyr::last() %>%
    { . }

  if (verbose)
    message(out)

  return(out)
}

experiment2index <-
e2i <- function(experiment) {
  experiments <- tar_read(e)
  # stopifnot(experiment %in% experiments)
  experiment <- stringr::str_replace(
    experiment, '(\\d{4}).*', '\\1')
  return(which(experiments == experiment))
}


get_filtering_settings <- function(experiment) {
  targets::tar_read(sc2f) %>%
    tibble::deframe() %>%
    purrr::pluck(experiment) %>%
    { tar_read(benchmark_sc_filtering)[[.]] }
}


column_present <- function(dtf, vn) {
  vn %in% colnames(dtf) && !is.null(dtf[[vn]]) &&
    !all(is.na(dtf[[vn]]))
}


column_present_nz <- function(dtf, vn) {
  column_present(dtf = dtf, vn = vn) && any(dtf[[vn]] > 0)
}


read_magic_targets <- function(sample_N = NULL) {
  objects <- tar_read_regex(
    'filtered_cleaned_magic_so_\\d+_.*',
    sample_N = sample_N
  )

  ## Read out meta data from objects
  meta <- purrr::imap_dfr(objects, ~
    tibble(
      name = .y,
      experiment = .x@meta.data$exp[1],
      min_avg_reads =
        as.integer(stringr::str_replace(.y,
            '.*so_(\\d+)_.{8}$', '\\1')),
      oname = glue::glue('{experiment} {min_avg_reads}')
    )
  )

  names(objects) <- meta$oname
  meta$object <- objects

  return(meta)
}


partial_pca <- function(M, N_PCs = 10L, weight_by_var = T) {
  pca.results <- irlba::irlba(M, nv = N_PCs)
  # total.variance <- sum(Seurat:::RowVar(x = t(M)))
  # sdev <- pca.results$d/sqrt(max(1, nrow(x = M) - 1))
  # if (weight_by_var) {
  #   feature.loadings <- pca.results$u %*% diag(pca.results$d)
  # } else{
  #   feature.loadings <- pca.results$u
  # }
  cell.embeddings <- t(pca.results$v)
  colnames(cell.embeddings) <- colnames(M)
  rownames(cell.embeddings) <- paste0('PC_', 1:N_PCs)
  return(cell.embeddings)
}


lib_normalize <- function(so, norm_method = 'CPM') {
  if (is.null(so)) return(NULL)
  if (norm_method == 'TMM') {
    ## Slow for large datasets, not optimal for sparse data such as
    ## scRNASeq
    so_P <- RunTMM(so)
  } else if (norm_method == 'CPM') {
    ## Data is stored in the 'data' field, rather than 'counts' or
    ## 'scale.data'
    # colSums(GetAssayData(so_P, slot = 'data', assay = 'SCT'))
    # colSums(GetAssayData(so, slot = 'data', assay = 'SCT'))
    so_P <- NormalizeData(so,
      normalization.method = 'RC', scale.factor = 1e6
    )
  } else if (norm_method == 'logCPM') {
    so_P <- NormalizeData(so, normalization.method = 'LogNormalize')
  } else if (norm_method == 'none') {
  }
  return(so)
}


#' Directly save object to targets object store
#'
#'
save_target <- function(obj, name) {
  saveRDS(obj, file.path('_targets', 'objects', name))
}


perform_drop_non_discriminatory <- function(dtf) {
  as.data.frame(dtf)[,
    map_lgl(dtf,
      ~!all(is.na(.x)) &&
        (data.table::uniqueN(setdiff(.x, 'Unknown')) > 1 ||
          is.numeric(.x))), drop = F]
}


get_grading_vars <- function(sc_experiment) {
  switch(
    sc_experiment,
    '6369' = c('ifn_duration', 'ifn_conc'),
    '6489' = c('sn_duration', 'sn_dilution'),
    # '5310' = c('tau_rank', 'tnf_rank', 'ifn_rank')
    '5310' = c('ifn_duration', 'ifn_conc', 
      'tnf_duration', 'tnf_conc'),
    NULL
  )
}


scale_so <- function(so, feature_weights) {
  ## Such that all included features will be used in the PCA
  ## embedding
  ## 2022-02-24 10:43 scale/prop all mess up harmony
  ## embeddings. Best is no ScaleData step in combination with SCT
  if (feature_weights == 'scale') {
    so <- ScaleData(so, do.center = T, do.scale = T, 
      assay = DefaultAssay(so))
  } else if (feature_weights == 'prop') {
    so <- ScaleData(so, do.center = T, do.scale = F,
      assay = DefaultAssay(so))
  } else if (feature_weights == 'experiment_spec_prop') {
    sample_list <- SplitObject(so, split.by = 'exp')

    sds <- purrr::map(sample_list, 
      function(so) {
        VariableFeatures(so) <- detected_genes(so)
        so <- ScaleData(so, do.center = T, do.scale = F, 
          assay = DefaultAssay(so))
        # so@assays$TMM
        # so@assays
        GetAssayData(so, slot = 'scale.data')
      })
    sep_scaled <- purrr::reduce(sds, cbind)
    so <- SetAssayData(
      object = so,
      slot = 'scale.data',
      new.data = as.matrix(sep_scaled),
      assay = DefaultAssay(so)
    )
  } else if (feature_weights %in% c('no_scale', 'none')) {
    ## Just copy over data to scale.data
    so <- SetAssayData(
      object = so,
      slot = 'scale.data',
      new.data = GetAssayData(object = so, slot = 'counts') %>%
        as.matrix(),
      assay = DefaultAssay(so)
    )
  }
  gc()
  return(so)
}


read_sample_annotation <- function(experiments = '6369') {
  purrr::map(experiments, function(experiment) {
    experiment <- stringr::str_replace(
      experiment, '(\\d{4}).*$', '\\1')
    if (experiment %in% sc_e) {
      sa <- tar_read_raw(
        glue::glue('sc_{experiment}_sample_annotation')
      )
    } else {
      sa <- tar_read_raw(
        glue::glue('sample_annotation_exp{experiment}')
      )
    }
    return(sa)
  })
}


merge_sa_on_rownames <- function(...) UseMethod('merge_sa_on_rownames')

 
merge_sa_on_rownames.data.frame <- function(
  dtf, experiment, 
  sa = read_sample_annotation(experiment)[[1]],
  i = 0) {

  sa$rn <- rownames(sa)
  if (all(c('mouse', 'stim_group') %in% colnames(sa))) {
    sa$condition_name_mouse <- 
      with(sa, glue::glue('{condition_name} - {mouse}'))
    sa$stim_group_mouse <- 
      with(sa, glue::glue('{stim_group} - {mouse}'))
  }
  colname_count <- 
    colnames(sa) %>%
    setdiff(c('HT1', 'HT2', 'sample_group')) %>%
    maartenutils::auto_name() %>%
    purrr::map_int(~sum(as.character(sa[[.x]]) %in% rownames(dtf))) %>%
    { . }

  if (all(colname_count == 0)) {
    if (i <= 0) {
      rownames(dtf) <- 
        stringr::str_replace(rownames(dtf), '^\\d{4}_', '')
      return(merge_sa_on_rownames(dtf, experiment, i = i + 1))
    }
    rlang::warn('No column found')
    return(dtf)
  }

  merge_cn <- 
    colname_count %>%
    { .[. == max(.)][1] } %>%
    names()
  
  tibble::rownames_to_column(dtf, merge_cn) %>%
    dplyr::left_join(sa, by = merge_cn)
}


merge_sa_on_rownames.matrix <- function(M, sa = NULL, experiment = NULL) {
  dtf <- as.data.frame(M)
  if (!is.null(experiment) && !is.na(experiment)) {
    dtf <-
      dtf %>%
      merge_sa_on_rownames(sa = sa, experiment = experiment) %>%
      dplyr::mutate(exp = experiment, experiment = experiment) %>%
      { . }
  } else {
    dtf <- dtf %>%
      dplyr::mutate(sample_name = rownames(dtf))
  }
  return(dtf)
}


#' Row bind items, harmonizing column types by coercing to the first
#' observed version of each column
#'
#'
harmonize_bind_rows <- function(...) {
  items <- list(...)

  # types <- purrr::imap(items, 
  #   ~tibble::enframe(purrr::map_chr(auto_name(.x), class)))
  # purrr::map(items, colnames)
  # for (i in 1:length(items)) {
  #   items[[i]]$stim_group <- as.character(items[[i]]$stim_group)
  # }

  ## All the R-targets stuff is going to my brain, I'm going mental.
  ## Why does this code generate the following warning!?!
  # Warning message:
  # Problem while computing `stim_group = order_stim_group(stim_group)`.
  # ℹ the condition has length > 1 and only the first element will be used
  items_non_factor <- 
    suppressWarnings(purrr::map(items, function(item) {
      as_tibble(purrr::map(item, function(.y) { 
        if (class(.y) %in% c('glue', 'factor')) {
          as.character(.y)
        } else {
          .y
        }
      }))
    }))

  types <- purrr::imap(items_non_factor, 
    ~purrr::map_chr(.x, ~class(.x)[1])) %>%
    purrr::map_dfr(~tibble::enframe(.x, 'cn', 'type')) %>%
    dplyr::distinct(cn, .keep_all = T)

  items_harmonized <- 
    purrr::map(items, function(.x) {
      as_tibble(purrr::imap(.x, function(.y, .z) { 
        desired_type <- types[['type']][types$cn == .z]
        if (desired_type != 'factor') {
          as(.y, desired_type) 
        } else {
          factor(.y)
        }
      }))
    })

  purrr::exec(dplyr::bind_rows, !!!items_harmonized)
}


first_non_NA <- function(v) {
  v[which(!is.na(v))[1]]
}
# first_non_NA(c(NA, NA, 1, 2))

gen_so_read_count_summary <- function(so) {
  experiment <- Project(so)
  lib_sizes <- colSums(so@assays$RNA)
  tibble::tibble(
    experiment = experiment,
    type = ifelse(experiment %in% bulk_e, 'bulk', 'SC'),
    N_samples = ncol(so),
    N_conditions = length(unique(so@meta.data$condition_name)),
    mean_library_size = mean(lib_sizes, na.rm = T),
    median_library_size = median(lib_sizes, na.rm = T),
    sd_library_size = sd(lib_sizes, na.rm = T),
    cv_library_size = sd_library_size / mean_library_size
  )
}


extract_experiment <- function(...) UseMethod('extract_experiment')


extract_experiment.NULL <- function(x) {
  return(NULL)
}


extract_experiment.default <- function(x) {
  if (is.null(x)) {
    return(NULL)
  } else {
    rlang::abort('Unexpected type')
  }
}


extract_experiment.Seurat <- function(so, expect_single = T) {
  out <- tryCatch(unique(so@meta.data$exp), 
      error = function(e) { NULL })
  if (length(out) > 1 && expect_single) {
    rlang::warn('Found multiple experiments in Seurat object')
  }
  return(out)
}


in2cm <- function(x) x / 2.54


force_numeric <- function(v) {
  if (is.factor(v)) {
    v <- factor_to_numeric(v)
  }
  return(as.numeric(v))
}


norm_columns <- function(M) {
  apply(M, 2, function(x) norm(as.matrix(x), type = 'F'))
}


# pacman::p_load('Rlibstree')
# remotes::install_github("omegahat/Rlibstree")
## Couldn't get Rlibstree to install; don't strictly need it
extract_common_prefix <- function(v, return_substring = TRUE) {
  if (is.null(v)) return(NULL)
  stopifnot(is.character(v))
  L <- purrr::map_int(v, stringr::str_length)
  minL <- min(L)
  ## Strings all equal up to the i-th character?
  i <- 1
  ae <- length(unique(substr(v, 1, i))) == 1
  while (ae && i <= minL) {
    i <- i + 1
    ae <- length(unique(substr(v, 1, i))) == 1
  }
  if (return_substring) {
    if (i == 1) {
      rlang::warn('No common prefix')
    }
    return(substr(v[1], 1, i-1))
  } else {
    return(i)
  }
}
# extract_common_prefix(c('aaaaaZZZ', 'aaaaaDKJDF')) == 'aaaaa'


extract_common_postfix <- function(v, return_substring = TRUE) {
  if (is.null(v)) return(NULL)
  stopifnot(is.character(v))
  L <- purrr::map_int(v, stringr::str_length)
  minL <- min(L)
  ## Strings all equal up to the i-th character?
  i <- 1
  ae <- length(unique(substr(v, L-i, L))) == 1
  while (ae && i <= minL) {
    i <- i + 1
    ae <- length(unique(substr(v, L-i, L))) == 1
  }
  if (return_substring) {
    if (i == 1) {
      rlang::warn('No common postfix')
    }
    return(substr(v[1], L[1]-i+1, L[1]))
  } else {
    return(L-i)
  }
}
# extract_common_postfix(c('ZZaaaaa', 'BBRRaaaaa')) == 'aaaaa'


diagnose_sample_list <- function(sample_list, slot = NULL, assay = NULL) {
  if (interactive()) {
    # map(sample_list, ~colSums(GetAssayData(.x)))
    # map(sample_list, ~dim(GetAssayData(.x)))
    # map(sample_list, ~compute_GDR(GetAssayData(.x)))
    # map(sample_list, ~compute_CDR(GetAssayData(.x)))
    mg <- function(v) cat('--------', v, '--------\n')
    mg('Lib size')
    print(map(sample_list, ~summary(colSums(GetAssayData(.x, assay = assay)))))
    mg('Dim')
    print(map(sample_list, ~dim(GetAssayData(.x, assay = assay))))
    mg('GDR')
    print(map(sample_list, ~summary(compute_GDR(GetAssayData(.x, assay = assay)))))
    mg('CDR')
    print(map(sample_list, ~summary(compute_CDR(GetAssayData(.x, assay = assay)))))
  }
}

q1 <- function(x) quantile(x, .1)
q9 <- function(x) quantile(x, .9)


test_NA <- function(x) {
  is.null(x) || is.na(x)
}


draw_table <- function(dtf, draw_immediately = T) {
  p <-
    gridExtra::tableGrob(
      dtf,
      rows = NULL,
      theme = gridExtra::ttheme_default(
        # base_family = 'Arial_MT',
        base_size = 8,
        # core = list(fg_params = list(parse=T), hjust = 1),
        padding = grid::unit(c(2, 2), 'mm'),
        colhead = list(fg_params = list(parse=F))
      )
    )
  if (draw_immediately) {
    grid.draw(p)
  } else {
    return(p)
  }
}


assay2dtf <- function(so, assay = 'GS_MM',
  s_vars = c('condition_name', 'sample_origin', 'ifn_conc',
    'sn_dilution', 'percent.mt', 'tnf_conc', 'duration')) {
  s_vars <- intersect(s_vars, colnames(so@meta.data))
  dtf <-
    cbind(t(as.matrix(so[[assay]][,])), FetchData(so, s_vars)) %>%
    as.data.frame() %>%
    order_stim_group() %>%
    order_condition_name() %>%
    order_duration()
  return(dtf)
}


invert_names <- function(v) {
  setNames(names(v), v)
}


line_split <- function(v) {
  stringr::str_split(v, '\n') %>% 
    unlist() %>%
    { .[. != ''] } %>%
    { .[!stringr::str_detect(., '^#')] } %>%
    stringr::str_replace('^ *', '') %>%
    { .[. != ''] } %>%
    { . }
}


#' Take in a named vector with objects as names (e.g. genes) and
#' cluster indices as values, like the output of cutree
#'
#' @param ct_obj Cutted tree object
clust2list <- function(ct_obj) {
  maartenutils::auto_name(unique(ct_obj)) %>%
    map(~names(ct_obj)[ct_obj == .x])
}

 
get_lMs <- function(M, genes = detected_genes(M)) {
  lM <- subset_feats(M, genes)
  lMs <- t(scale(t(lM)))
  return(lMs)
}


format_6743 <- function(...) UseMethod('format_6743')
format_6743.data.frame <- function(dtf, rev_labels = F) {
  dtf$stim_group <- 
    factor(
      dtf$stim_group, 
      levels = 
        c('Mix + T', 'Ag-GAS + T', 'Mix PBS') %>%
        { if (rev_labels) rev(.) else . }
    )
  return(dtf)
}
format_6743.Seurat <- function(so, rev_labels = F) {
  so@meta.data <- format_6743(so@meta.data, rev_labels = rev_labels)
  return(so)
}


robust_scale <- function(M, drop_NA = TRUE) {
  dn <- dimnames(M)
  cM <- apply(M, 2, median)
  M <- t(apply(M, 1, function(x) x - cM))
  out <- M %*% diag(1/apply(M, 2, IQR))
  if (drop_NA) {
    idxs <- which(apply(out, 2, function(x) all(is.finite(x))))
    dn[[2]] <- dn[[2]][idxs]
    out <- out[, idxs]
  }
  dimnames(out) <- dn
  return(out)
}
robust_scale(cbind(rnorm(10, 10, 1), rnorm(10, 1000, 1)))


#' Compute the cosine matriX of the columns of two matrices
#'
#'
matrix_cosine <- function(X, Y) {
  stopifnot(nrow(X) == nrow(Y))
  co <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
  dimnames(co) = list(colnames(X), colnames(Y))
  for (i in 1:ncol(X)) {
    for (j in 1:ncol(Y)) {
      co[i,j] <- crossprod(X[, i], Y[, j]) /
        sqrt(crossprod(X[, i]) * crossprod(Y[, j]))
    }
  }
  return(co)
}


format_gs_names <- function(dtf) {
  colnames(dtf) <-
    stringr::str_replace(colnames(dtf),
      '.vanilla.unweighted.sum', '')
  colnames(dtf) <-
    stringr::str_replace_all(colnames(dtf),
      '-', '_')
  return(dtf)
}


call_func <- function(f, args, verbose = F) {
  if (missing(args) || is.null(args) || !is.list(args)) 
    stop('No valid args')
  overlapping_args <- intersect(names(formals(f)), names(args))
  # ignored_args <- setdiff(names(formals(f)), names(args))
  ignored_args <- setdiff(names(args), names(formals(f)))
  if (verbose && length(ignored_args) > 0) {
    message('Following list items from call to', ' are ignored:',
      paste0(ignored_args, collapse = ', '))
  }
  purrr::exec(f, !!!args[overlapping_args])
}



weighted_t_ci <- function(x, weights, conf.level = 0.95) {
  require(Hmisc)
  nx <- length(x)
  df <- nx - 1
  weights[is.na(weights)] <- 0
  vx <- wtd.var(x, weights, normwt = TRUE) ## From Hmisc
  mx <- weighted.mean(x, weights)
  stderr <- sqrt(vx/nx)
  tstat <- mx/stderr ## not mx - mu
  alpha <- 1 - conf.level
  cint <- qt(1 - alpha/2, df)
  cint <- tstat + c(-cint, cint)
  cint * stderr
}


append_name <- function(v, n) {
  setNames(v, paste0(names(v), n))
}


M_bool_to_int <- function(M) {
  M_bin <- matrix(0, nrow = nrow(M), ncol = ncol(M))
  M_bin[M] <- 1
  dimnames(M_bin) <- dimnames(M)
  return(M_bin)
}

M_bool_to_chr <- function(M, labels = c('Included', 'Excluded')) {
  M_bin <- matrix(labels[2], nrow = nrow(M), ncol = ncol(M))
  M_bin[M] <- labels[1]
  dimnames(M_bin) <- dimnames(M)
  return(M_bin)
}


compile_all_pairs <- function(v = c('a', 'b', 'c')) {
  L <- length(v)
  tidyr::expand_grid(cn1 = 1:L, cn2 = 1:L) %>%
    dplyr::filter(cn2 > cn1) %>%
    dplyr::mutate(v1 = v[cn1]) %>%
    dplyr::mutate(v2 = v[cn2]) %>%
    dplyr::select(v1, v2) %>%
    purrr::pmap(function(v1, v2) {
      c(v1, v2)
    }) %>%
    { . }
}
