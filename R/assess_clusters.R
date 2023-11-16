pielou_evenness <- function(v) {
  v <- v / sum(v, na.rm = T)
  H <- -sum(v * log(v), na.rm = T)
  H / log(length(v))
}
# pielou_evenness(c(.1, .1, .1, .1)) == 1
# pielou_evenness(c(.1, .1, .1, .1, .2)) < 1


compute_cluster_evenness <- function(
  so, search_stim_group = 'Exposed to T-cells in vivo') {
  cluster_props <- so@meta.data %>% 
    group_by(seurat_clusters) %>%
    dplyr::summarise(
      p = sum(stim_group == search_stim_group, na.rm = T) / n())
  pielou_evenness(cluster_props$p)
}


#' Compare two sets of objects
#'
#'
set_analysis <- function(a, b) {
  list('union' = union(a, b), 'intersection' = intersect(a, b),
       'unique_a' = setdiff(a, b), 'unique_b' = setdiff(b, a))
}


#' Compute AUC of gene vs. unstimulated for each of the stim groups
#' (:= combination of stimulus and duration), in order to decide
#' which stimulus most strongly affects the gene and to what degree
#'
#'
compute_gene_AUC <- result_cacher(
  f = function(goi = 'CD274', lookup_data = load_EE_feature_mat()) {
    if (goi %nin% rownames(lookup_data))
      return(NULL)

    p_dat <- lookup_data[which(rownames(lookup_data) == goi), ] %>%
      named_vec_to_dt('sample_name', 'gexp') %>%
      merge_sample_annotation(exp = '5892') %>%
      { . }

    if (!grepl('^PC\\d+$', goi)) {
      p_dat <- dplyr::mutate(p_dat, gexp = log2(gexp + 1))
    }

    max_conc_stim_groups <- c('Unstimulated',
      '10 ng/ml TNFa', '100 ng/ml IFNy',
      '100 ng/ml IFNy 10 ng/ml TNFa')

    ## Find time point where biggest difference occurs for this gene
    stim_group_AUC_data <-
      p_dat %>%
      dplyr::group_by(duration) %>%
      dplyr::summarize(
        stim_group, gexp,
        diff_vs_us = gexp - gexp[stim_group == 'Unstimulated'],
        .groups = 'drop') %>%
      # debug_pipe %>%
      # dplyr::filter(stim_group != 'Unstimulated') %>%
      # dplyr::filter(!grepl('SN', stim_group)) %>%
      dplyr::filter(stim_group %in% max_conc_stim_groups) %>%
      # debug_pipe %>%
      dplyr::arrange(stim_group, duration) %>%
      dplyr::mutate(
        duration = as.integer(as.character(duration))) %>%
      dplyr::group_by(stim_group) %>%
      dplyr::summarize(
        auc = DescTools::AUC(x = duration, y = gexp),
        .groups = 'drop') %>%
      dplyr::mutate(
        auc_diff = auc - auc[stim_group == 'Unstimulated']) %>%
      dplyr::mutate(auc_norm = auc_diff / sum(auc_diff)) %>%
      # dplyr::mutate(auc_norm = 2^(exp_mod*auc)) %>%
      # dplyr::mutate(auc_norm = auc_norm / sum(auc_norm)) %>%
      # dplyr::mutate(stim_group[!grepl('IFNy.*TNFa$', stim_group)])
      { . }

    if (F) {
      scores_wo_combo <- stim_group_AUC_data %>%
        dplyr::filter(!grepl('IFNy.*TNFa$', stim_group)) %>%
        dplyr::mutate(auc_norm_woc = auc_diff / sum(auc)) %>%
        # dplyr::mutate(auc_norm_woc = 2^(exp_mod*auc)) %>%
        # dplyr::mutate(auc_norm_woc = auc_norm_woc / sum(auc_norm_woc)) %>%
        dplyr::select(stim_group, auc_norm_woc) %>%
        { . }
      stim_group_AUC_data <- stim_group_AUC_data %>%
        dplyr::full_join(scores_wo_combo, by = c('stim_group'))
    }

    ## AUC works as expected
    # DescTools::AUC(x = c(0, 5), y = c(0, 3))
    # DescTools::AUC(x = c(0, 5), y = c(-1, -3))
    return(stim_group_AUC_data)
  },
  filename = function() {
    file.path(rds_dir, 'compute_gene_AUC', glue::glue('{goi}.rds'))
  },
  min_mod_time = '2021-06-21 09:52'
)


classify_gene_AUC <- function(goi, redo = F,
  lookup_data = load_EE_feature_mat()) {

  dtf <- compute_gene_AUC(goi = goi, redo = redo,
    lookup_data = lookup_data)

  res <- dtf %>%
    dplyr::select(stim_group, auc, auc_diff) %>%
    dplyr::summarize(
      tnf_auc = auc[stim_group == '10 ng/ml TNFa'],
      ifn_auc = auc[stim_group == '100 ng/ml IFNy'],
      combo_auc = auc[stim_group == '100 ng/ml IFNy 10 ng/ml TNFa'],
      norm_tnf_auc = auc_diff[stim_group == '10 ng/ml TNFa'],
      norm_ifn_auc = auc_diff[stim_group == '100 ng/ml IFNy'],
      norm_combo_auc = auc_diff[stim_group == '100 ng/ml IFNy 10 ng/ml TNFa'],
      ## (Observed - expected) / expected
      synergy_score =
        (auc[stim_group == '100 ng/ml IFNy 10 ng/ml TNFa'] -
         auc[stim_group == '100 ng/ml IFNy'] -
         auc[stim_group == '10 ng/ml TNFa']) /
        (auc[stim_group == '100 ng/ml IFNy'] +
         auc[stim_group == '10 ng/ml TNFa'])
    )
  if (is.na(res$synergy_score)) res$synergy_score <- 0

  if (res$synergy_score >= 0.2) {
    res$class <- 'synergy'
  } else if (
    abs(res$norm_tnf_auc) > abs(res$norm_ifn_auc) &&
    abs(res$norm_tnf_auc) >= 10) {
    res$class <- 'tnf'
  } else if (
    abs(res$norm_ifn_auc) > abs(res$norm_tnf_auc) &&
    abs(res$norm_ifn_auc) >= 10) {
    res$class <- 'ifn'
  } else {
    res$class <- 'none'
  }

  res$gene <- goi

  return(res)
}


#' Determine what cytokine achieves the largest difference with the
#' unstimulated setting, across time points
#' This somewhat works but is far from perfect.
#'
#'
classify_gene_using_bulk <- result_cacher(
  f = function(goi = 'CD274',
               lookup_data = load_EE_feature_mat(),
               diff_vs_winning_threshold = NULL) {
    if (goi %nin% rownames(lookup_data))
      return(NULL)

    p_dat <- lookup_data[which(rownames(lookup_data) == goi), ] %>%
      named_vec_to_dt('sample_name', 'gexp') %>%
      merge_sample_annotation(exp = '5892')

    if (!grepl('^PC\\d+$', goi)) {
      p_dat <- dplyr::mutate(p_dat, gexp = log2(gexp + 1))
    }

    ## Find time point where biggest difference occurs for this gene
    max_diff_timepoint <- p_dat %>%
      dplyr::group_by(duration) %>%
      dplyr::summarize(max_diff = max(gexp, na.rm = T) -
        min(gexp, na.rm = T), .groups = 'drop') %>%
      arrange(-abs(max_diff)) %>%
      pull(duration) %>% pluck(1)

    # if (strict)
    ## Subset to relevant data
    p_dat <- p_dat %>%
      dplyr::filter(duration == max_diff_timepoint) %>%
      dplyr::filter(sn_rank == 1)

    ## Extract unstimulated gene expression (potential time effect
    ## only)
    gexp_us <- p_dat %>%
      dplyr::filter(grepl('unstim', sample_name)) %>%
      pull(gexp)

    # winning_stim <- p_dat %>%
    #   mutate(exp_diff = abs(gexp - gexp_us)) %>%
    #   arrange(-exp_diff) %>% pull(stim_group) %>% pluck(1)

    FCs <- p_dat %>%
      dplyr::mutate(exp_diff = abs(gexp - gexp_us)) %>%
      dplyr::arrange(-exp_diff)

    winning_stim_ranks <- FCs %>%
      dplyr::top_n(1) %>%
      dplyr::select(tnf_rank, ifn_rank) %>% unlist

    winning_fc <- FCs %>%
      dplyr::filter(tnf_rank == winning_stim_ranks[1] &
                    ifn_rank == winning_stim_ranks[2]) %>%
      pull(exp_diff)

    if (!is.null(diff_vs_winning_threshold)) {
      diff_vs_winning <- FCs %>%
        dplyr::filter(!(tnf_rank > winning_stim_ranks[1] &
                        ifn_rank > winning_stim_ranks[2])) %>%
        dplyr::mutate(exp_diff_vs_winning = winning_fc / exp_diff) %>%
        dplyr::filter(exp_diff_vs_winning <= diff_vs_winning_threshold)
    }

    if (winning_stim_ranks[1] == 4 && winning_stim_ranks[2] == 4) {
      label <- 'comb'
      ## Test for synergy
      ifn_mono <-
        p_dat %>% filter(tnf_rank == 1 & ifn_rank == winning_stim_ranks[2]) %>%
        pull(gexp)
      tnf_mono <-
        p_dat %>% filter(ifn_rank == 1 & tnf_rank == winning_stim_ranks[1]) %>%
        pull(gexp)
      comb_exp <- p_dat %>% filter(tnf_rank == winning_stim_ranks[1] &
                                   ifn_rank == winning_stim_ranks[2]) %>%
        pull(gexp)
      if ((1.10 * (ifn_mono + tnf_mono)) < comb_exp) {
        label <- 'synergy'
      }
    } else if (winning_stim_ranks[1] > winning_stim_ranks[2]) {
      label <- 'TNFa'
    } else if (winning_stim_ranks[1] < winning_stim_ranks[2]) {
      label <- 'IFNy'
    } else {
      label <- NA_character_
    }
    return(label)
  },
  filename = function() {
    file.path(rds_dir, 'gene_classification_using_bulk',
      glue::glue('{goi}{make_flag(diff_vs_winning_threshold)}.rds'))
  },
  min_mod_time = '2021-01-06 14:23'
)
# classify_gene_using_bulk(goi = 'CD274')


#' Find genes that are diff exp between any two clusters
#'
#'
identify_contrasting_genes <-
  function(seurat_object, features = VariableFeatures(seurat_object)) {
  clusters <- 
    setdiff(levels(seurat_object@meta.data$seurat_clusters), 0) %>%
    as.numeric
  all_comps <- expand.grid(ident.1 = clusters, ident.2 = clusters) %>%
    dplyr::filter(ident.2 > ident.1) %>%
    dplyr::arrange(ident.1, ident.2) %>%
    plyr::adply(1, function(r) {
      message('Cluster ',  r[[1]], ' vs. ', r[[2]])
      markers <- tryCatch(FindMarkers(pbmc,
                                      ident.1 = r[['ident.1']],
                                      ident.2 = r[['ident.2']],
                                      features = features),
                          error = function(e) { NULL })
      if (null_dat(markers)) {
        return(NULL)
      } else {
        markers %>%
          tibble::rownames_to_column('gene') %>%
          dplyr::mutate(
            ident.1 = r[['ident.1']], 
            ident.2 = r[['ident.2']]) %>%
          dplyr::filter(p_val_adj <= .25) %>%
          return
      }
  })

  ## Order results identically as in queried features
  all_comps <- all_comps[order(match(all_comps$gene, features)), ]
}


# read_genelist <- function(name = '5092_5310_shared_informativeV3') {
#   fn <- file.path(data_raw_dir, 'gene_lists', 
#     glue::glue('{name}_genes.txt'))
#   genes <- readLines(fn)
#   genes <- genes[!grepl('^#', genes)]
#   genes <- gsub('\\s*#.*$', '', genes)
#   return(genes)
# }
# read_genelist(name = '5092_5310_shared_informativeV3')
# read_genelist(name = 'weird')
