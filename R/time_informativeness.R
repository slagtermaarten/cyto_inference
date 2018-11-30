#' Determine whether a gene is informative with respect to time for at
#' least one stimulus
#'
#'
determine_time_informativeness <- function(
  gene = 'CD274',
  unconvincing_lfc = 0.15,
  epsilon = .1375035,
  lookup_data = targets::tar_read('kallisto_5029')) {

  if (gene %nin% rownames(lookup_data))
    return(NULL)

  e_stats <- prepare_max_diff_stats(
    gene = gene,
    min_lfc = min_lfc,
    min_exp_diff = min_exp_diff,
    epsilon = epsilon,
    lookup_data = lookup_data
  )

  stats <- e_stats$stats
  p_dat <- e_stats$p_dat
  winning_stats <- e_stats$winning_stats

  if (T) {
    ts <- stats %>%
      dplyr::filter(stim_group == winning_stats$stim_group) %>%
      dplyr::arrange(duration) %>%
      dplyr::pull(lfc_from_us) %>%
      { . }
  } else {
    ts <- p_dat %>%
      dplyr::filter(stim_group == winning_stats$stim_group) %>%
      dplyr::arrange(duration) %>%
      pull(gexp) %>%
      { log2(. + 1) }
  }

  ## lfc between successive timepoints
  seq_diff_ts <- set_names(sort(diff(ts)),
    paste0('diff_', c('lo', 'mid', 'hi')))

  ## Log fold change after timepoint
  lfc_a_tp <- diff(ts)
  subs_change <- as.integer(abs(lfc_a_tp) > unconvincing_lfc)
  if (F && any(subs_change)) {
    # unconvincing_hops <- which(abs(diff(ts)) <= unconvincing_lfc)
    # convincing_hops <- setdiff(1:(length(ts)-1), unconvincing_hops)
    # change_rle <- rle(as.integer(subs_change))
    # group_idxs <- map(seq_along(change_rle$lengths), function(i) {
    #   if (i == 1) return(1:change_rle$lengths[i])
    #   (sum(change_rle$lengths[1:(i-1)])+1):sum(change_rle$lengths[1:i])
    # })

    ## Identify time regions for which this gene is informative
    subs_change_idxs <- which(subs_change == 1)
    ## Idxs refers to timepoints (i.e. left and right demarcations of
    ## intervals)
    ## Start points of informative intervals are identical to idx of
    ## interval
    start_idxs <- unique(c(1, subs_change_idxs))
    end_idxs <- start_idxs + 1
    ss_rle <- rle(sign(lfc_a_tp))
    ## Are there any consequent intervals with similar lfcs? Then
    ## combine them
    if (any(ss_rle$lengths > 1)) {
      ## i iterates over sign streaks
      for (i in seq_along(ss_rle$values)) {
        if (ss_rle$lengths[i] == 1) next()
        ## Any start point within this sign streak will be 'ingested'
        ## by the first start point of this sign streak, so first
        ## infer the coordinates of this sign streak
        start_ss <- if (i == 1) 1 else sum(ss_rle$lengths[1:(i-1)])+1
        end_ss <- sum(ss_rle$lengths[1:(i)])+1
        if (end_ss - start_ss <= 1) browser()
        ## These are start sites to be 'ingested' by this streak
        ingest_ss <- seq(start_ss+1, end_ss)
        start_idxs <- setdiff(start_idxs, ingest_ss)
        ## These are end sites to be 'ingested' by this streak
        ingest_ss <- seq(start_ss, end_ss-1)
        end_idxs <- setdiff(end_idxs, ingest_ss)
      }
    }
    stopifnot(length(end_idxs) == length(start_idxs))

    ## We now have a set of non-redundant, informative intervals,
    ## what's the fc over each of these intervals?
    int_lfc <- map_dbl(seq_along(start_idxs), ~ts[end_idxs[.x]] - ts[start_idxs[.x]])

    # group_idxs <- map(seq_along(start_idxs),
    #   ~start_idxs[.x]:end_idxs[.x])
    # ts
    # group_med <- map_dbl(group_idxs, ~median(ts[.x]))
    # group_lfc <- diff(group_med)
    # # 2^group_lfc
    # ## Timepoints after which 'group' changes occur
    # fc_tps <- setdiff(map_int(group_idxs, ~max(.x)), 4)

    # grp_lfc <- tibble(
    #   stp = tpi_to_tp(fc_tps),
    #   etp = tpi_to_tp(fc_tps + 1),
    #   lfc = group_lfc,
    #   # sign_label = sign(lfc),
    #   ## Map -1 to 1 and 1 to 2
    #   sign_label = glue('{c(\'-\',\'+\')[1.5+sign(lfc)*.5]}'),
    #   label = glue('{stp}h-{etp}h{sign_label}'),
    #   convincing_lfc = abs(lfc) >= unconvincing_lfc
    # ) %>%
    # dplyr::select(-sign_label) %>%
    # dplyr::mutate(N_groupings = n()+1)

    grp_lfc <- tibble::tibble(
      int_s = tpi_to_tp(start_idxs),
      int_e = tpi_to_tp(end_idxs),
      lfc = int_lfc,
      sign_label = glue('{c(\'-\',\'+\')[1.5+sign(int_lfc)*.5]}'),
      label = glue('{int_s}h-{int_e}h{sign_label}'),
      convincing_lfc = abs(lfc) >= unconvincing_lfc
    ) %>%
    dplyr::select(-sign_label) %>%
    dplyr::mutate(N_intervals = n())

    if (nrow(grp_lfc) > 1) {
      grp_lfc <- grp_lfc %>%
        # dplyr::filter(convincing_lfc == T) %>%
        dplyr::filter(abs(lfc) == max(abs(lfc)))
      # if (nrow(grp_lfc) != 1) browser()
    }
  } else {
    # grp_lfc <- tibble(
    #   stp = NA, etp = NA, lfc = NA, label = NA, convincing_lfc = F,
    #   N_groupings = NA)
    grp_lfc <- tibble::tibble(convincing_lfc = F)
  }

  ## lfc between all timepoints
  diff_ts <-
    tidyr::expand_grid(t1 = seq_along(ts), t2 = seq_along(ts)) %>%
    dplyr::filter(t1 < t2) %>%
    dplyr::mutate(diff_ts = ts[t2] - ts[t1]) %>%
    dplyr::mutate(len = t2 - t1) %>%
    { . }

  for (i in 1:nrow(diff_ts)) {
    all_i_tps <- diff_ts[i, ][['t1']]:diff_ts[i, ][['t2']]
    ## Determine whether the change between the indicated timepoints
    ## is monotonic
    diff_ts[i, 'mono'] <- 
      map_dbl(all_i_tps[-1], ~sign(ts[.x] - ts[.x-1])) %>%
      data.table::uniqueN() %>%
      { . == 1 }
  }

  max_t <- which.max(ts)
  if (F) {
    while (max_t != 1) {
      if ((ts[max_t] - ts[(max_t-1)]) <= unconvincing_lfc) {
        max_t <- max_t - 1
      } else {
        break
      }
    }

    max_duration <-
      p_dat %>%
        dplyr::filter(stim_group == winning_stats$stim_group) %>%
        dplyr::pull(duration) %>%
        pluck(max_t)
    max_duration <- levels(max_duration)[max_duration] %>%
      as.numeric()

    out <- list(
      'gene' = gene,
      'max_duration' = max_duration)
    names(out)[2] <- glue('max_duration_{unconvincing_lfc}')
  } else {
    time_class <-
      if (all(abs(diff_ts$diff_ts) < unconvincing_lfc)) {
        'no_time_informativeness'
      } else {
        if (all((ts[max_t] - ts[-max_t]) > unconvincing_lfc)) {
          'dominant_tp'
        } else if (any((ts[max_t] - ts[-max_t]) > unconvincing_lfc)) {
          'non_dominant_informative_tp'
        } else if (any(grp_lfc$convincing_lfc)) {
          grp_lfc$label
        } else {
          'no_dominant_tp'
        }
      }

    out <- tibble(
        'gene' = gene,
        'max_timepoint' = as.integer(levels(p_dat$duration)[max_t]),
        'max_stimgroup' = winning_stats$stim_group,
        'time_class' = time_class
      ) %>%
      dplyr::bind_cols(as_tibble(t(seq_diff_ts))) %>%
      dplyr::bind_cols(grp_lfc)
  }

  return(out)
}
